"""Simulation: Zylinderumströmung (Flow around Cylinder).

Dieses Skript:
1. Lädt/generiert Mesh für Zylinder
2. Setzt Parameter (Reynolds-Zahl basiert)
3. Löst transiente Navier–Stokes
4. Berechnet Kräfte (Drag/Lift) am Zylinder
5. Visualisiert Ergebnisse
"""
import os
import sys
import numpy as np

# Headless backend for servers/CI
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dolfin import *

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.parameters import Params
from src.solver import NavierStokesSolver


def compute_forces(w, mesh, nu, rho=1.0, cylinder_center=(0.5, 0.2), cylinder_radius=0.05):
    """Berechne Drag- und Lift-Kräfte am Zylinder.
    
    Verwendet Cauchy-Spannungstensor: σ = -pI + ν*ρ*(∇u + ∇u^T)
    Kraft: F = ∫_∂Ω σ·n ds
    
    Args:
        w: Mixed solution (u, p)
        mesh: Domain mesh
        nu: kinematische Viskosität
        rho: Dichte
        cylinder_center, cylinder_radius: Geometrie
    """
    u, p = w.split()
    
    # Konstanten
    nu_const = Constant(nu)
    rho_const = Constant(rho)
    mu_const = rho_const * nu_const  # dynamische Viskosität
    
    # Normale und Spannungstensor
    n = FacetNormal(mesh)
    I = Identity(u.geometric_dimension())
    # σ = -pI + ρν(∇u + ∇u^T) [mit impliziter ρ in Viskosität]
    sigma = -p * I + mu_const * (grad(u) + grad(u).T)
    
    # Zylinder-Boundary markieren (SubDomain)
    class CylinderBoundary(SubDomain):
        def inside(self, x, on_boundary):
            cx, cy = cylinder_center
            r = cylinder_radius
            return on_boundary and ((x[0] - cx)**2 + (x[1] - cy)**2 < (r * 1.1)**2)
    
    cylinder_boundary = CylinderBoundary()
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)
    cylinder_boundary.mark(boundaries, 1)
    
    ds_cyl = ds(subdomain_data=boundaries)(1)
    
    # Kraft = ∫ σ·n ds (integriert über Zylinder-Oberfläche)
    F = dot(sigma, n)
    drag = assemble(F[0] * ds_cyl)   # x-Komponente (Widerstand)
    lift = assemble(F[1] * ds_cyl)   # y-Komponente (Auftrieb)
    
    return drag, lift


def run_cylinder_simulation(transient=True, T=5.0, dt=0.01, Re=100, perturb=True):
    """Hauptsimulation für Zylinderumströmung.
    
    Args:
        transient: Zeitabhängig (True) oder stationär (False)
        T: Endzeit
        dt: Zeitschritt
        Re: Reynolds-Zahl
        perturb: Initiale Störung hinzufügen (triggert Vortex Shedding)
    """
    print("=" * 60)
    print("ZYLINDERUMSTRÖMUNG SIMULATION")
    print("=" * 60)
    
    # Get project root
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    mesh_file = os.path.join(project_root, "meshes", "cylinder_flow.xml")
    if not os.path.exists(mesh_file) or os.path.getsize(mesh_file) < 100:
        print("\nMesh existiert nicht oder ist leer. Generiere...")
        sys.path.insert(0, project_root)
        from src.mesh_generator import generate_cylinder_mesh
        mesh = generate_cylinder_mesh(mesh_file, resolution=64)
    else:
        mesh = Mesh(mesh_file)
        print(f"\nMesh geladen: {mesh_file}")
        print(f"  Vertices: {mesh.num_vertices()}, Cells: {mesh.num_cells()}")
    
    # Parameter basierend auf Reynolds-Zahl
    # Re = U*D/nu, D=2r=0.1, U=inlet_velocity
    D = 2 * 0.05  # Zylinderdurchmesser
    U_mean = 1.0  # Mittlere Einströmgeschwindigkeit
    nu = U_mean * D / Re
    
    params = Params(
        inlet_velocity=U_mean,
        nu=nu,
        steady=not transient,
        dt=dt,
        t_end=T,
        tol=1e-3,  # Gelockerte Toleranz für schnellere Konvergenz
        picard_maxiter=10,  # Weniger Iterationen (schneller)
        verbose=True,
        add_perturbation=perturb
    )
    
    print(f"\nParameter:")
    print(f"  Reynolds-Zahl: Re = {Re}")
    print(f"  Viskosität: ν = {nu:.6f}")
    print(f"  Einlass-Geschwindigkeit: U = {U_mean}")
    print(f"  Modus: {'Transient' if transient else 'Stationär'}")
    if transient:
        print(f"  Zeitschritte: dt={dt}, T={T} ({int(T/dt)} steps)")
    
    # Solver initialisieren
    solver = NavierStokesSolver(mesh, params)
    
    if perturb and transient:
        print("\n⚠ Perturbation aktiviert: Verwende reduzierte Toleranz für numerisches Rauschen")
    
    # Lösen
    output_path = os.path.join(project_root, "results", "cylinder_flow.xdmf")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    if transient:
        print("\n" + "=" * 60)
        print("TRANSIENTE SIMULATION STARTET...")
        print("=" * 60)
        timesteps = solver.solve(save_path=output_path)
        
        # Kräfte berechnen
        times = []
        drags = []
        lifts = []
        
        print("\nBerechne Kräfte (Drag/Lift)...")
        for t, w in timesteps[::5]:  # Jeder 5. Schritt
            drag, lift = compute_forces(w, mesh, nu)
            times.append(t)
            drags.append(drag)
            lifts.append(lift)
        
        # Visualisierung
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        ax1.plot(times, drags, 'b-', linewidth=2, label='Drag (Widerstand)')
        ax1.set_xlabel('Zeit [s]')
        ax1.set_ylabel('Drag-Kraft [N]')
        ax1.set_title(f'Widerstandskraft am Zylinder (Re={Re})')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        ax2.plot(times, lifts, 'r-', linewidth=2, label='Lift (Auftrieb)')
        ax2.set_xlabel('Zeit [s]')
        ax2.set_ylabel('Lift-Kraft [N]')
        ax2.set_title(f'Auftriebskraft am Zylinder (Re={Re})')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        force_plot = os.path.join(project_root, "results", "cylinder_forces.png")
        plt.savefig(force_plot, dpi=150)
        print(f"\n✓ Kraft-Plot gespeichert: {force_plot}")
        
        # Statistik
        print(f"\nKraft-Statistik:")
        print(f"  Drag: mean={np.mean(drags):.6f}, max={np.max(drags):.6f}, min={np.min(drags):.6f}")
        print(f"  Lift: mean={np.mean(lifts):.6f}, max={np.max(lifts):.6f}, min={np.min(lifts):.6f}")
        
        # Letzte Lösung visualisieren
        _, w_final = timesteps[-1]
        u_final, p_final = w_final.split()
        
        print(f"\nFinal solution (t={times[-1]:.2f}):")
        print(f"  ||u||_L2 = {norm(u_final, 'L2'):.6e}")
        print(f"  ||p||_L2 = {norm(p_final, 'L2'):.6e}")
        
    else:
        print("\n" + "=" * 60)
        print("STATIONÄRE SIMULATION STARTET...")
        print("=" * 60)
        w = solver.solve(save_path=output_path)
        u, p = w.split()
        
        drag, lift = compute_forces(w, mesh, nu)
        print(f"\nKräfte (stationär):")
        print(f"  Drag = {drag:.6f}")
        print(f"  Lift = {lift:.6f}")
    
    print("\n" + "=" * 60)
    print(f"✓ SIMULATION ABGESCHLOSSEN")
    print(f"✓ Ergebnisse: {output_path}")
    print("=" * 60)


if __name__ == "__main__":
    # Hohe Reynolds-Zahl Simulation: Re=200, T=3.0s mit größerem Zeitschritt
    run_cylinder_simulation(transient=True, T=8.0, dt=0.01, Re=200, perturb=True)
    
    # Alternative: Kürzere Simulation für schnellere Ergebnisse
    # run_cylinder_simulation(transient=True, T=1.0, dt=0.01, Re=200, perturb=True)
    
    # Beispiel 1: Transient, Re=100 (erwartet: Periodische Wirbelablösung)
    # run_cylinder_simulation(transient=True, T=8.0, dt=0.01, Re=100)
    
    # Beispiel 2: Stationär, Re=20 (stabil, keine Wirbelablösung)
    # run_cylinder_simulation(transient=False, Re=20)
