from dolfin import *
from src.boundary_conditions import apply_bcs_cylinder_flow
import numpy as np


class NavierStokesSolver:
    def __init__(self, mesh_file, params):
        """Initialisiere Navier-Stokes Löser

        Implementiert: Stokes (linear, `params.linearized_only=True`) oder
        volle Navier–Stokes mit Picard-Iteration (nichtlinear, default).
        Randbedingungen: Einlass (Dirichlet), Wände (no-slip), Auslass (Druck-Dirichlet).
        """
        # Accept either a path to a mesh file or a dolfin.Mesh instance
        try:
            MeshClass = Mesh
        except NameError:
            MeshClass = None

        if MeshClass is not None and isinstance(mesh_file, MeshClass):
            self.mesh = mesh_file
        else:
            self.mesh = Mesh(mesh_file)
        self.params = params
        self.W = None
        self.bcs = []

    def define_function_space(self):
        """Definiere gemischte Funktionenräume für Geschwindigkeit und Druck"""
        # Taylor-Hood Elemente (P2-P1) mit MixedElement
        Ve = VectorElement("P", self.mesh.ufl_cell(), 2)
        Qe = FiniteElement("P", self.mesh.ufl_cell(), 1)
        W_elem = MixedElement([Ve, Qe])
        self.W = FunctionSpace(self.mesh, W_elem)

    def define_boundary_conditions(self):
        """Definiere Randbedingungen für Zylinderumströmung."""
        if self.W is None:
            raise RuntimeError("Function space not defined. Call define_function_space() first.")
        
        # Verwende die zentrale BC-Definition
        self.bcs = apply_bcs_cylinder_flow(
            self.W, 
            self.mesh, 
            self.params.inlet_velocity, 
            self.params.outlet_pressure
        )
        
        # Füge Zylinder-BC hinzu (no-slip am Zylinder)
        coords = self.mesh.coordinates()
        xmin = float(np.min(coords[:, 0]))
        eps = 1e-8 * max(1.0, abs(np.max(coords[:, 0]) - xmin))
        
        class Cylinder(SubDomain):
            def inside(self, x, on_boundary):
                cx, cy, r = 0.5, 0.2, 0.05
                return on_boundary and ((x[0] - cx)**2 + (x[1] - cy)**2 < (r * 1.1)**2)
        
        self.bcs.append(DirichletBC(self.W.sub(0), Constant((0.0, 0.0)), Cylinder()))

    def _setup_problem(self):
        """Setup Funktionsraum und Randbedingungen falls nötig."""
        if self.W is None:
            self.define_function_space()
        if not self.bcs:
            self.define_boundary_conditions()

    def _picard_solve(self, w, u_prev, nu, f, dt=None, u_n=None):
        """Führe Picard-Iteration durch.
        
        Returns: (converged, iterations)
        """
        (u, p) = TrialFunctions(self.W)
        (v, q) = TestFunctions(self.W)
        
        for i in range(self.params.picard_maxiter):
            a = (nu * inner(grad(u), grad(v)) * dx +
                 inner(dot(grad(u), u_prev), v) * dx -
                 div(v) * p * dx - q * div(u) * dx)
            L = inner(f, v) * dx
            
            if dt is not None:  # Transient
                a += (1.0/dt) * inner(u, v) * dx
                L += (1.0/dt) * inner(u_n, v) * dx
            
            w_old = Function(self.W)
            w_old.assign(w)
            solve(a == L, w, self.bcs, solver_parameters={"linear_solver": "lu"})
            
            u_prev.assign(w.split(deepcopy=True)[0])
            
            dw = Function(self.W)
            dw.assign(w)
            dw.vector()[:] -= w_old.vector()
            res = norm(dw, "l2")
            
            if self.params.verbose and dt is None:
                print(f"  Picard iter {i + 1}: residuum = {res:.6e}")
            
            if res < self.params.tol:
                return True, i + 1
        
        return False, self.params.picard_maxiter

    def solve_steady(self, save_path=None):
        """Löse stationäre Navier–Stokes."""
        self._setup_problem()
        
        nu = Constant(self.params.nu)
        f = Constant((0.0, 0.0))
        w = Function(self.W)
        u_prev = Function(self.W.sub(0).collapse())
        
        if self.params.linearized_only:
            # Stokes
            (u, p) = TrialFunctions(self.W)
            (v, q) = TestFunctions(self.W)
            a = nu * inner(grad(u), grad(v)) * dx - div(v) * p * dx - q * div(u) * dx
            L = inner(f, v) * dx
            solve(a == L, w, self.bcs, solver_parameters={"linear_solver": "lu"})
            if self.params.verbose:
                print(f"Stokes: ||u||_L2 = {norm(w.split()[0], 'L2'):.6e}")
        else:
            # Navier-Stokes
            if self.params.verbose:
                print("Starte Picard-Iteration für Navier–Stokes...")
            conv, iters = self._picard_solve(w, u_prev, nu, f)
            if self.params.verbose:
                print(f"{'Konvergiert' if conv else 'Max. Iterationen'} nach {iters} Iterationen")
                print(f"Navier–Stokes: ||u||_L2 = {norm(w.split()[0], 'L2'):.6e}")
        
        if save_path:
            self._save_solution(w, save_path)
        return w

    def _save_solution(self, w, save_path, t=0.0):
        """Speichere Lösung als XDMF für Post-Processing."""
        try:
            file = XDMFFile(save_path)
            file.write(w.sub(0), t)
            file.write(w.sub(1), t)
            if self.params.verbose:
                print(f"Lösung gespeichert: {save_path} (t={t:.4f})")
        except Exception as e:
            if self.params.verbose:
                print(f"Warnung: Konnte Lösung nicht speichern: {e}")

    def solve(self, save_path=None):
        """Wrapper: ruft solve_steady() oder solve_transient() je nach params.steady."""
        if self.params.steady:
            return self.solve_steady(save_path)
        else:
            return self.solve_transient(save_path)

    def solve_transient(self, save_path=None):
        """Löse zeitabhängige Navier–Stokes mit implizitem Euler + Picard."""
        self._setup_problem()
        
        nu, dt, f = Constant(self.params.nu), Constant(self.params.dt), Constant((0.0, 0.0))
        w_n, w = Function(self.W), Function(self.W)
        u_prev = Function(self.W.sub(0).collapse())
        
        # Anfangslösung
        if self.params.verbose:
            print("Berechne Anfangslösung (Stokes)...")
        w_n.assign(self.solve_steady(save_path=None))
        
        # Störung
        if hasattr(self.params, 'add_perturbation') and self.params.add_perturbation:
            self._add_perturbation(w_n)
        
        # Output-Setup
        xdmf_file, u_out, p_out = None, None, None
        if save_path:
            xdmf_file, u_out, p_out = self._setup_output(save_path)
        
        if self.params.verbose:
            print(f"\n=== Transiente Simulation: t=[0, {self.params.t_end}], dt={self.params.dt} ===")
        
        # Zeit-Loop
        t, timesteps = 0.0, []
        while t < self.params.t_end:
            t += self.params.dt
            u_prev.assign(w_n.split(deepcopy=True)[0])
            
            _, iters = self._picard_solve(w, u_prev, nu, f, dt, w_n.split(deepcopy=True)[0])
            
            w_n.assign(w)
            w_copy = Function(self.W)
            w_copy.assign(w)
            timesteps.append((t, w_copy))
            
            if save_path:
                u_sol, p_sol = w.split(deepcopy=True)
                assign(u_out, u_sol)
                assign(p_out, p_sol)
                xdmf_file.write(u_out, t)
                xdmf_file.write(p_out, t)
            
            if self.params.verbose and int(t / self.params.dt) % 10 == 0:
                print(f"  t={t:.4f}: ||u||_L2={norm(w.split()[0], 'L2'):.6e}, Picard: {iters} iter")
        
        if save_path:
            xdmf_file.close()
        
        if self.params.verbose:
            print(f"\n✓ Transiente Simulation: {len(timesteps)} Zeitschritte")
        
        return timesteps
    
    def _add_perturbation(self, w_n):
        """Füge Störung zur Anfangslösung hinzu."""
        if self.params.verbose:
            print("Füge asymmetrische Störung hinzu (Random Noise)...")
        V_vel = self.W.sub(0).collapse()
        u_pert = Function(V_vel)
        np.random.seed(42)
        u_pert.vector()[:] = np.random.randn(u_pert.vector().size()) * (0.005 * self.params.inlet_velocity)
        u_current = w_n.split(deepcopy=True)[0]
        u_current.vector()[:] += u_pert.vector()[:]
        assign(w_n.sub(0), u_current)
        if self.params.verbose:
            print(f"  Störung: ||u_pert||_L2 = {norm(u_pert, 'L2'):.4e}")
    
    def _setup_output(self, save_path):
        """Setup XDMF Output für transiente Simulation."""
        xdmf_file = XDMFFile(save_path)
        xdmf_file.parameters["flush_output"] = True
        xdmf_file.parameters["rewrite_function_mesh"] = False
        xdmf_file.parameters["functions_share_mesh"] = True
        
        V_out = FunctionSpace(self.mesh, VectorElement("P", self.mesh.ufl_cell(), 2))
        Q_out = FunctionSpace(self.mesh, FiniteElement("P", self.mesh.ufl_cell(), 1))
        u_out, p_out = Function(V_out), Function(Q_out)
        u_out.rename("velocity", "velocity")
        p_out.rename("pressure", "pressure")
        
        return xdmf_file, u_out, p_out


if __name__ == "__main__":
    print("Navier-Stokes (Stokes) Solver initialisiert")