"""Generiere Mesh mit besser positioniertem Zylinder für Vortex Shedding.

Problem mit Original-Mesh:
- Zylinder bei x=0.2 (zu nah am Einlass)
- Strömung nicht voll entwickelt
- Instabilität wird unterdrückt

Lösung:
- Zylinder bei x=1.0 (Mitte des Kanals)
- Längerer Einlaufbereich für Strömungsentwicklung
"""
from dolfin import *
import mshr
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def generate_improved_cylinder_mesh(output_path="meshes/cylinder_improved.xml", resolution=64):
    """Generiere verbessertes Mesh für Vortex Shedding.
    
    Geometrie (angepasst für bessere Instabilität):
    - Kanal: Länge L=3.0, Höhe H=0.41
    - Zylinder: Mittelpunkt (1.0, 0.205), Radius r=0.05
    
    Args:
        output_path: Pfad für XML-Mesh
        resolution: Auflösung (höher = feiner)
    """
    print("=" * 60)
    print("VERBESSERTES MESH FÜR VORTEX SHEDDING")
    print("=" * 60)
    
    # Geometrie-Parameter (optimiert)
    L = 3.0   # Längerer Kanal für Nachlauf
    H = 0.41  # Standardhöhe (DFG Benchmark)
    c_x = 1.0 # Zylinder bei 1/3 (mehr Einlauf)
    c_y = H/2 # Zylinder in Kanalmitte
    r = 0.05  # Zylinderradius
    
    print(f"\nGeometrie:")
    print(f"  Kanal: L={L}m, H={H}m")
    print(f"  Zylinder: Position ({c_x}, {c_y}), Radius r={r}")
    print(f"  Einlaufstrecke: {c_x}m ({c_x/L*100:.1f}% der Gesamtlänge)")
    print(f"  Nachlaufstrecke: {L-c_x}m")
    print(f"  Durchmesser/Höhe: D/H = {2*r/H:.1%}")
    
    # Geometrie definieren mit mshr
    channel = mshr.Rectangle(Point(0, 0), Point(L, H))
    cylinder = mshr.Circle(Point(c_x, c_y), r)
    
    # Kanal minus Zylinder
    domain = channel - cylinder
    
    # Mesh generieren
    print(f"\nGeneriere Mesh (Auflösung={resolution})...")
    mesh = mshr.generate_mesh(domain, resolution)
    
    # Speichern
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    File(output_path) << mesh
    
    print(f"\n✓ Mesh gespeichert: {output_path}")
    print(f"  Vertices: {mesh.num_vertices()}")
    print(f"  Cells: {mesh.num_cells()}")
    print("\n" + "=" * 60)
    
    return mesh


if __name__ == "__main__":
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    mesh_file = os.path.join(project_root, "meshes", "cylinder_improved.xml")
    
    mesh = generate_improved_cylinder_mesh(mesh_file, resolution=64)
    print(f"\nMesh erfolgreich generiert!")
    print(f"Verwende in Simulation mit:")
    print(f"  run_cylinder_simulation(mesh_path='{mesh_file}', ...)")
