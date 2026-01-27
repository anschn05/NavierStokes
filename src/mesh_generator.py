"""Mesh-Generator für Zylinderumströmung (Flow around Cylinder).

Erzeugt ein 2D-Rechteck-Kanal mit einem Zylinder-Hindernis.
Klassisches Benchmark-Problem (DFG-Benchmark 2D-1, 2D-2, 2D-3).
"""
from dolfin import *
import mshr
import os


def generate_cylinder_mesh(output_path="meshes/cylinder_flow.xml", resolution=64):
    """Generiere Mesh für Zylinderumströmung.
    
    Geometrie (modifiziert für Vortex Shedding):
    - Kanal: Länge L=2.2, Höhe H=0.41
    - Zylinder: Mittelpunkt (0.5, 0.2), Radius r=0.05 (weiter vom Einlass!)
    
    Args:
        output_path: Pfad für XML-Mesh
        resolution: Auflösung (höher = feiner, langsamer)
    """
    # Geometrie-Parameter (optimiert für Vortex Shedding)
    L = 2.2   # Kanallänge
    H = 0.41  # Kanalhöhe
    c_x = 0.5 # Zylinder x-Position (5D vom Einlass für Entwicklung)
    c_y = 0.2 # Zylinder y-Position (zentriert)
    r = 0.05  # Zylinderradius
    
    # Geometrie definieren mit mshr
    channel = mshr.Rectangle(Point(0, 0), Point(L, H))
    cylinder = mshr.Circle(Point(c_x, c_y), r)
    
    # Kanal minus Zylinder
    domain = channel - cylinder
    
    # Mesh generieren
    mesh = mshr.generate_mesh(domain, resolution)
    
    # Speichern
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    File(output_path) << mesh
    
    print(f"Mesh generiert: {output_path}")
    print(f"  Vertices: {mesh.num_vertices()}")
    print(f"  Cells: {mesh.num_cells()}")
    print(f"  Geometrie: L={L}, H={H}, Zylinder bei ({c_x},{c_y}), r={r}")
    
    return mesh


def generate_channel_mesh(output_path="meshes/channel.xml", nx=100, ny=20):
    """Generiere einfaches Kanal-Mesh (Rechteck ohne Hindernis).
    
    Args:
        output_path: Pfad für XML-Mesh
        nx, ny: Anzahl Elemente in x/y
    """
    L = 2.2
    H = 0.41
    
    mesh = RectangleMesh(Point(0, 0), Point(L, H), nx, ny)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    File(output_path) << mesh
    
    print(f"Kanal-Mesh generiert: {output_path}")
    print(f"  Vertices: {mesh.num_vertices()}")
    print(f"  Cells: {mesh.num_cells()}")
    
    return mesh


if __name__ == "__main__":
    print("=== Mesh-Generierung ===\n")
    
    # Zylinder-Mesh (benötigt mshr)
    try:
        mesh_cyl = generate_cylinder_mesh(resolution=64)
        print("✓ Zylinder-Mesh erfolgreich\n")
    except Exception as e:
        print(f"✗ Zylinder-Mesh fehlgeschlagen: {e}\n")
    
    # Kanal-Mesh (immer möglich)
    try:
        mesh_ch = generate_channel_mesh(nx=80, ny=16)
        print("✓ Kanal-Mesh erfolgreich")
    except Exception as e:
        print(f"✗ Kanal-Mesh fehlgeschlagen: {e}")
