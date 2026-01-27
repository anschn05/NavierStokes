"""Geometry utilities and mesh handling.

Provides geometry classes, mesh bounds calculation, and flow profile functions
for standard benchmarks like cylinder flow (DFG 2D-1/2D-2/2D-3).
"""
from dolfin import *
import numpy as np


class CylinderGeometry:
    """Geometrie-Parameter für Zylinderumströmung (DFG-Benchmark).
    
    Standard 2D cylinder flow benchmark geometry:
    - Rectangular channel with inlet/outlet
    - Cylinder obstacle in center
    """
    
    # Modified for Vortex Shedding (cylinder further from inlet)
    L = 2.2          # Kanallänge [m]
    H = 0.41         # Kanalhöhe [m]
    cx = 0.5         # Zylinder x-Position [m] (5D vom Einlass)
    cy = 0.2         # Zylinder y-Position [m] (zentriert)
    r = 0.05         # Zylinderradius [m]
    D = 2 * r        # Durchmesser = 0.1 [m]
    
    @classmethod
    def info(cls):
        """Drucke Geometrie-Info."""
        print(f"Cylinder Flow Geometry (DFG Benchmark):")
        print(f"  Channel: L={cls.L}, H={cls.H}")
        print(f"  Cylinder: center=({cls.cx}, {cls.cy}), r={cls.r}, D={cls.D}")


def get_mesh_bounds(mesh):
    """Extrahiere Bounding Box des Meshes.
    
    Args:
        mesh: FEniCS Mesh-Objekt
    
    Returns:
        dict: {'xmin', 'xmax', 'ymin', 'ymax'} mit Koordinaten-Grenzen
    """
    coords = mesh.coordinates()
    return {
        'xmin': float(np.min(coords[:, 0])),
        'xmax': float(np.max(coords[:, 0])),
        'ymin': float(np.min(coords[:, 1])),
        'ymax': float(np.max(coords[:, 1]))
    }


def get_cylinder_center(mesh=None):
    """Finde oder gebe Zylinder-Mittelpunkt zurück.
    
    Args:
        mesh: Optional FEniCS Mesh (für zukünftige auto-detect)
    
    Returns:
        tuple: (cx, cy) Mittelpunkt
    """
    return (CylinderGeometry.cx, CylinderGeometry.cy)


def reynolds_to_viscosity(Re, D=CylinderGeometry.D, U=1.0):
    """Berechne kinematische Viskosität aus Reynolds-Zahl.
    
    Re = U * D / ν  =>  ν = U * D / Re
    
    Args:
        Re: Reynolds-Zahl
        D: Charakteristische Länge (Standard: Zylinderdurchmesser)
        U: Charakteristische Geschwindigkeit (Standard: Einlass)
    
    Returns:
        float: Kinematische Viskosität ν
    """
    return U * D / Re


def viscosity_to_reynolds(nu, D=CylinderGeometry.D, U=1.0):
    """Berechne Reynolds-Zahl aus kinematischer Viskosität.
    
    Re = U * D / ν
    
    Args:
        nu: Kinematische Viskosität
        D: Charakteristische Länge
        U: Charakteristische Geschwindigkeit
    
    Returns:
        float: Reynolds-Zahl
    """
    return U * D / nu


def velocity_profile_parabolic(x, U_mean, H):
    """Parabolisches Geschwindigkeitsprofil (Poiseuille).
    
    Für 2D Kanal (oder 3D Rohr):
    u(y) = 4 * U_mean * y * (H - y) / H²
    
    Maximum bei y = H/2: u_max = U_mean
    
    Args:
        x: Punkt-Koordinaten (array-like)
        U_mean: Mittlere Geschwindigkeit
        H: Kanalhöhe
    
    Returns:
        float: Geschwindigkeit an Position
    """
    y = x[1] if len(x) > 1 else 0
    if 0 <= y <= H:
        return 4 * U_mean * y * (H - y) / (H * H)
    return 0


class ChannelFlow:
    """Einfacher Kanal-Strömung (ohne Hindernis)."""
    
    L = 2.2
    H = 0.41
    
    @classmethod
    def info(cls):
        print(f"Channel Flow Geometry:")
        print(f"  L={cls.L}, H={cls.H}")
