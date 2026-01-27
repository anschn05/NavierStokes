"""Boundary condition definitions and utilities."""
from dolfin import *
import numpy as np


class BoundaryMarker:
    """Markiert verschiedene Rand-Regionen basierend auf Koordinaten.
    
    Beispiel:
        marker = BoundaryMarker(mesh)
        marker.mark_inlet(Inlet())
        marker.mark_walls(Walls())
        ds_inlet, ds_outlet, ds_wall = marker.get_measures()
    """
    
    def __init__(self, mesh):
        """Initialisiere Boundary-Marker.
        
        Args:
            mesh: FEniCS Mesh-Objekt
        """
        self.mesh = mesh
        self.boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
        self.boundaries.set_all(0)
        self.inlet_id = 1
        self.outlet_id = 2
        self.wall_id = 3
    
    def mark_inlet(self, condition):
        """Markiere Einlass-Boundary.
        
        Args:
            condition: SubDomain-Objekt für Einlass
        """
        condition.mark(self.boundaries, self.inlet_id)
    
    def mark_outlet(self, condition):
        """Markiere Auslass-Boundary.
        
        Args:
            condition: SubDomain-Objekt für Auslass
        """
        condition.mark(self.boundaries, self.outlet_id)
    
    def mark_walls(self, condition):
        """Markiere Wand-Boundary.
        
        Args:
            condition: SubDomain-Objekt für Wände
        """
        condition.mark(self.boundaries, self.wall_id)
    
    def get_measures(self):
        """Rückgabe von Integrations-Maßen für markierte Boundaries.
        
        Returns:
            tuple: (ds_inlet, ds_outlet, ds_wall)
        """
        ds = ds(subdomain_data=self.boundaries)
        return ds(self.inlet_id), ds(self.outlet_id), ds(self.wall_id)


def apply_bcs_cylinder_flow(W, mesh, inlet_velocity, outlet_pressure):
    """Standard Randbedingungen für Zylinderumströmung.
    
    - Einlass (xmin): Dirichlet u = (inlet_velocity, 0)
    - Wände (ymin, ymax): no-slip u = (0, 0)
    - Auslass (xmax): Druck-Dirichlet p = outlet_pressure
    
    Args:
        W: Gemischter Funktionsraum (V×Q)
        mesh: FEniCS Mesh
        inlet_velocity: Einlass-Geschwindigkeit (float)
        outlet_pressure: Auslass-Druck (float)
    
    Returns:
        list: Liste von DirichletBC-Objekten
    """
    coords = mesh.coordinates()
    xmin = float(np.min(coords[:, 0]))
    xmax = float(np.max(coords[:, 0]))
    ymin = float(np.min(coords[:, 1]))
    ymax = float(np.max(coords[:, 1]))
    eps = 1e-8 * max(1.0, abs(xmax - xmin), abs(ymax - ymin))
    
    class Inlet(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0] < xmin + eps
    
    class Outlet(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0] > xmax - eps
    
    class Walls(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (x[1] < ymin + eps or x[1] > ymax - eps)
    
    bcs = [
        DirichletBC(W.sub(0), Constant((inlet_velocity, 0.0)), Inlet()),
        DirichletBC(W.sub(0), Constant((0.0, 0.0)), Walls()),
        DirichletBC(W.sub(1), Constant(outlet_pressure), Outlet())
    ]
    return bcs
