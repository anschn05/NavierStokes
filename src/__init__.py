"""FEniCS Navier-Stokes Solver Package.

A modular package for solving Stokes and Navier–Stokes equations using FEniCS.
Includes solvers with Picard iteration, mesh generation, and postprocessing utilities.
"""

__version__ = "0.1.0"
__author__ = "Martin"

from .parameters import Params, default_params
from .solver import NavierStokesSolver
from .geometry import CylinderGeometry, get_mesh_bounds, reynolds_to_viscosity
from .boundary_conditions import BoundaryMarker, apply_bcs_cylinder_flow

__all__ = [
    "Params",
    "default_params",
    "NavierStokesSolver",
    "CylinderGeometry",
    "get_mesh_bounds",
    "reynolds_to_viscosity",
    "BoundaryMarker",
    "apply_bcs_cylinder_flow",
]
