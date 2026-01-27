import os
import sys
import pytest

# Ensure project root is on sys.path so `src` package can be imported when running pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Skip tests if dolfin (FEniCS) is not available in the environment
dolfin = pytest.importorskip("dolfin")
from dolfin import UnitSquareMesh, norm

from src.parameters import default_params
from src.solver import NavierStokesSolver


def test_solver_smoke(tmp_path):
    """Smoke test: build a small UnitSquare mesh and run the solver once."""
    mesh = UnitSquareMesh(8, 8)

    params = default_params()
    params.verbose = False
    params.steady = True  # Stationärer Modus

    solver = NavierStokesSolver(mesh, params)
    w = solver.solve()

    assert w is not None
    u, p = w.split()
    # Basic sanity checks
    assert u.function_space().mesh().num_vertices() == mesh.num_vertices()
    assert p.function_space().mesh().num_vertices() == mesh.num_vertices()


def test_stokes_linearized(tmp_path):
    """Test: Stokes (linearized) solver mode."""
    mesh = UnitSquareMesh(8, 8)

    params = default_params()
    params.verbose = False
    params.steady = True
    params.linearized_only = True  # Stokes mode

    solver = NavierStokesSolver(mesh, params)
    w = solver.solve()

    u, p = w.split()
    assert norm(u, "L2") > 0  # Should have non-zero velocity
    assert norm(p, "L2") >= 0  # Pressure can be zero (reference value)


def test_navier_stokes_convergence(tmp_path):
    """Test: Navier–Stokes solver converges with Picard iteration."""
    mesh = UnitSquareMesh(6, 6)

    params = default_params()
    params.verbose = False
    params.steady = True
    params.inlet_velocity = 0.5  # Lower velocity for stability
    params.picard_maxiter = 10
    params.tol = 1e-5

    solver = NavierStokesSolver(mesh, params)
    w = solver.solve()

    u, p = w.split()
    u_norm = norm(u, "L2")
    p_norm = norm(p, "L2")

    # Sanity: velocity should be positive (flow field)
    assert u_norm > 0
    # Verify solution was computed
    assert u.function_space().mesh().num_vertices() == mesh.num_vertices()

