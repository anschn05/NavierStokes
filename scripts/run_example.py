"""Kleines Beispielskript zum Starten des Solvers mit einem mitgelieferten Mesh."""
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.parameters import default_params
from src.solver import NavierStokesSolver


def main():
    params = default_params()
    params.verbose = True

    # Get project root (parent of scripts/)
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    mesh_file = os.path.join(project_root, "meshes", "cylinder_flow.xml")
    if not os.path.exists(mesh_file):
        mesh_file = os.path.join(os.path.dirname(__file__), "meshes", "channel.xml")

    print(f"Using mesh: {mesh_file}")

    solver = NavierStokesSolver(mesh_file, params)
    try:
        w = solver.solve(save_path=os.path.join(project_root, "results", "solution.xdmf"))
        print("Solver finished successfully.")
    except Exception as e:
        print("Fehler beim Lösen:", e)


if __name__ == "__main__":
    main()
