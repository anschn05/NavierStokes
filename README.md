# FEniCS Navier–Stokes Solver

Ein Projekt zur Lösung der stationären Navier–Stokes-Gleichungen mit FEniCS (DOLFIN), Python und Picard-Iteration.

## Features

- **Stokes-Solver** (linear, für Referenzlösungen): schnelle, stabile Grundimplementierung
- **Navier–Stokes-Solver** (nichtlinear) mit **Picard-Linearisierung**: Konvektion $(\mathbf{u} \cdot \nabla) \mathbf{u}$ iterativ gelöst
- **Taylor–Hood Elemente** (P2–P1): stabiles Geschwindigkeit-Druck-Paar für inkompressible Strömungen
- **Einfache Randbedingungen**: Einlass (Dirichlet), Wände (no-slip), Auslass (Druck)
- **Automatisierte Tests**: Smoke-Tests + Konvergenzprüfung (pytest)

## Verwendung

### Installation

**Voraussetzung**: FEniCS/DOLFIN (klassisch) oder DOLFINx (modern). Empfohlene Installation:

```bash
# Conda (schnellste Option)
conda create -n fenics -c conda-forge fenics=2019.1 matplotlib
conda activate fenics

# Oder Docker
docker pull dolfinx/dolfinx:stable
```

**Python-Abhängigkeiten**:
```bash
pip install -r requirements.txt
```

### Schnellstart

#### 1. Zylinderumströmung (empfohlen für schnellen Einstieg)
```bash
# Mesh generieren
python3 src/mesh_generator.py

# Simulation starten (transient, Re=100, Wirbelablösung)
python3 simulate_cylinder.py
```
Dies simuliert die **Strömung um einen Zylinder** mit:
- Transienter Navier–Stokes-Lösung (Zeitschritte)
- Berechnung von Drag- und Lift-Kräften
- Visualisierung der Kräfte über die Zeit
- Speicherung in `results/cylinder_flow.xdmf` (mit Paraview öffnen)

**Erwartetes Ergebnis bei Re=100**: Periodische **Wirbelablösung** (Kármán-Wirbelstraße)

#### 2. Einfaches Beispiel (Stokes, stationär)
```bash
python3 run_example.py
```
Löst die Stokes-Gleichung auf vorhandenem Mesh und speichert die Lösung in `results/solution.xdmf`.

#### 2. Navier–Stokes mit Picard-Iteration

```python
from src.parameters import Params
from src.solver import NavierStokesSolver

# Parameter setzen
params = Params(
    inlet_velocity=1.0,
    nu=1e-3,           # kinematische Viskosität
    picard_maxiter=20, # Picard-Iterationen
    tol=1e-6,          # Konvergenztolerantz
    verbose=True
)

# Solver initialisieren und lösen
solver = NavierStokesSolver("meshes/cylinder_flow.xml", params)
w = solver.solve(save_path="results/solution.xdmf")

# Ergebnis extrahieren
u, p = w.split()
print(f"Velocity norm: {norm(u, 'L2'):.6e}")
print(f"Pressure norm: {norm(p, 'L2'):.6e}")
```

#### 3. Stokes-Modus (linear)

```python
params.linearized_only = True  # Nur Stokes
solver = NavierStokesSolver(mesh_file, params)
w = solver.solve()
```

## Projektstruktur

```
fenics_navier-stokes/
├── src/
│   ├── solver.py           # Hauptsolver: Stokes + Navier–Stokes (Picard)
│   ├── parameters.py       # Params-Dataclass mit Defaults
│   ├── boundary_conditions.py  # (für zukünftige Erweiterungen)
│   ├── geometry.py         # (für zukünftige Erweiterungen)
│   └── utils.py            # Hilfsfunktionen
├── tests/
│   └── test_solver.py      # Unit-Tests (Stokes, NS Konvergenz)
├── meshes/
│   ├── cylinder_flow.xml   # Beispiel-Mesh: Zylinder-Umströmung
│   └── channel.xml         # Beispiel-Mesh: Kanalströmung
├── notebooks/
│   ├── 01_problem_setup.ipynb     # (zu befüllen)
│   ├── 02_simulation.ipynb        # (zu befüllen)
│   └── 03_visualization.ipynb     # (zu befüllen)
├── results/                # Ausgabeverzeichnis (wird erstellt)
├── run_example.py          # Beispielskript
├── requirements.txt        # Python-Abhängigkeiten
└── README.md              # Diese Datei
```

## Tests

Alle Tests mit pytest ausführen:

```bash
pytest -q
```

Spezifische Tests:

```bash
# Smoke-Test: Basis-Funktionalität
pytest tests/test_solver.py::test_solver_smoke -v

# Stokes-Modus (linear)
pytest tests/test_solver.py::test_stokes_linearized -v

# Navier–Stokes Konvergenz
pytest tests/test_solver.py::test_navier_stokes_convergence -v
```

**Aktuelle Test-Status**: ✅ 3/3 bestanden

## Mathematischer Hintergrund

### Stokes-Gleichungen (linear)
$$
-\nu \Delta \mathbf{u} + \nabla p = \mathbf{f}
$$
$$
\nabla \cdot \mathbf{u} = 0
$$

### Navier–Stokes-Gleichungen (nichtlinear)
$$
-\nu \Delta \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} + \nabla p = \mathbf{f}
$$
$$
\nabla \cdot \mathbf{u} = 0
$$

**Lösung mit Picard-Linearisierung**:
Iterative Lösung der bilinearen Form
$$
a_n(\mathbf{u}, p; \mathbf{v}, q) = \nu (\nabla \mathbf{u}, \nabla \mathbf{v}) + ((\mathbf{u}_{n-1} \cdot \nabla) \mathbf{u}, \mathbf{v}) - (\nabla \cdot \mathbf{v}, p) - (q, \nabla \cdot \mathbf{u})
$$
wobei $\mathbf{u}_{n-1}$ aus der vorherigen Iteration stammt.

## Parameter

Wichtige Einstellungen in `src/parameters.py`:

| Parameter | Default | Bedeutung |
|-----------|---------|-----------|
| `inlet_velocity` | 1.0 | Einlass-Geschwindigkeit [m/s] |
| `nu` | 1e-3 | kinematische Viskosität [m²/s] |
| `picard_maxiter` | 20 | Max Picard-Iterationen |
| `tol` | 1e-6 | Konvergenztolerantz |
| `linearized_only` | False | Stokes-Modus (True) oder Navier–Stokes (False) |
| `verbose` | True | Ausgabepegel |

## Beispiele

### Zylinderumströmung (Flow around Cylinder)

Das Hauptbeispiel `simulate_cylinder.py` demonstriert:

- **Geometrie**: Rechteck-Kanal (2.2×0.41) mit Zylinder (Radius 0.05 bei (0.2, 0.2))
- **Reynolds-Zahlen**:
  - Re=20: Stationär, symmetrisch (keine Wirbelablösung)
  - Re=100: **Instabil, periodische Wirbelablösung** (Kármán-Wirbelstraße)
  - Re>200: Turbulenter Nachlauf
- **Messung**: Drag- und Lift-Koeffizienten über Zeit

**Beispiel-Kommandos**:
```python
# Transient, Re=100, 8 Sekunden
python3 -c "from simulate_cylinder import run_cylinder_simulation; run_cylinder_simulation(transient=True, T=8.0, dt=0.01, Re=100)"

# Stationär, Re=20
python3 -c "from simulate_cylinder import run_cylinder_simulation; run_cylinder_simulation(transient=False, Re=20)"
```

**Ausgaben**:
- `results/cylinder_flow.xdmf`: Strömungsfeld (öffnen mit Paraview)
- `results/cylinder_forces.png`: Drag/Lift-Zeitreihe


## Anforderungen

- Python 3.8+
- FEniCS 2019.1 oder neuer (mit DOLFIN)
- numpy, scipy
- pytest (für Tests)
- Jupyter (optional, für Notebooks)

## Lizenz

MIT License — siehe `LICENSE` für Details.

## Autor

Martin (Seminar FEniCS Navier–Stokes)
python poisson_test.py
```

Dies wird:
- Die Poisson-Gleichung auf einem 2D-Mesh lösen
- Statistiken ausgeben
- Eine Visualisierung als `poisson_solution.png` speichern

## Weitere Ressourcen

- [FEniCS Tutorial](https://fenicsproject.org/tutorial/)
- [FEniCS Documentation](https://fenicsproject.org/documentation/)
- [FEniCSx Documentation](https://docs.fenicsproject.org/)
