"""Standard-Parameter für das Navier-Stokes Projekt.

Dieses Modul hält eine simple `Params`-Datenstruktur (als dict-kompatibel),
die in den Beispielskripten und dem Solver verwendet wird.
"""
from dataclasses import dataclass


@dataclass
class Params:
    # Physikalische Parameter
    nu: float = 1e-3       # kinematische Viskosität [m²/s]

    # Solver / numerische Parameter
    steady: bool = True    # stationär (True) oder zeitabhängig (False)
    add_perturbation: bool = False  # Störung für Vortex Shedding hinzufügen
    dt: float = 0.01       # Zeitschritt (nur bei transient)
    t_end: float = 1.0     # Endzeit (nur bei transient)
    tol: float = 1e-6      # Residuen-Toleranz für Picard
    picard_maxiter: int = 20  # Max-Iterationen für Picard-Linearisierung

    # Randbedingungs-Parameter
    inlet_velocity: float = 1.0
    outlet_pressure: float = 0.0

    # Ausgabe-Optionen
    verbose: bool = True
    linearized_only: bool = False  # Nur Stokes lösen (Testing)


def default_params() -> Params:
    return Params()
