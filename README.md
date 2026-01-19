# NavierStokes_OpenFOAM — 2D-Kanalströmung mit Zylinder (OpenFOAM 13)

OpenFOAM-13-Case für eine inkompressible Kanalströmung um einen kreisförmigen Zylinder (Hindernis) mit dem modularen
Solver `incompressibleFluid` und `snappyHexMesh`. Für die Visualisierung wird zusätzlich die Wirbelstärke (Vorticity)
ausgegeben, um den Nachlauf (von-Kármán-Wirbelstraße) gut darstellen zu können.

## Ordnerstruktur
- `0/` : Anfangs- und Randbedingungen (z.B. `U`, `p`)
- `constant/` : physikalische Eigenschaften und Geometrie (`transportProperties`, `triSurface/`)
- `system/` : Numerik- und Laufzeiteinstellungen (`controlDict`, `fvSchemes`, `fvSolution`, Meshing-Dictionaries)
- `Allrun`, `Allclean` : Hilfsskripte zum Ausführen und Aufräumen

## Voraussetzungen
- OpenFOAM 13 (getestet unter Linux/WSL)
- ParaView / `paraFoam` für die Visualisierung

## Ausführen

Nur Mesh erzeugen:
```bash
./Allclean
./Allrun -mesh
```

Mesh + Solver:
```bash
./Allclean
./Allrun
```
### Postprocessing (ParaView)

Case öffnen:
```bash
paraFoam
```

Empfohlen:

- `mag(vorticity)` auf einer Mittelschnitt-Ebene (Slice) anzeigen (dünnes 3D / quasi-2D)
- Optional `Contour` hinzufügen, um „Landkarten-ähnliche“ Isolinien zu erhalten


## Ressourcen

- OpenFOAM 13 Dokumentation (User Guide)
- ParaView Dokumentation (Filters: Slice, Contour)
- snappyHexMesh: grundlegender Workflow (castellated → snap → layers)
