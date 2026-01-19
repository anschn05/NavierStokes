# Devlog / Notizen — NavierStokes_OpenFOAM

## Ziel
2D (quasi-2D) Kanalströmung um Zylinder mit OpenFOAM 13 (`incompressibleFluid`) und `snappyHexMesh`.
Fokus: reproduzierbarer Case + Visualisierung (vorticity, Contours).

## Setup-Entscheidungen
- Solver: `incompressibleFluid`
- Geometrie: Kanal (L x H x dünne Dicke), Zylinder via STL in `constant/triSurface/`
- Meshing: `blockMesh` + `snappyHexMesh` (runder Zylinder)
- Postprocessing: Ausgabe von `vorticity` und `mag(vorticity)` über functionObjects

## Probleme & Lösungen (chronologisch)

### 1) snappyHexMesh: “Mesh provided is not fully 3D”
**Symptom:** snappy bricht ab, weil das Mesh als 2D (empty patches) erkannt wird.  
**Ursache:** Snapping/Relaxation braucht ein echtes 3D-Mesh.  
**Fix:** Kanal als dünnes 3D mit kleiner z-Dicke bauen (z.B. ±0.005), erst nach snappy ggf. auf quasi-2D umstellen.

### 2) snappyHexMesh: STL / Dictionary Fehler
**Symptom:** “keyword file undefined … cylinder.stl” (oder ähnlich).  
**Fix:** `geometry{ cylinder{ type triSurfaceMesh; file "cylinder.stl"; } }` korrekt setzen und Pfad `constant/triSurface/`.

### 3) changeDictionary / boundary patch types
**Symptom:** deprecated Tools / falsche patch types / (symmetryPlane nicht planar etc.).  
**Fix:** Patch-Typen konsistent halten (z.B. `patch` vs `symmetry` vs `empty`) und erst am Ende gezielt anpassen.

### 4) Solver-Start: transportProperties keys fehlen
**Symptom:** `viscosityModel` oder `nu` undefined.  
**Fix:** `constant/transportProperties` auf OF13-Format bringen:
- `viscosityModel constant;`
- `nu [0 2 -1 0 0 0 0] <Wert>;`

### 5) Laufzeit extrem langsam / deltaT sehr klein
**Symptom:** viele Schritte, `deltaT` klebt am `maxCo`.  
**Fix:** `maxCo`, `maxDeltaT`, PIMPLE-Korrekturen und `writeInterval` sinnvoll wählen.  
**Hinweis:** Für „schöne“ Wirbelstraße ist Re-Ziel (z.B. 80–200) wichtiger als nur U zu ändern.

### 6) ParaView UI/WSL: “Apply/Properties”/Warnfenster hängt
**Symptom:** Properties in eigenem Fenster blockiert, Warning Dialog nicht klickbar.  
**Fix:** Layout reset / Window-Manager Shortcuts (Alt+F4, Alt+Space → Move/Close), ggf. ParaView config reset.

## Prompt-/Arbeitsverlauf (Kurzform)
- Case-Struktur erstellt (0/, constant/, system/, scripts)
- blockMesh + snappyHexMesh Workflow (inkl. Allrun/Allclean)
- Debugging: 2D/3D snappy Issue, boundary patch types
- Solver-Fehler: transportProperties Format
- Performance-Tuning: maxCo/maxDeltaT/PIMPLE/writeInterval
- Postprocessing: vorticity + mag(vorticity), Slice + Contours

## Offene TODOs
- “schöne” ParaView State-Datei (`.pvsm`) speichern
- ggf. Grid refinement study (Qualität vs Laufzeit)
