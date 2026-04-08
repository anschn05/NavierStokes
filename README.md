# Navier-Stokes Simulation – NGSolve (FEM)

Dieser Branch enthält den Code der Navier-Stokes-Simulation mithilfe von **NGSolve**, einem auf der Finiten-Elemente-Methode (FEM) basierenden Framework. 

---

## Problembeschreibung

Simuliert wird eine instationäre, inkompressible Strömung um einen Kreiszylinder – ein klassisches Benchmark-Problem der Strömungsmechanik, bekannt für die dabei entstehende **Kármánsche Wirbelstraße**.

**Geometrie:**
- Rechteckiger Kanal: $2.0 \times 0.41$ m
- Kreisförmiges Hindernis: Mittelpunkt (0.5m, 0.2m), Radius 0.5m $(0.5 \mid 0.2)$ m, Radius $R = 0.05$ m

**Randbedingungen:**
- **Einlass:** konstante Geschwindigkeit, mit $U_\mathrm{max} = 5{,}0$ m/s
- **Wände & Zylinder:** No-Slip ($\mathbf{u} = \mathbf{0}$)
- **Auslass:** Do-Nothing-Randbedingung

**Physikalische Parameter:**
- Kinematische Viskosität: $\nu = 10^{-3}$ m²/s
- Reynoldszahl: $\mathrm{Re} = U_\mathrm{max} \cdot 2R / \nu$

---

## Numerische Methode

Die schwache Formulierung der Navier-Stokes-Gleichungen wird mit **Taylor-Hood-Elementen** (P2/P1) diskretisiert:
- Geschwindigkeit: quadratische Ansatzfunktionen (P2)
- Druck: lineare Ansatzfunktionen (P1)

Die Nichtlinearität wird pro Zeitschritt über eine **Picard-Iteration** aufgelöst. Das Mesh wird mit **Netgen** erzeugt.
---

## Voraussetzungen

```bash
pip install ngsolve numpy jupyter
```

NGSolve enthält Netgen bereits — keine separate Installation nötig.

---

## Ausführung

```bash
jupyter notebook navierstokes_chat_claude_ENDVERSION.ipynb
```

> Die Visualisierung (`Draw(...)`) funktioniert nur im Jupyter-Notebook! Ein Python-Interpreter reicht dafür nicht aus.

---

## Ergebnisse

Die Simulation visualisiert:
- Das **Geschwindigkeitsfeld** der zeitabhängigen Strömung
- Die **Wirbelstärke (Vorticity)** zur Darstellung der Kármánschen Wirbelstraße

---

## Autoren

Anna Sonnleitner, Emanuel Steininger, Martin Strobl
