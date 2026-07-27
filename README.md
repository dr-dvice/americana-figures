# Figure scripts — Optic Lobe Evolution and GenoMiss paper

R and Python scripts used to generate the main-text figures (Figures 1–5) of the
paper. These scripts consume the outputs of the analysis pipelines (single-cell
atlas + SAMap, and pySCENIC regulon inference) plus other intermediates created 
during the analysis process and render the publication-formatted panels. 

Part of the project hub: https://dr-dvice.github.io/americana/

## Repository layout

```
figure1/   Intron-length distributions, GenoMiss figures, GO enrichment
figure2/   Schistocerca americana single cell atlas figures
figure3/   Neuromodulator / neurotransmitter panels for americana
figure4/   Other species UMAPs, conserved cell types
figure5/   pySCENIC analysis results
```

Figures are formatted for *Nature Ecology and Evolution* (Helvetica, 5–7 pt,
EPS + PNG output).

## Important: data paths

These scripts are provided **for transparency and reproducibility of the figure
code**. They contain **hardcoded absolute paths** to input objects on the authors' system, 
including:

- the labelled atlas Seurat objects (`.rds`),
- per-cell AUCell matrices (`.loom`) and regulon tables (`.csv`),
- per-species intron-length tables,
- cross-species homolog tables.

The large input objects are produced by the analysis pipelines (linked from the
project hub) and are **not redistributed here** (gonna update this when we upload our data). 

Small mapping/lookup tables that the scripts depend on **are** included
(`harmonization*.csv`, `cluster_mappings_centralized.xlsx`).

DISCLAIMER: These scripts do not represent the final published figures, only a starting point. Additional changes
were made in Adobe Illustrator.
