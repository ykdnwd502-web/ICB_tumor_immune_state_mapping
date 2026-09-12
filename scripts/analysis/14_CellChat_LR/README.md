# CellChat and CosMx ligand–receptor analysis

This module reproduces the CellChat/CosMx ligand–receptor analysis used for Supplementary Figures S24–S26 and Supplementary Tables S24–S26. The final step compares the regenerated results with the frozen reference outputs.

## Running the module

From the project root, run:

```r
source(
  "scripts/analysis/14_CellChat_LR/00_RUN_CellChat_LR_SEQUENTIAL.R",
  echo = TRUE
)
```

If the frozen reference files are stored outside the project tree, set `ICB_FROZEN_ROOT` before running. The project root can likewise be supplied through `ICB_PROJECT_DIR` when the repository is not located at the default path.

The runner executes ten steps: CellChat candidate nomination; CosMx spatial ligand–receptor analysis at `k = 20`; the `k = 10` and `k = 30` sensitivity analyses; consolidation of reporting sources; generation of S24, S25 and S26; generation of Tables S24–S26; and the final current-versus-frozen reproduction check.

## Analysis settings

The CellChat analysis uses seed `20260504`, with at most 800 cells per group and at least 30 cells per group. `CellChatDB.human` is used in full, communication probabilities are estimated with `truncatedMean` (`trim = 0.1`), and `filterCommunication(min.cells = 30)` is applied. Twelve predefined axes are carried through the analysis.

For the CosMx spatial analysis, `k = 20` is the primary neighborhood size and `k = 10`/`30` are sensitivity settings. Each `k` uses 999 permutations. High-expression cells are defined as the top 25% among expressing cells; when fewer than 50 cells are positive, all positive cells are treated as high. At least 100 eligible edges are required. The analysis does not cap the cell number used for kNN construction (`max_cells_for_knn = Inf`). Sixty tests are evaluated per `k` (12 axes × 5 contrasts), with BH correction applied separately within each `k`-specific testing universe.

## Frozen inputs and checks

The files under `locked_input/` preserve the historical CellChat cell selection, while `frozen_reference/` contains the summary outputs used for the final comparison. These files are inputs to the reproduction check and should not be edited.

A complete public run should finish with the module-level reproduction audit passing and should regenerate S24–S26 and Tables S24–S26 without blocking errors.
