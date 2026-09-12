# CosMx spatial validation

This module evaluates the predefined tumor–immune state programs in the CosMx melanoma spatial dataset. It starts from the Slide 4 CosMx expression matrix, reconstructs cell-level state scores, adds spatial coordinates and metadata/QC information, and then regenerates Supplementary Figures S21 and S22.

The workflow separates the scoring, spatial validation, metadata sensitivity, and reporting steps so that each layer can be inspected independently.

## Running the module

Use the sequential runner:

```r
source(
  "D:/ICB_resistance_project/scripts/analysis/11_CosMx/run_11_CosMx_full_sequential_reproduction.R",
  echo = TRUE
)
```

The runner executes:

1. `01_CosMx_state_scoring.R`
2. `02_CosMx_cell_level_spatial_validation.R`
3. `03_CosMx_metadata_QC_and_technical_sensitivity.R`
4. `04_CosMx_make_S21.R`
5. `05_CosMx_make_S22.R`
6. `11_CosMx_reproduction_check.R`

The project root defaults to `D:/ICB_resistance_project` and can be overridden with `ICB_PROJECT_DIR`.

## Input data

The public workflow uses Slide 4 of the Dryad CosMx release, specifically the `Run5611_MK3` files prepared by the data-acquisition scripts:

```text
data_raw/CosMx_melanocytic_tumors_Dryad/
└── Slide_4/
    └── Run5611_MK3/
        ├── Run5611_MK3_exprMat_file.csv
        ├── Run5611_MK3_metadata_file.csv
        ├── Run5611_MK3_fov_positions_file.csv
        └── Run5611_MK3_tx_file.csv
```

The expression matrix is used to reconstruct the four state scores. Spatial coordinates and the CosMx metadata file are then merged with the scored cells for the spatial and technical-sensitivity analyses.

## Analysis flow

`01_CosMx_state_scoring.R` scores the four predefined state programs at cell level and writes the state-score table and gene-overlap audit.

`02_CosMx_cell_level_spatial_validation.R` adds the spatial information and evaluates dominant-state distributions, dual-high categories, and pooled/per-FOV threshold sensitivity. The primary dual-high calculation uses a top-25% threshold, with 20%, 25%, and 30% thresholds retained for sensitivity analysis.

`03_CosMx_metadata_QC_and_technical_sensitivity.R` joins the CosMx metadata and evaluates whether state-score relationships are robust to available technical/QC variables. Its spatially resolved output is the reporting source used by S21 and S22.

`04_CosMx_make_S21.R` reconstructs Supplementary Figure S21 from the metadata-QC merged table. `05_CosMx_make_S22.R` reconstructs Supplementary Figure S22, including the pooled versus per-FOV threshold-sensitivity panel.

The main tables are written under:

```text
results/tables/CosMx_spatial_validation/
results/tables/CosMx_metadata_QC/
```

and the canonical S21/S22 reporting products are written under the corresponding canonical/figure output directories used by the frozen project.

## Reproduction check

The clean sequential rerun completed all six steps successfully. The final log recorded 86,600 scored CosMx cells, of which 86,572 had the spatial/metadata information required for the downstream spatial analyses.

The reproduction audit checks four points:

- the Step 01 raw scored-cell universe is 86,600 cells;
- state scores are numerically preserved from Step 01 to Step 02;
- state scores are numerically preserved after the Step 03 metadata/QC merge;
- the six S22 threshold-sensitivity combinations (pooled/per-FOV × 20%/25%/30%) reproduce within the predefined reporting tolerance.

In the final run, the Step 01→02 maximum difference was `0`, the Step 02→03 maximum difference was approximately `2.22e-16`, and all six S22 threshold combinations matched. The maximum reporting-level odds-ratio difference was `0.0151`, below the predefined tolerance of `0.02`. There were no blocking failures.

The runner writes:

```text
results/audit/11_CosMx_full_sequential_runner/
    11_CosMx_full_sequential_runner_log.csv
```

and saves the environment record to:

```text
logs/sessionInfo_11_CosMx_full_sequential_reproduction.txt
```
