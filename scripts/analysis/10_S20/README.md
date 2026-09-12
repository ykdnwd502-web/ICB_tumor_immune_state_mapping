# External spatial recurrence (Supplementary Figure S20)

This module rebuilds the external spatial recurrence analysis used for Supplementary Figure S20. It brings two independent melanoma spatial transcriptomics sources—GSE250636 Visium and the Thrane et al. 2018 legacy-ST dataset—through the historical 18J analysis, reconstructs the source-locked S20 tables, and then regenerates the combined reporting figure.

The analysis is intended as an **external spatial recurrence audit** of the predefined tumor–immune state programs. These datasets do not provide matched immune-checkpoint-blockade response labels, so this module should not be interpreted as an independent ICB-response or predictive-biomarker validation.

## Running the module

The preferred entry point is:

```r
source(
  "D:/ICB_resistance_project/scripts/analysis/10_S20/run_10_S20_full_sequential_reproduction.R",
  echo = TRUE
)
```

The project root can also be supplied through the `ICB_PROJECT_DIR` environment variable.

The runner executes the following sequence:

1. `01_GSE250636_build_18J.R`
2. `02_Thrane2018_build_18J.R`
3. `03_build_historical_S28_from_18J.R`
4. `04_GSE250636_build_S20_canonical.R`
5. `05_Thrane2018_build_S20_canonical.R`
6. `06_verify_S20_locked_state_genes.R`
7. `07_S20_aggregate.R`
8. `08_S20_make_figure.R`
9. `09_S20_external_spatial_recurrence_reproduction_check.R`

Each step is run in sequence from the project root, and the runner writes its session information to `logs/sessionInfo_10_S20_full_sequential_reproduction.txt`.

## Source data and locked inputs

GSE250636 is read from:

```text
data_external_spatial_18J/GSE250636/
```

with the public workflow using the unpacked GEO supplementary files under:

```text
unpacked_recursive/GSE250636_RAW/
```

The Thrane et al. 2018 data are read from:

```text
data_external_spatial_18J/Thrane2018_legacyST/
```

The public acquisition layer prepares the original Zenodo archive and its extracted count matrices before this module is run.

The four-state gene definition used for the final S20 rebuild is locked in:

```text
public_release/S20_external_spatial_recurrence/config/S20_locked_state_genes.csv
```

`06_verify_S20_locked_state_genes.R` checks that the state-gene tables reconstructed independently from GSE250636 and Thrane2018 agree with this locked source. It does not select, filter, or redefine state genes.

## What the scripts produce

Steps 01–03 reproduce the historical 18J/S28 layer retained for provenance and comparison. Steps 04–05 rebuild the canonical S20 outputs for the two datasets, and Step 07 combines them into the five reporting inputs used by Supplementary Figure S20.

The main canonical tables are written below:

```text
results/tables/revision_external_spatial_recurrence_18J/
```

The combined S20 reporting sources include dataset-level summaries, gene coverage, Spearman co-variation, dual-high enrichment, and bivariate Moran summaries. `08_S20_make_figure.R` uses these sources to generate the PNG, JPG, and PDF versions of Supplementary Figure S20 together with panel-specific source tables and a figure audit file.

## Reproduction check

The final checker compares the historical and rebuilt layers rather than simply testing for file existence. In the final sequential rerun, all nine scripts completed successfully. The audit reported:

- historical 18J outputs: 20/20 present;
- historical S28 core outputs: 10/10 present;
- canonical dataset outputs: 21/21 present;
- locked state-gene sources: 5/5 identical;
- historical-to-canonical numerical comparisons: 10/10 passed, with a maximum absolute numerical difference of `0`;
- aggregate outputs: 10/10 present;
- aggregate consistency gates: 7/7 passed;
- Supplementary Figure S20/reporting outputs: 9/9 present;
- figure audit consistency: passed.

The final decision was `READY_FOR_S20_ANALYSIS_FREEZE`.

The runner-level audit is written under:

```text
results/audit/10_S20_full_sequential_runner/
```

and the detailed reproduction audit under the S20 audit directory created by `09_S20_external_spatial_recurrence_reproduction_check.R`.
