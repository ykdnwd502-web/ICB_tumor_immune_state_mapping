# Additional robustness analyses

This folder contains the analyses added during revision to document patient-level structure, cross-cohort comparability and the project-wide statistical testing inventory. The module uses the already frozen outputs from the earlier analysis stages; it does not redefine the four tumor–immune states or modify the analytical thresholds used in modules 01–14.

## Running the module

From the project root:

```r
source(
  "scripts/analysis/15_Additional_robustness/run_15_Additional_robustness_full_sequential.R",
  echo = TRUE
)
```

A successful run ends with `Blocking failures: 0`, `Overall: PASS` and `Decision: READY_FOR_PUBLIC_FREEZE`.

## What the scripts do

The analysis scripts are:

- `01_GSE244983_patient_level_summary.R`, which derives the patient-level source tables used in S29 from the 26,053-cell GSE244983 analysis;
- `02_ICB_cross_cohort_comparability.R`, which evaluates GSE78220 and GSE91061 in their shared expression space;
- `03_statistical_testing_inventory.R`, which builds the statistical-testing inventory used in Table S29.

Reporting scripts then generate S29, S30 and Tables S29–S30:

- `01_make_S29_patient_level_structure.R`
- `02_make_S30_cross_cohort_robustness.R`
- `03_make_Table_S29_statistical_analysis_summary.R`
- `04_make_Table_S30_cross_cohort_comparability.R`

`15_Additional_robustness_reproduction_check.R` performs the final module audit.

## Figure conventions

S29 and S30 use the same publication-scale typography as the rest of the revised figures. The four state colors remain `#4DBBD5`, `#00A087`, `#E64B35` and `#3C5488` when state identity is represented categorically. The patient-state heatmap in S29 intentionally uses a continuous diverging scale because the fill encodes a numerical mean score rather than a categorical state.

Public S29/S30 filenames use the manuscript-facing names. Development suffixes such as `CANONICAL`, `FINAL` and `FIXED` are not added to regenerated public outputs. Historical upstream filenames are left unchanged when they are part of an already frozen dependency.
