# Primary Visium analysis

This folder contains the public analysis workflow for the primary 10x Visium melanoma dataset. The easiest way to reproduce the module is to use the sequential runner; the individual scripts can still be run separately when a particular step needs to be inspected.

## Quick start

From the project root:

```r
source(
  "scripts/analysis/12_Primary_Visium/run_12_Primary_Visium_full_sequential_reproduction.R",
  echo = TRUE
)
```

The runner follows the same order used in the final clean reproduction:

1. `01_Visium_build_object_and_state_scoring.R` builds the spatial object and scores the four tumor–immune states.
2. `02_Visium_spatial_QC_Moran_colocalization.R` performs spatial QC, QC-residualized analyses, and Moran-type spatial autocorrelation analyses, and provides the source results used for Figure 8A.
3. `03_Visium_dual_high_niche_analysis.R` evaluates the dual-high niche candidates used for Figure 8B–D.
4. `03B_Visium_make_Figure8.R` assembles the validated source tables and generates the final Figure 8 panels.
5. `04_Visium_threshold_independent_covariation.R` reproduces the historical 18D threshold-independent spatial co-variation analysis.
6. `05_Visium_composition_aware_analysis.R` runs the historical 18E composition-aware mixed-spot analysis.
7. `06_Visium_humoral_TLS_context.R` generates the historical 18F humoral/TLS-like spatial context used by the later supplementary analyses.
8. `07_Visium_primary_canonical_analysis.R` assembles the canonical statistics and source tables for Supplementary Figures S15–S19.
9. `08_Visium_make_S15_S19_figures.R` generates Supplementary Figures S15–S19.
10. `09_Visium_compare_with_frozen.R` compares the regenerated scientific tables with the frozen reference results.
11. `09A_Visium_diagnose_reproduction_mismatches.R` provides the final diagnostic comparison for any residual mismatch.

The 18D, 18E, and 18F analyses remain explicit steps because their outputs are genuine inputs to the final reporting workflow; they are not copied in as hidden intermediates.

## Inputs

The raw 10x files are expected under:

```text
data_raw/spatial_10x_human_melanoma_IF_FFPE/
```

The four-state gene sets are read from:

```text
results/tables/ICBcomb_final_input_gene_sets_CLEAN/
└── 02_final_ICBcomb_gene_sets/
    └── final_ICBcomb_gene_sets_CLEAN.gmt
```

The processed spatial object created in Step 01 is used by the downstream Visium analyses. The canonical processed Visium object is also distributed through the associated Zenodo processed-object release so that users can enter the workflow at a validated intermediate stage if they do not need to rebuild it from raw data.

## Figures and reporting

The four manuscript-facing state labels are:

- immune-defective/cold
- myeloid–Treg immunosuppressive
- tumor-dedifferentiation/stromal-remodeling
- melanocytic differentiation

Their categorical colors are `#4DBBD5`, `#00A087`, `#E64B35`, and `#3C5488`, respectively. Technical or QC panels do not reuse these colors where doing so could imply a biological state identity.

Final Figure 8 panels are written to:

```text
results/figures/spatial_melanoma_validation/
```

Supplementary Figures S15–S19 are written to:

```text
results/figures/primary_visium_spatial/
```

Intermediate diagnostic figures are kept under:

```text
results/diagnostics/12_Primary_Visium/
```

The publication figures use the same visual scale adopted across the revised project, approximately matching the proportions used for Supplementary Figure S20.

## Reproduction checks

A successful run retains the 3,458-spot primary Visium universe and produces all four state-score columns. The canonical reporting sources contain 3,458 spots in the S15 dominant-state table and 3,458 unique spots in the S19 source table.

The regenerated core scientific tables are compared with the frozen references using a numerical tolerance of `1e-10`, with exact matching for nonnumeric fields. Figure 8A–D and Supplementary Figures S15–S19 must also be present.

`09_Visium_compare_with_frozen.R` performs the formal frozen-output comparison, while `09A_Visium_diagnose_reproduction_mismatches.R` provides the diagnostic layer used to investigate any discrepancy. In the final clean sequential run, all Primary Visium steps completed successfully and no blocking mismatch was identified.

Public outputs use manuscript-facing filenames. A small number of legacy names are retained only where renaming them would break an upstream frozen dependency or a frozen-reference comparison; these naming exceptions do not change the statistical analysis.
