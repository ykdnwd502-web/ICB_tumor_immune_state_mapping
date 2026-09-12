# IF/mIF marker-signal validation

This module links the CosMx-derived tumor–immune state scores to the protein/marker signals available in the matched IF/mIF metadata. It rebuilds the historical Step14B cell-level state scores, reproduces the Step14A matched-cell and correlation analysis, and generates Supplementary Figure S23 and Supplementary Table S23.

The intended interpretation is deliberately limited: this analysis provides **protein/marker-signal association support** for the predefined state programs. It is not a functional validation experiment.

## Running the module

The complete public chain is run with:

```r
source(
  "D:/ICB_resistance_project/scripts/analysis/13_IF_mIF/00_RUN_IF_mIF_SEQUENTIAL.R",
  echo = TRUE
)
```

The runner executes:

1. `01_IF_mIF_build_CosMx_state_scores.R`
2. `02_IF_mIF_match_and_correlation.R`
3. `03_IF_mIF_make_S23.R`
4. `04_IF_mIF_make_Supplementary_Table_S23.R`
5. `05_IF_mIF_compare_with_frozen.R`

Each child script is sourced in a separate environment. The runner records a console log and execution table under:

```text
results/audit/13_IF_mIF_runner/
```

and saves the session information as:

```text
logs/sessionInfo_13_IF_mIF_full_sequential_reproduction.txt
```

## Inputs and historical reference

The state-score producer uses the same CosMx Slide 4 expression source used by the CosMx module. It reconstructs the four cell-level state scores with the historical Step14B logic: the composite cell key is based on `fov + cell_ID`, expression is log1p transformed, genes are z-scaled before averaging, and the locked state direction is applied when available.

The matched IF/mIF analysis then audits candidate metadata files by FOV/cell overlap. The metadata source with the largest overlap is selected, and the raw versus normalized composite-key strategy is chosen by overlap rather than manually. The final matched-cell analysis requires at least 1,000 cells.

The frozen comparison files are stored under:

```text
results/tables/IF_mIF_validation/frozen_reference/
```

and include the clean historical Step14B/Step14A references used by the final comparison script.

## Correlation analysis and reporting

`02_IF_mIF_match_and_correlation.R` merges the state scores with the selected IF/mIF marker metadata and computes Spearman correlations using `cor.test(..., exact = FALSE)`. Benjamini–Hochberg correction is applied across the complete protein-feature × state comparison set.

The resulting tables are written under:

```text
results/tables/IF_mIF_validation/
```

The main outputs include the matched state/marker table and the protein-state Spearman summary.

`03_IF_mIF_make_S23.R` builds Supplementary Figure S23 from the fixed correlation source. The final reporting source contains 32 rows representing 8 marker features across 4 states, based on 86,572 matched cells.

`04_IF_mIF_make_Supplementary_Table_S23.R` exports the corresponding workbook as a clean XLSX file with the sheets:

```text
README
IF_mIF_Spearman
S23_harmonized
Source_manifest
```

## Frozen-output comparison

`05_IF_mIF_compare_with_frozen.R` performs the final release audit against three historical reference layers:

- Step14B CosMx cell-state scores, keyed by `.merge_id`;
- Step14A IF/state merged table, keyed by `.key_raw`;
- Step14A protein-state correlation summary.

For each layer, the current rebuild must agree with the frozen reference to a numerical tolerance of `1e-12`. The comparison is intended to establish that the public scripts reproduce the historical analytical results while using the cleaned public file names and project structure.

The final validated release reproduced the frozen state-score, merged-table, and correlation layers with no blocking differences, and the reporting products for Supplementary Figure/Table S23 were present. The module was therefore cleared for public freeze.

The detailed comparison table is written to the IF/mIF audit directory as:

```text
13_IF_mIF_frozen_comparison_v1.5_FINAL.csv
```
