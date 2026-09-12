# ICBcomb analysis

`ICBcomb/` is a self-contained analysis module that sits alongside the main top-level `scripts/` directory.

It is kept separate because this part of the study depends on historically fixed ICBcomb web-platform query results rather than on a raw sequencing dataset that can be downloaded and regenerated from scratch. The public release therefore preserves the exact query inputs and platform outputs used in the study and rebuilds the manuscript-facing ICBcomb figures and tables from that frozen source.

The module does **not** submit new queries to the current ICBcomb website.

## Directory layout

```text
ICBcomb/
├── README.md
├── run_ICBcomb_full_sequential.R
├── inputs/
│   └── ICBcomb_INPUT_FREEZE_v1.0/
├── scripts/
│   ├── audit/
│   ├── analysis/
│   ├── figures/
│   └── tables/
└── results/
```

`inputs/ICBcomb_INPUT_FREEZE_v1.0/` contains the frozen query-gene-set bundle, the three raw platform result files, historical reference tables, and the corresponding integrity/audit records. These files define the historical ICBcomb input boundary and should not be edited.

The `results/` subdirectories are created or populated by the runner and do not need to contain pre-generated outputs before a clean run.

## Running the module

From the project root, run:

```r
source(
  "ICBcomb/run_ICBcomb_full_sequential.R",
  echo = TRUE
)
```

The runner resolves `ICBcomb/` as its own module root and then executes:

1. `scripts/audit/01_ICBcomb_source_lock_audit.R`
2. `scripts/analysis/02_ICBcomb_build_source_tables.R`
3. `scripts/figures/03_make_Figure9_ICBcomb.R`
4. `scripts/figures/04_make_Supplementary_Figure_S27.R`
5. `scripts/figures/05_make_Supplementary_Figure_S28.R`
6. `scripts/tables/06_make_Supplementary_Tables_S27A_S27B.R`
7. `scripts/audit/07_ICBcomb_post_rebuild_audit.R`

The runner sets the module-root environment variables used by the child scripts, so the ICBcomb scripts remain independent of the main `scripts/analysis/` hierarchy.

## Frozen input boundary

The frozen input package contains exactly three raw ICBcomb platform result files, corresponding to the three resistance-state queries:

- immune-defective/cold — RESTORE
- myeloid–Treg immunosuppressive — SUPPRESS
- tumor-dedifferentiation/stromal-remodeling — SUPPRESS

The melanocytic-differentiation program is retained as the reference/control state and is not submitted as an independent reversal query.

The historical access date and query configuration are preserved inside the frozen input package. Re-running the live web platform is intentionally excluded because later changes in the platform or its underlying databases could alter the returned results.

## Outputs

The module rebuilds the source tables used for:

- Figure 9
- Supplementary Figure S27
- Supplementary Figure S28
- Supplementary Table S27A
- Supplementary Table S27B

Generated files are written below `ICBcomb/results/`, while the main project-level session information for the runner is archived with the release metadata.

## Reproduction status

The final clean run completed all seven steps successfully. The source-lock audit passed, all three raw platform files were recovered as 45 total rows, the rebuilt numerical reference comparisons passed, and the post-rebuild audit passed all 6/6 gates.

The final module decision was:

```text
READY_FOR_MODULE_FREEZE
```

The ICBcomb results should be interpreted as exploratory database-derived perturbation-prioritization signals. They are not experimental validation, therapeutic efficacy evidence, clinical actionability, or causal proof.
