# GSE244983 doublet and CNV sensitivity analysis

This module contains the revision-stage sensitivity analyses used to check whether the main GSE244983 results are robust to doublet removal and to obtain exploratory CNV-like support with inferCNV. It starts from the current state-localized GSE244983 Seurat object (26,053 cells and 11,616 features). The historical 25,972-cell branch is not used here, and none of the analyses in this folder replace the primary 26,053-cell analysis.

## Running the module

The module is normally run after `08b_GSE244983_malignant_subclustering_v1.4.2_FINAL.R` and before the spatial/CosMx modules. From the project root, run:

```r
source(
  "scripts/analysis/09A_GSE244983_doublet_CNV_sensitivity/run_09A_GSE244983_doublet_CNV_sensitivity.R",
  echo = TRUE
)
```

The runner executes the following steps in order:

1. `01_input_lineage_audit.R`
2. `02_pANN_doublet_detection.R`
3. `03_finalize_doublet_filter.R`
4. `04_primary_doublet_robustness.R`
5. `05_doublet_rate_sensitivity.R`
6. `06_generate_inferCNV_gene_order.R`
7. `07_run_inferCNV_support.R`
8. `08_final_module_audit.R`

## Doublet analysis

Doublets are assessed separately within each sample using a manual pANN implementation. The analysis uses `pN = 0.25`, `pK = 0.09`, PCA dimensions 1–20, 2,000 variable features, a maximum of 300 neighbors, and base seed `1234`. Homotypic adjustment is based on the current `SeuratCluster` annotation. The primary expected doublet rate is 7.5%, with 5% and 10% included as sensitivity settings.

The primary singlet object is created by subsetting only; it is not renormalized or reclustered. This is deliberate so that the sensitivity analysis tests cell removal rather than rebuilding the primary single-cell analysis from scratch.

## inferCNV support analysis

Gene order is derived from `TxDb.Hsapiens.UCSC.hg38.knownGene` and `org.Hs.eg.db`. T/NK, T/NK/Treg-like, B/Plasma and Myeloid cells form the reference pool; CAF/stromal, endothelial, malignant and cycling malignant cells are excluded from that reference. Sampling is capped at 750 reference cells per immune major cell type and 600 non-reference cells per current Seurat cluster.

inferCNV is run with cutoff `0.1`, `denoise = TRUE`, `HMM = FALSE` and `cluster_by_groups = TRUE`. The default thread count is up to eight and can be changed with `ICB_INFERCNV_THREADS`. These results are treated as exploratory CNV-like support, not as an independent malignant-cell classifier.

## Outputs

Tables are written to `results/tables/GSE244983/doublet_CNV_sensitivity/`, R objects to `results/objects/GSE244983/doublet_CNV_sensitivity/`, and inferCNV runtime files to `results/inferCNV/GSE244983_doublet_CNV_sensitivity/`. `D08_S7_source_manifest.csv` identifies the source tables used for the revised Supplementary Table S7.

The scripts do not install packages automatically or search recursively for alternative inputs. Required dependencies are checked before the time-consuming pANN/inferCNV steps begin.
