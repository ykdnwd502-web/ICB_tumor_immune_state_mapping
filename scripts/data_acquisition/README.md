# Raw-data acquisition

This folder contains the data-acquisition layer used to reconstruct the public/raw inputs expected by the released analysis workflow.

It does **not** run the scientific analyses or change any of the frozen analytical definitions, thresholds, tumor–immune state definitions, statistical procedures, or interpretation rules. Its only job is to turn publicly available source data into the directory structure expected by the downstream scripts.

## Quick start

Place this folder at:

```text
D:/ICB_resistance_project/scripts/data_acquisition/
```

The default project root is:

```text
D:/ICB_resistance_project
```

If the project lives somewhere else, set the root before running the scripts:

```r
Sys.setenv(
  ICB_PROJECT_DIR = "D:/your_project_path/ICB_resistance_project"
)
```

Then run:

```r
source(
  "D:/ICB_resistance_project/scripts/data_acquisition/00_run_data_acquisition.R",
  echo = TRUE
)
```

The runner calls, in order:

```text
01_download_GEO_inputs.R
02_prepare_primary_Visium_10x.R
03_prepare_CosMx_Dryad_Slide4.R
04_prepare_external_spatial_inputs.R
05_validate_raw_data_contract.R
```

The scripts are intended to be rerunnable. Files that are already present and pass the expected checks are reused rather than downloaded again, and existing source data are not deleted.

One exception is the CosMx Dryad archive: depending on the current Dryad access mechanism, `Slide_4.zip` may need to be downloaded manually. In that case the runner can stop at Step 03 with `MANUAL_DOWNLOAD_REQUIRED`; this is a data-access condition, not a failure of the scientific analysis.

## Public datasets

### GEO datasets

The GEO acquisition step prepares the following inputs.

**GSE244982**

```text
data_raw/GSE244982/
└── GSE244982_ProcessedData_bulkRNAseq.txt.gz
```

This is the bulk RNA-seq expression input used by the discovery state-scoring analysis.

**GSE244983**

```text
data_raw/GSE244983/
├── GSE244983_RawCounts_scRNAseq.txt.gz
├── GSE244983_SingleCellAnnotations.txt.gz
└── GSE244983_NormalizedData_scRNAseq.txt.gz
```

The canonical single-cell reconstruction uses the raw-count matrix together with the author-provided cell annotation table. The normalized matrix is retained for source completeness but is not the primary input used to rebuild the canonical GSE244983 Seurat object.

**GSE78220**

```text
data_raw/external_validation/GSE78220/
└── GSE78220_PatientFPKM.xlsx
```

The script also rebuilds:

```text
GSE78220_pData_raw_from_GEO.csv
```

through `GEOquery`.

**GSE91061**

```text
data_raw/external_validation/GSE91061/
├── GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz
├── GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz
├── GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz
└── GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt.gz
```

The frozen external-response analysis uses the RLD-transformed expression data. The phenotype table is rebuilt as:

```text
GSE91061_pData_raw_from_GEO.csv
```

again using `GEOquery`.

The GEO supplementary `.gz` files are intentionally left compressed because the validated analysis scripts read them directly.

## Primary Visium dataset

The primary spatial dataset is the 10x Genomics:

```text
CytAssist_FFPE_Human_Skin_Melanoma
```

The acquisition script downloads the filtered feature-barcode matrix archive, the spatial archive, and the metrics summary into:

```text
data_raw/spatial_10x_human_melanoma_IF_FFPE/
```

The downstream Visium builder expects the matrix files in both of the following locations:

```text
filtered_feature_bc_matrix/
```

and:

```text
filtered_feature_bc_matrix_extracted/filtered_feature_bc_matrix/
```

This duplication is intentional. The validated Visium workflow consumes the first location while retaining the second as the historical audit location.

A typical expected layout is:

```text
data_raw/spatial_10x_human_melanoma_IF_FFPE/
├── filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── filtered_feature_bc_matrix_extracted/
│   └── filtered_feature_bc_matrix/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
└── spatial/
    ├── tissue_positions.csv
    ├── scalefactors_json.json
    └── ...
```

The processed Visium Seurat object is **not** created here. It is generated later by the Primary Visium analysis workflow.

## CosMx Slide 4

The CosMx source is the Dryad dataset:

```text
10.5061/dryad.ksn02v7b1
```

Only `Slide_4.zip` is required by the released CosMx workflow.

Automated Dryad downloads are not always available, so the most reliable route is to download `Slide_4.zip` from the DOI landing page and place it at:

```text
data_raw/CosMx_melanocytic_tumors_Dryad/Slide_4.zip
```

Then rerun the acquisition runner.

If a current direct download URL is available, it can be supplied through:

```r
Sys.setenv(
  COSMX_SLIDE4_URL = "CURRENT_DIRECT_DOWNLOAD_URL"
)
```

After extraction, the CosMx and IF/mIF modules expect:

```text
data_raw/CosMx_melanocytic_tumors_Dryad/
└── Slide_4/
    └── Run5611_MK3/
        ├── Run5611_MK3_exprMat_file.csv
        ├── Run5611_MK3_metadata_file.csv
        ├── Run5611_MK3_fov_positions_file.csv
        └── Run5611_MK3_tx_file.csv
```

The historical `Slide_4_unzipped/` directory is not recreated because it was a redundant copy and is not needed by the released analysis chain.

## External spatial recurrence datasets

### GSE250636

The script downloads:

```text
GSE250636_RAW.tar
```

and extracts it under:

```text
data_external_spatial_18J/
└── GSE250636/
    └── unpacked_recursive/
        └── GSE250636_RAW/
```

The final validator expects nine complete sample-specific Visium inputs. For each retained sample, this includes:

```text
*_filtered_feature_bc_matrix.h5
*_tissue_positions_list.csv.gz
```

The nine sample prefixes used by the frozen workflow are:

```text
GSM7983358_sample2
GSM7983359_sample4
GSM7983360_sample7
GSM7983361_sample8
GSM7983362_sample12
GSM7983363_sample13
GSM7983364_sample14
GSM7983365_sample15
GSM7983366_sample16
```

### Thrane et al. 2018

The legacy spatial-transcriptomics source is downloaded from Zenodo as:

```text
Thrane_et_al_2018_CAN_RES.zip
```

and stored under:

```text
data_external_spatial_18J/
└── Thrane2018_legacyST/
    └── raw_zenodo/
```

The archive is extracted to:

```text
data_external_spatial_18J/
└── Thrane2018_legacyST/
    └── unpacked/
        └── Thrane_et_al_2018_CAN_RES/
```

The public recurrence workflow uses eight count matrices:

```text
ST_mel1_rep1_counts.tsv
ST_mel1_rep2_counts.tsv
ST_mel2_rep1_counts.tsv
ST_mel2_rep2_counts.tsv
ST_mel3_rep1_counts.tsv
ST_mel3_rep2_counts.tsv
ST_mel4_rep1_counts.tsv
ST_mel4_rep2_counts.tsv
```

Objects under `standardized/*.rds` are not downloaded because they are regenerated by the public analysis workflow.

## Frozen external-query inputs

### ICBcomb

ICBcomb is deliberately **not** re-queried during reproduction.

The historical web-platform results used in the study are already included in the prepared reproduction workspace under:

```text
ICBcomb/inputs/ICBcomb_INPUT_FREEZE_v1.0/
```

The acquisition scripts therefore do not download or regenerate these files. A separate `Frozen_Input_Archive_v1.0` is retained as an archival and integrity-preservation copy of historically fixed inputs, but restoring that archive is not an additional step when the prepared reproduction workspace is used.

This preserves the exact external-query state used in the study and avoids drift caused by later platform or database updates.

## Data not used by the released workflow

Two development-era sources are intentionally excluded from the current acquisition layer:

```text
data_raw/ICBcomb_reversal/
```

which is a legacy duplicate of the frozen ICBcomb query material, and:

```text
Gide2019 / PRJEB23709
```

which is not part of the current frozen manuscript/data-availability chain.

Neither is required to reproduce the released analysis.

## Validation

The final acquisition step writes its audit and provenance files under:

```text
release_metadata/raw_data/
```

The main records are:

```text
RAW_DATA_CONTRACT.csv
RAW_DATA_CONTRACT_SUMMARY.txt
RAW_DATA_ACQUISITION_RUN_LOG.csv
sessionInfo_data_acquisition.txt
```

Dataset-specific checks may also produce:

```text
GEO_DOWNLOAD_MANIFEST.csv
VISIUM_10X_RAW_CONTRACT.csv
COSMX_SLIDE4_RAW_CONTRACT.csv
COSMX_SLIDE4_ACQUISITION_PROVENANCE.csv
GSE250636_RAW_CONTRACT.csv
THRANE2018_RAW_CONTRACT.csv
```

A normal acquisition run may end in one of several states:

```text
PASS
MANUAL_DOWNLOAD_REQUIRED
BLOCKED_RAW_DATA_MISSING
READY_FOR_ANALYSIS_CHAIN
```

`MANUAL_DOWNLOAD_REQUIRED` means that a repository requires a browser download before the workflow can continue. It should not be interpreted as a failure of the analysis code.

The analysis chain should begin only after the required raw-data contract has been satisfied. The expected final state is:

```text
Decision: READY_FOR_ANALYSIS_CHAIN
```

## How this fits into the full release

The acquisition layer sits at the beginning of the public reproduction workflow:

```text
Public source repositories
        ↓
scripts/data_acquisition/
        ↓
validated raw-data directory structure
        ↓
prepared reproduction workspace
(frozen inputs and configuration already in place)
        ↓
analysis scripts and module runners
        ↓
regenerated figures, tables and processed outputs
        ↓
reproduction checks and release metadata
```

The released reproduction workspace already contains the historically fixed inputs and configuration files at the locations expected by the analysis scripts. No separate frozen-input restoration step is required when that prepared workspace is used.

`Frozen_Input_Archive_v1.0` is retained separately as an archival and integrity-preservation copy. It is useful for backup, verification, or reconstruction of the workspace if needed, but it is not part of the normal execution sequence.

External repositories may change URLs, redirects, authentication requirements, or API policies over time. For that reason, the acquisition layer relies on stable accessions and DOIs wherever possible and validates the local filenames and directory structure expected by the downstream analysis.

In practice, reproducibility depends on using the same public source datasets, the same frozen inputs and configuration, and the same released analysis scripts—not on any particular transient download URL.

## Public sources used by the project

The released workflow draws on:

- GSE244982
- GSE244983
- GSE78220
- GSE91061
- 10x Genomics CytAssist Human Melanoma IF-stained FFPE Visium dataset
- Dryad CosMx melanoma RNA-SMI dataset (DOI: `10.5061/dryad.ksn02v7b1`)
- GSE250636
- Thrane et al. 2018 legacy spatial-transcriptomics dataset
- frozen historical ICBcomb web-platform query outputs
