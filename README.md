# ICB tumor immune state mapping

This repository contains the reproducible analysis workflow, frozen input definitions, and release metadata for the ICB tumor immune state mapping project.

## Repository organization

```
ICB_tumor_immune_state_mapping/
├── ICBcomb/
├── scripts/
├── public_release/
├── release_metadata/
└── README.md
```

## Reproduction workflow

The recommended workflow is:

1. Restore frozen inputs.
2. Verify input integrity using SHA256 manifests.
3. Install the required computational environment.
4. Execute sequential workflows.
5. Compare regenerated outputs with frozen references and release metadata.

The workflow manifest is:

```
release_metadata/RUN_MANIFEST.csv
```

Execution records are summarized in:

```
release_metadata/REPRODUCTION_LOG_INDEX.csv
```

## Frozen inputs

Frozen input definitions are provided in:

```
ICBcomb/inputs/ICBcomb_INPUT_FREEZE_v1.0/
```

Integrity validation files include:

```
SHA256SUMS.csv
FREEZE_STATUS.txt
validate_ICBcomb_INPUT_FREEZE_v1.0.R
```

## Computational environment

Environment information is provided in:

```
release_metadata/INSTALLED_PACKAGES_BEFORE_RUN.csv
release_metadata/sessionInfo/
```

## Release validation

The repository includes script, metadata, and repository SHA256 manifests together with execution logs for reproducible verification.

## Data availability

Source code and reproducibility metadata are available through GitHub. Large frozen archives are deposited through Zenodo.

## Citation

Please cite the associated publication and archived release version when using this repository.
