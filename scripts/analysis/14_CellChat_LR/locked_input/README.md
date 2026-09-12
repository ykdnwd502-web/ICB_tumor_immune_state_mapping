# Historical CellChat cell-selection manifest

`CellChat_historical_sampling_group_manifest.csv` records the exact 3,219 GSE244983 cells and CellChat groups used in the historical Step16 analysis. It is kept with the public workflow so that the CellChat results can be reproduced without distributing the earlier ~2.5 GB Seurat object from which the selection was originally made.

The historical sampling used seed `20260504`, a maximum of 800 cells per group and a minimum of 30 cells per group after broad-cell-type standardization. All 3,219 cell IDs were confirmed to be present in the rebuilt public GSE244983 object. Running CellChat on that object with this manifest reproduced the frozen S24 summary across all audited fields (12/12 rows per field; maximum numeric difference < `4e-16`).

SHA256 for the manifest:

```text
b4e1dd0cd5f3da4808d9284e23de7491dd4779186021530971c7dad43f579620
```

The manifest should therefore be treated as part of the historical reproducibility record and should not be regenerated or edited.
