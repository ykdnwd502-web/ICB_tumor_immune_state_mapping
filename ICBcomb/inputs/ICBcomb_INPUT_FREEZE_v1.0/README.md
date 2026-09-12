# ICBcomb input freeze v1.0

## Status

**Decision: READY_FOR_ICBCOMB_REBUILD**

Blocking failures: **0**

This package freezes the exact ICBcomb input boundary used by the historical analysis.
It does not live-requery the ICBcomb platform and does not reinterpret the ICBcomb
results as therapeutic validation.

## Final perturbation-query contract

- immune-defective/cold: RESTORE hybrid, 98 genes
- myeloid–Treg immunosuppressive: SUPPRESS hybrid, 73 genes
- tumor-dedifferentiation/stromal-remodeling: SUPPRESS hybrid, 81 genes
- melanocytic differentiation: hybrid reference/control only, 32 genes

The three resistance-state query TXT files are gene-for-gene identical to their
corresponding entries in `final_ICBcomb_gene_sets_CLEAN.gmt`.

## Platform result boundary

Exactly three raw platform CSVs are frozen. Each has 15 rows and the required
14-column ICBcomb result schema. Historical reconstruction yields:

- combined raw rows: 45
- positive prioritization rows (`NES > 0`): 42
- state × strategy summary rows: 10

The historical audit contract records the platform access date as **2025-04-30**.
This package preserves that contract verbatim; filesystem timestamps and text embedded
inside platform-generated path strings are not treated as independent access-date authority.

## Directory roles

- `data/raw_platform_results/`: immutable raw platform downloads
- `data/query_gene_set_bundle/`: exact frozen query/config/provenance snapshot
- `data/historical_reference/`: read-only numerical/reporting reference tables
- `audit/`: SHA256 inventories and freeze gates

The nested `state_associated_genes_set.zip` inside the query bundle is retained only
because it was present in the frozen working-directory snapshot. It is redundant for
public execution and may be omitted later from the streamlined final public release,
provided its extracted provenance files are retained.

## Interpretation boundary

ICBcomb results are exploratory database-derived perturbation-prioritization signals.
They are not experimental validation, therapeutic efficacy evidence, clinical
actionability, or causal proof.

## Next step

Run the ICBcomb source-locked rebuild/audit chain against this frozen input boundary,
then freeze Figure 9, Supplementary Figures S27/S28, and Supplementary Tables S27A/B.

## Historical-table compatibility note

The historical combined/positive CSV tables reproduce the raw-source values exactly up to ordinary CSV numeric round-trip precision (maximum observed difference < 1e-12); all nonnumeric source fields and row identities are exact. This 1e-12 threshold is an audit-serialization tolerance only, not a scientific/statistical tolerance.
