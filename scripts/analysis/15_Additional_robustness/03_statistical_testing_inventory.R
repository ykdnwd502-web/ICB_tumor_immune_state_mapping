############################################################
## 03_statistical_testing_inventory.R
##
## Reviewer-requested project-wide inventory of statistical tests,
## uncertainty procedures, multiplicity handling, and testing burden.
##
## IMPORTANT:
##   This is a reporting inventory of already-frozen analyses.
##   Frozen source scripts refer to the public-release execution names.
##   It does not perform or reinterpret any scientific analysis.
############################################################

options(stringsAsFactors=FALSE)
if (!requireNamespace("readr", quietly=TRUE)) stop("Package readr is required.",call.=FALSE)
PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR", unset="D:/ICB_resistance_project")
PROJECT_DIR <- normalizePath(PROJECT_DIR,winslash="/",mustWork=TRUE)
OUT_DIR <- file.path(PROJECT_DIR,"results","tables","additional_robustness","statistical_testing_inventory")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)

inventory <- data.frame(
  ID=sprintf("T%02d",1:20),
  Analysis_module=c(
    "Bulk state-score visualization",
    "Primary ICB response association",
    "GSE78220 per-state permutation fragility",
    "GSE78220 max-AUC state-selection permutation",
    "GSE78220 leave-one-sample-out sensitivity",
    "Exploratory survival association",
    "Compact composition-aware bulk sensitivity",
    "Formal composition-adjusted bulk ICB sensitivity grid",
    "Cross-cohort harmonized ICB comparability",
    "Single-cell source attribution",
    "Patient-level single-cell structure audit",
    "CAF/stromal support analysis",
    "Malignant-cell program heterogeneity",
    "Primary Visium dual-high enrichment",
    "Primary Visium continuous/spatial permutation analyses",
    "External melanoma spatial recurrence audits",
    "CosMx dual-high threshold sensitivity",
    "IF/mIF marker-state correlation audit",
    "CellChat-nominated CosMx ligand-receptor proximity",
    "ICBcomb perturbation prioritization"
  ),
  Dataset_resource=c(
    "GSE244982","GSE78220","GSE78220","GSE78220","GSE78220","GSE78220",
    "GSE244982/GSE78220","GSE78220/GSE91061","GSE78220/GSE91061","GSE244983","GSE244983","GSE244983","GSE244983",
    "Primary Visium melanoma section","Primary Visium melanoma section","GSE250636 Visium + Thrane 2018 legacy-ST","CosMx RNA-SMI","CosMx IF/mIF + RNA-state matches","GSE244983 CellChat + CosMx RNA-SMI","ICBcomb"
  ),
  Unit_of_analysis=c(
    "Bulk sample","Bulk sample","Bulk sample","Bulk sample","Bulk sample","Patient/bulk sample",
    "Bulk sample","Bulk sample","Bulk sample","Single cell","Patient summary derived from single cells","Single cell","Malignant single cell",
    "Visium spot","Visium spot / spatial graph","Spatial spot within each external dataset","CosMx cell","Matched CosMx cell","CosMx cell neighborhood / axis-contrast","ICBcomb-reported strategy signal"
  ),
  Statistical_test_or_procedure=c(
    "Descriptive standardized-score heatmap/PCA; no null-hypothesis test",
    "Fixed-direction ROC AUC for non-response; stratified bootstrap 95% CI",
    "Response-label permutation of AUC for each predefined state",
    "Permutation of the maximum AUC across the four predefined states",
    "LOSO recomputation of AUC after deleting one sample at a time",
    "Univariable Cox proportional-hazards model per standardized state score",
    "State-composition Spearman audit plus residualized target-state AUC sensitivity",
    "Raw/residualized AUC, LOSO, seeded bootstrap and response-label permutation across prespecified methods",
    "Common-gene cohort-only ComBat; fixed-direction AUC with DeLong CI; stratified bootstrap; response permutation within cohort",
    "Six prespecified rank-biserial compartment contrasts; 24 state-marker-module Spearman correlations",
    "Per-patient cell counts, four-state means and target-state distributions; descriptive only",
    "CAF/stromal marker-expression summaries and z-standardized core-score localization; descriptive support",
    "Program-score summaries and Spearman correlations among five malignant-cell programs; descriptive",
    "Fisher exact odds ratio for operational dual-high overlap; effect-size summary only",
    "Spearman/partial Spearman, symmetric bivariate Moran-type statistic, and spatial permutation sensitivity",
    "Primary-pair Spearman, descriptive Fisher OR and bivariate Moran-type recurrence summaries",
    "Fisher OR across pooled/per-FOV top 20%, 25%, 30% operational cutoffs",
    "Spearman correlation for each marker/intensity feature × state pair",
    "Permutation-based ligand-high/receptor-high neighborhood enrichment after CellChat nomination",
    "Database-derived median prioritization score and support count; no formal hypothesis test"
  ),
  Testing_family_or_endpoint=c(
    "Four predefined states visualized; no inferential family",
    "Four prespecified state scores vs binary non-response",
    "Four state-specific AUC null distributions",
    "One selection-aware family comprising the four state AUCs",
    "Four states × 27 deletions",
    "Four prespecified state-score Cox models",
    "Target-state composition sensitivity plus prespecified composition modules",
    "Prespecified cohort × state × adjustment grid",
    "Four predefined states in one harmonized 60-sample response space",
    "Prespecified cell-compartment and canonical marker-module audit",
    "Four contributing patients; four predefined state means",
    "Locked 21-gene CAF/stromal panel and core score",
    "Five predefined malignant-cell programs",
    "One prespecified myeloid-Treg vs tumor-dedifferentiation/stromal pair per locked score definition",
    "Same prespecified primary state pair across locked raw/residualized and graph definitions",
    "One prespecified primary state pair per external dataset and locked score definition",
    "2 cutoff modes × 3 fractions",
    "8 protein/intensity features × 4 predefined RNA-derived states",
    "12 nominated ligand-receptor axes × 5 contrasts, evaluated separately for each k",
    "Three selected resistance-relevant state queries summarized by strategy class"
  ),
  Number_of_hypotheses_or_recomputations=c(
    "0 formal hypotheses","4 AUC endpoints","4 permutation tests","1 max-statistic family over 4 AUCs","108 AUC recomputations (4 × 27); no new hypothesis tests","4 Cox models",
    "6 target-state AUC definitions in the compact sensitivity set; composition correlations reported separately","56 model-level rows in the frozen formal AUC grid","4 AUC endpoints",
    "6 prespecified contrasts + 24 descriptive Spearman correlations","0 formal hypotheses; 4 patient summaries","0 formal hypotheses; 21 locked markers","10 unique off-diagonal program-pair correlations (25 matrix entries including diagonal)",
    "1 primary pair per locked definition; OR treated descriptively","1 primary pair per locked spatial definition; sensitivity grid not treated as independent confirmatory hypotheses","1 primary pair per dataset/definition; recurrence audit only","6 Fisher summaries","32 correlations","60 tests per k; 180 across k=10/20/30","0 formal hypotheses"
  ),
  Resampling_or_permutation_N=c(
    "N/A","10,000 bootstrap resamples","10,000 permutations per state","10,000 max-AUC permutations","N/A","N/A",
    "5,000 bootstrap resamples for residualized AUC CIs","10,000 bootstrap + 10,000 permutation iterations per frozen model","2,000 stratified bootstrap resamples per state; 5,000 within-cohort permutations per state","N/A","N/A","N/A","N/A","N/A","999 permutations in locked spatial permutation analyses","N/A for descriptive recurrence summaries","N/A for Fisher threshold summaries; neighborhood sensitivity elsewhere uses 999 permutations","N/A","999 permutations per axis-contrast test","N/A"
  ),
  Multiple_testing_correction=c(
    "Not applicable","None; four states are prespecified and results are descriptive","Per-state empirical P values reported; selection-aware max-AUC permutation is reported separately","Max-statistic permutation accounts for selection among four states","Not applicable","None; exploratory/descriptive models",
    "No across-model multiplicity claim; prespecified sensitivity analysis","No global correction across heterogeneous sensitivity models; model-level permutation results are interpreted as robustness diagnostics","None; descriptive cross-cohort comparability sensitivity","Not applicable; effect-size/source-attribution audit","Not applicable","Not applicable","Not applicable","Not applicable to descriptive Fisher OR","Permutation-based spatial null; no claim that sensitivity definitions are independent tests","No global multiplicity claim across external recurrence readouts","None; threshold grid is sensitivity analysis","Benjamini-Hochberg FDR in the frozen authoritative 32-row correlation table","Benjamini-Hochberg FDR separately within each k (60 tests per k)","Not applicable"
  ),
  Adjusted_significance_threshold=c(
    "N/A","N/A — no confirmatory multiplicity-adjusted threshold","Empirical P ≤ 0.05 if referenced; interpreted with max-AUC audit","Empirical max-statistic P ≤ 0.05","N/A","N/A — exploratory",
    "N/A — sensitivity analysis","N/A — no global confirmatory threshold","N/A — descriptive sensitivity; permutation P reported","N/A","N/A","N/A","N/A","N/A — OR is descriptive","Permutation P ≤ 0.05 when used as spatial-null support; interpreted with Moran/sensitivity results","N/A — recurrence audit","N/A — threshold sensitivity","FDR < 0.05 where inferential flags are displayed; qualitative interpretation remains cautious","FDR < 0.05 within each k","N/A"
  ),
  Confidence_interval_or_uncertainty=c(
    "N/A","95% bootstrap CI","Empirical permutation distribution","Empirical max-AUC null distribution","LOSO range/median","95% Cox CI",
    "Bootstrap 95% AUC CI","Bootstrap 95% CI + permutation distribution + LOSO","DeLong 95% CI; bootstrap percentile 95% interval; exact permutation P","Rank-biserial effect size and Spearman rho","Across-patient descriptive values","Distributional summaries","Spearman rho matrix","Exact 95% CI shown for descriptive Fisher OR","Permutation null distributions; Moran-type statistic","Fisher 95% CI shown descriptively where available","Fisher exact 95% CI","Spearman rho, P and FDR","Empirical permutation P and BH FDR","Median score and support count"
  ),
  Inferential_status=c(
    "Descriptive","Descriptive association; not predictive validation","Fragility/sensitivity audit","Selection-bias sensitivity audit","Influence sensitivity audit","Exploratory/descriptive survival association",
    "Composition-aware sensitivity","Formalized robustness grid but not biomarker validation","Cross-cohort comparability sensitivity; not independent validation","Source-attribution/descriptive","Reviewer-requested descriptive structure audit","Supportive descriptive","Supportive descriptive","Descriptive effect-size summary","Spatial inferential robustness support","External recurrence audit only","Threshold sensitivity","Contextual association; mostly modest correlations","Spatial proximity-constrained hypothesis support only","Exploratory prioritization only"
  ),
  Figure_table_reference=c(
    "Figure 2; S2-S4","Figure 7; Table S10","Figure S11; Table S11","Figure S11; Table S11","Figure S11/S13","Figure 7","Figure S12","Figure S13; Table S13","Figure S30; Table S30","Figure 3; S7-S8; Tables S5-S6","Figure S29","Figure 4; S9; Table S8","Figures 5-6; S10; Table S9","Figure 8B","Figures S14-S19","Figure S20; Table S28-related recurrence audit","Figure S22","Figure S23; Table S23","Figures S24-S26; Tables S24-S26","Figure 9; S27-S28; Tables S27A/B"
  ),
  Frozen_source_contract=c(
    "01_GSE244982_bulk_state_scoring.R","02_GSE78220_state_scoring_response_survival.R","05_bulk_response_sensitivity_analysis.R","05_bulk_response_sensitivity_analysis.R","05_bulk_response_sensitivity_analysis.R","02_GSE78220_state_scoring_response_survival.R",
    "05_bulk_response_sensitivity_analysis.R","05_bulk_response_sensitivity_analysis.R","15/02_ICB_cross_cohort_comparability.R","07_GSE244983_scRNA_state_localization_source_attribution_v1.2.3_FINAL.R","15/01_GSE244983_patient_level_summary.R","08_GSE244983_CAF_stromal_analysis.R","08b_GSE244983_malignant_subclustering_v1.4.2_FINAL.R","12_Primary_Visium canonical source","12_Primary_Visium canonical source","10_S20 external spatial recurrence","11_CosMx","13_IF_mIF","14_CellChat_LR","ICBcomb frozen reporting lineage"
  ),
  stringsAsFactors=FALSE
)

readme <- data.frame(
  Item=c("Purpose","Multiplicity principle","Interpretation","Scope boundary"),
  Value=c(
    "Reviewer-requested summary of statistical procedures and testing burden across the manuscript.",
    "Correction is applied within scientifically coherent testing families when defined; heterogeneous sensitivity analyses are not pooled into a single artificial global family.",
    "N/A in the adjusted-threshold column means that the analysis is descriptive, an effect-size/sensitivity audit, or has no confirmatory multiplicity-adjusted threshold.",
    "This inventory does not change any frozen statistical method, endpoint, state definition, numerical result, or significance decision."
  ),
  stringsAsFactors=FALSE
)

readr::write_csv(inventory,file.path(OUT_DIR,"Table_S29_statistical_testing_inventory.csv"),na="")
readr::write_csv(readme,file.path(OUT_DIR,"Table_S29_README.csv"),na="")

required_cols <- c("Statistical_test_or_procedure","Number_of_hypotheses_or_recomputations","Multiple_testing_correction","Adjusted_significance_threshold")
pass <- nrow(inventory)==20L && all(required_cols %in% colnames(inventory)) && all(nzchar(inventory$Statistical_test_or_procedure)) && all(nzchar(inventory$Number_of_hypotheses_or_recomputations)) && all(nzchar(inventory$Multiple_testing_correction)) && all(nzchar(inventory$Adjusted_significance_threshold))
gate <- data.frame(Gate="Table S29 inventory completeness",Observed=paste0(nrow(inventory)," rows"),Expected="20 rows with explicit test/burden/correction/threshold fields",Pass=pass,stringsAsFactors=FALSE)
readr::write_csv(gate,file.path(OUT_DIR,"Table_S29_inventory_gate.csv"))
if (!pass) stop("Table S29 inventory completeness gate failed.",call.=FALSE)
cat("Table S29 statistical-testing inventory: PASS (20 rows)\n")
