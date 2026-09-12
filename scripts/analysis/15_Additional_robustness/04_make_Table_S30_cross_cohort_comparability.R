############################################################
## 04_make_Table_S30_cross_cohort_comparability.R
##
## Reporting-only reconstruction of Supplementary Table S30 from
## canonical cross-cohort and S30 reporting sources.
## Historical 18H tables are intentionally not read.
############################################################

options(stringsAsFactors=FALSE)
required <- c("readr","dplyr","openxlsx")
miss <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)>0) stop("Missing package(s): ",paste(miss,collapse=", "),call.=FALSE)

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR",unset="D:/ICB_resistance_project")
PROJECT_DIR <- normalizePath(PROJECT_DIR,winslash="/",mustWork=TRUE)
CROSS_DIR <- file.path(PROJECT_DIR,"results","tables","cross_cohort_harmonization")
S30_DIR <- file.path(PROJECT_DIR,"results","reporting","additional_robustness","S30_cross_cohort","tables")
OUT_DIR <- file.path(PROJECT_DIR,"results","tables")
AUDIT_DIR <- file.path(PROJECT_DIR,"results","tables","additional_robustness","S30_table")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
dir.create(AUDIT_DIR,recursive=TRUE,showWarnings=FALSE)

readc <- function(path) {
  if(!file.exists(path)) stop("Missing file: ",path,call.=FALSE)
  readr::read_csv(path,show_col_types=FALSE,progress=FALSE) |> as.data.frame(check.names=FALSE)
}

panelA <- readc(file.path(S30_DIR,"S30_PanelA_raw_AUC.csv"))
panelC <- readc(file.path(S30_DIR,"S30_PanelC_composition_adjusted.csv"))
panelD <- readc(file.path(S30_DIR,"S30_PanelD_cross_cohort.csv"))
report_gates <- readc(file.path(S30_DIR,"S30_reporting_gates.csv"))
boot <- readc(file.path(CROSS_DIR,"cross_cohort_seeded_bootstrap_AUC_summary.csv"))
perm <- readc(file.path(CROSS_DIR,"cross_cohort_seeded_within_cohort_permutation_summary.csv"))
method <- readc(file.path(CROSS_DIR,"cross_cohort_method_manifest.csv"))
cross_gates <- readc(file.path(CROSS_DIR,"cross_cohort_reproducibility_gates.csv"))

if(!all(as.logical(report_gates$Pass))) stop("S30 reporting gates are not all PASS.",call.=FALSE)
if(!all(as.logical(cross_gates$Pass))) stop("Cross-cohort producer gates are not all PASS.",call.=FALSE)

state_order <- c("Immune_defective_Cold","Myeloid_Treg_Immunosuppressive","Tumor_dedifferentiation_Stromal_remodeling","Melanocytic_Differentiation")
requiredA <- c("State","StateDisplay","AUC","N","Responders","NonResponders")
requiredD <- c("State","StateDisplay","N","Responders","NonResponders","AUC","CI_low","CI_high","Permutation_P","Permutation_N","N_extreme_exact")
requiredB <- c("State","Mean_AUC","Median_AUC","Percentile_low","Percentile_high","N_boot","BootstrapStrata")
stopifnot(all(requiredA %in% colnames(panelA)),all(requiredD %in% colnames(panelD)),all(requiredB %in% colnames(boot)))

compact <- panelD |>
  dplyr::select(State,StateDisplay,Harmonized_N=N,Harmonized_Responders=Responders,Harmonized_NonResponders=NonResponders,Harmonized_AUC=AUC,DeLong_CI_low=CI_low,DeLong_CI_high=CI_high,Permutation_P,Permutation_N,N_extreme_exact) |>
  dplyr::left_join(panelA |> dplyr::select(State,Discovery_AUC=AUC,Discovery_N=N,Discovery_Responders=Responders,Discovery_NonResponders=NonResponders),by="State") |>
  dplyr::left_join(boot |> dplyr::select(State,Bootstrap_mean=Mean_AUC,Bootstrap_median=Median_AUC,Bootstrap_CI_low=Percentile_low,Bootstrap_CI_high=Percentile_high,Bootstrap_N=N_boot,Bootstrap_strata=BootstrapStrata),by="State") |>
  dplyr::mutate(
    AUC_change=Harmonized_AUC-Discovery_AUC,
    Interpretation=dplyr::case_when(
      State=="Tumor_dedifferentiation_Stromal_remodeling" ~ "Discovery association attenuated after common-space harmonization; harmonized AUC remained modest and the exact within-cohort permutation P was not below 0.05. This is not cross-cohort predictive validation.",
      TRUE ~ "Cross-cohort harmonized estimate reported as descriptive comparability/sensitivity evidence; no predictive-validation claim."
    )
  )
compact <- compact[match(state_order,compact$State),,drop=FALSE]

TARGET <- "Tumor_dedifferentiation_Stromal_remodeling"
target <- compact[compact$State==TARGET,,drop=FALSE]
stopifnot(nrow(target)==1L)
expected_auc <- 0.6057142857142858
expected_p <- 0.1569686062787443
expected_boot <- 0.606762857142857
expected_boot_low <- 0.468571428571429
expected_boot_high <- 0.7406
TOL <- 1e-12

gates <- data.frame(
  Gate=c("Four states","Combined N=60","Target harmonized AUC","Target exact permutation P","Target bootstrap mean","Target bootstrap lower","Target bootstrap upper","Historical 18H excluded"),
  Observed=c(nrow(compact),paste(sort(unique(compact$Harmonized_N)),collapse=";"),target$Harmonized_AUC,target$Permutation_P,target$Bootstrap_mean,target$Bootstrap_CI_low,target$Bootstrap_CI_high,"TRUE"),
  Expected=c(4,"60",expected_auc,expected_p,expected_boot,expected_boot_low,expected_boot_high,"TRUE"),
  Pass=c(nrow(compact)==4L,all(compact$Harmonized_N==60L),abs(target$Harmonized_AUC-expected_auc)<=TOL,abs(target$Permutation_P-expected_p)<=TOL,abs(target$Bootstrap_mean-expected_boot)<=TOL,abs(target$Bootstrap_CI_low-expected_boot_low)<=TOL,abs(target$Bootstrap_CI_high-expected_boot_high)<=TOL,TRUE),
  stringsAsFactors=FALSE
)
readr::write_csv(gates,file.path(AUDIT_DIR,"Table_S30_reporting_gates.csv"))
readr::write_csv(compact,file.path(AUDIT_DIR,"Table_S30_compact_source.csv"))
if(any(!gates$Pass)) stop("Table S30 canonical gate failure: ",paste(gates$Gate[!gates$Pass],collapse=" | "),call.=FALSE)

readme <- data.frame(
  Item=c("Table","Purpose","GSE78220 universe","GSE91061 universe","Combined universe","Common genes","Harmonization","Outcome","Direction","Uncertainty","Bootstrap","Permutation","Multiplicity","Interpretive boundary"),
  Value=c(
    "Supplementary Table S30","Quantitative cross-cohort comparability analysis requested during revision","27 strict responders/non-responders","33 pretreatment strict CR/PR-vs-PD samples","60 samples","19,799 exact common HGNC symbols","Cohort-only ComBat; response status not included; gene-wise z-standardization after ComBat","Non-response","Higher state score -> greater non-response tendency","DeLong 95% CI","2,000 seeded resamples stratified by cohort × response","5,000 seeded response-label permutations within cohort; exact discrete two-sided distance from AUC=0.5","No multiplicity-adjusted confirmatory threshold; four states were prespecified and this is a descriptive sensitivity analysis","Not independent validation, not a calibrated predictive model, and not evidence of cohort-invariant biology"
  ),stringsAsFactors=FALSE
)

wb <- openxlsx::createWorkbook(creator="15_Additional_robustness")
for(nm in c("README","Table S30","Composition sensitivity","Bootstrap","Permutation","Method contract","Reporting gates")) openxlsx::addWorksheet(wb,nm)
openxlsx::writeData(wb,"README",readme)
openxlsx::writeData(wb,"Table S30",compact)
openxlsx::writeData(wb,"Composition sensitivity",panelC)
openxlsx::writeData(wb,"Bootstrap",boot)
openxlsx::writeData(wb,"Permutation",perm)
openxlsx::writeData(wb,"Method contract",method)
openxlsx::writeData(wb,"Reporting gates",gates)
header <- openxlsx::createStyle(textDecoration="bold",fgFill="#D9EAF7",border="Bottom",halign="center",valign="center",wrapText=TRUE)
body <- openxlsx::createStyle(valign="top",wrapText=TRUE)
sheet_data <- list(
  "README"=readme,
  "Table S30"=compact,
  "Composition sensitivity"=panelC,
  "Bootstrap"=boot,
  "Permutation"=perm,
  "Method contract"=method,
  "Reporting gates"=gates
)
for(nm in names(sheet_data)) {
  dat <- sheet_data[[nm]]
  nr <- nrow(dat)
  nc <- ncol(dat)
  openxlsx::addStyle(wb,nm,header,rows=1,cols=seq_len(nc),gridExpand=TRUE)
  if(nr>0) openxlsx::addStyle(wb,nm,body,rows=2:(nr+1),cols=seq_len(nc),gridExpand=TRUE)
  openxlsx::freezePane(wb,nm,firstRow=TRUE)
  openxlsx::setColWidths(wb,nm,cols=seq_len(nc),widths="auto")
}
out <- file.path(OUT_DIR,"Supplementary Table S30. Quantitative cross-cohort comparability analysis of predefined tumor-immune state associations between GSE78220 and GSE91061.xlsx")
openxlsx::saveWorkbook(wb,out,overwrite=TRUE)
cat("Table S30 canonical rebuild PASS\n")
cat("Target raw AUC: ",format(target$Discovery_AUC,digits=15),"\n",sep="")
cat("Target harmonized AUC: ",format(target$Harmonized_AUC,digits=15),"\n",sep="")
cat("Target permutation P: ",format(target$Permutation_P,digits=15),"\n",sep="")
cat("Saved: ",normalizePath(out,winslash="/",mustWork=FALSE),"\n",sep="")
