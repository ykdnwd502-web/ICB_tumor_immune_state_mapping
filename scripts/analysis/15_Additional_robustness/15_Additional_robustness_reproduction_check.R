############################################################
## 15_Additional_robustness_reproduction_check.R
## Final closure audit for Reviewer-requested additional robustness.
############################################################
options(stringsAsFactors=FALSE)
required <- c("readr","openxlsx")
miss <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)>0) stop("Missing package(s): ",paste(miss,collapse=", "),call.=FALSE)
PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR",unset="D:/ICB_resistance_project")
PROJECT_DIR <- normalizePath(PROJECT_DIR,winslash="/",mustWork=TRUE)
AUDIT_DIR <- file.path(PROJECT_DIR,"results","audit","15_Additional_robustness_reproduction_check")
dir.create(AUDIT_DIR,recursive=TRUE,showWarnings=FALSE)

readc <- function(path){if(!file.exists(path)) stop("Missing file: ",path,call.=FALSE); readr::read_csv(path,show_col_types=FALSE,progress=FALSE)|>as.data.frame(check.names=FALSE)}
TOL <- 1e-12

s29_src <- file.path(PROJECT_DIR,"results","tables","additional_robustness","S29_patient_level")
s29_rep <- file.path(PROJECT_DIR,"results","reporting","additional_robustness","S29_patient_level")
s30_cross <- file.path(PROJECT_DIR,"results","tables","cross_cohort_harmonization")
s30_rep <- file.path(PROJECT_DIR,"results","reporting","additional_robustness","S30_cross_cohort","tables")
s29_inv <- file.path(PROJECT_DIR,"results","tables","additional_robustness","statistical_testing_inventory")
s30_tab_audit <- file.path(PROJECT_DIR,"results","tables","additional_robustness","S30_table")
fig_dir <- file.path(PROJECT_DIR,"results","figures","supplementary")
s30_fig_dir <- fig_dir

s29_g <- readc(file.path(s29_src,"S29_lineage_gates.csv"))
s29_counts <- readc(file.path(s29_src,"S29_patient_counts.csv"))
s29_means <- readc(file.path(s29_src,"S29_patient_state_means.csv"))
s29_rg <- readc(file.path(s29_rep,"S29_reporting_gates.csv"))

expected_counts <- c(Pat_ICBnaive2=9081L,Pat_ICBnaive1=7414L,Pat42=7188L,Pat5=2370L)
expected_means <- rbind(
  Pat_ICBnaive2=c(0.3325608419486311,-0.006194506141643475,0.1759250596840239,0.29067310298207455),
  Pat_ICBnaive1=c(0.39356295782044515,-0.2730601224886434,-0.22499937648125515,0.16393561070755508),
  Pat42=c(-0.2653437647613925,-0.19702662352376465,0.026172049120634644,0.037924389913327694),
  Pat5=c(-1.7006623603001734,1.4755052440050156,-0.04960117250536579,-1.7416124812924185)
)
mean_cols <- c("Immune_defective_Cold_mean","Myeloid_Treg_Immunosuppressive_mean","Tumor_dedifferentiation_Stromal_remodeling_mean","Melanocytic_Differentiation_mean")
s29_counts <- s29_counts[match(names(expected_counts),s29_counts$Sample),,drop=FALSE]
s29_means <- s29_means[match(rownames(expected_means),s29_means$Sample),,drop=FALSE]
s29_mean_diff <- max(abs(as.matrix(s29_means[,mean_cols])-expected_means),na.rm=TRUE)

cross_g <- readc(file.path(s30_cross,"cross_cohort_reproducibility_gates.csv"))
s30_g <- readc(file.path(s30_rep,"S30_reporting_gates.csv"))
s30_d <- readc(file.path(s30_rep,"S30_PanelD_cross_cohort.csv"))
t30_g <- readc(file.path(s30_tab_audit,"Table_S30_reporting_gates.csv"))
t29_g <- readc(file.path(s29_inv,"Table_S29_inventory_gate.csv"))
t29 <- readc(file.path(s29_inv,"Table_S29_statistical_testing_inventory.csv"))

target <- s30_d[s30_d$State=="Tumor_dedifferentiation_Stromal_remodeling",,drop=FALSE]

expected_files <- c(
  file.path(fig_dir,"Supplementary_Figure_S29_patient_level_structure.pdf"),
  file.path(fig_dir,"Supplementary_Figure_S29_patient_level_structure.png"),
  file.path(s30_fig_dir,"Supplementary_Figure_S30_integrated_fragility_crosscohort.pdf"),
  file.path(s30_fig_dir,"Supplementary_Figure_S30_integrated_fragility_crosscohort.png"),
  file.path(PROJECT_DIR,"results","tables","Supplementary Table S29. Statistical analysis and multiple-testing summary.xlsx"),
  file.path(PROJECT_DIR,"results","tables","Supplementary Table S30. Quantitative cross-cohort comparability analysis of predefined tumor-immune state associations between GSE78220 and GSE91061.xlsx")
)

checks <- data.frame(
  Section=c(
    "S29 primary-lineage gates","S29 exact patient counts","S29 frozen patient-state means","S29 reporting gates",
    "Cross-cohort producer gates","S30 reporting gates","S30 target harmonized AUC","S30 target permutation P",
    "Table S29 inventory gate","Table S29 has 20 rows","Table S30 reporting gates","Final S29/S30/Table S29/Table S30 artifacts exist"
  ),
  Observed=c(
    paste0(sum(as.logical(s29_g$Pass)),"/",nrow(s29_g)),
    paste(s29_counts$n_cells,collapse=";"),
    format(s29_mean_diff,scientific=TRUE),
    paste0(sum(as.logical(s29_rg$Pass)),"/",nrow(s29_rg)),
    paste0(sum(as.logical(cross_g$Pass)),"/",nrow(cross_g)),
    paste0(sum(as.logical(s30_g$Pass)),"/",nrow(s30_g)),
    if(nrow(target)==1) format(target$AUC,digits=15) else "not unique",
    if(nrow(target)==1) format(target$Permutation_P,digits=15) else "not unique",
    paste0(sum(as.logical(t29_g$Pass)),"/",nrow(t29_g)),
    nrow(t29),
    paste0(sum(as.logical(t30_g$Pass)),"/",nrow(t30_g)),
    paste0(sum(file.exists(expected_files)),"/",length(expected_files))
  ),
  Expected=c(
    "all PASS","9081;7414;7188;2370","<=1e-12","all PASS","all PASS","all PASS","0.605714285714286","0.156968606278744","all PASS","20","all PASS",paste0(length(expected_files),"/",length(expected_files))
  ),
  Pass=c(
    all(as.logical(s29_g$Pass)),
    identical(as.integer(s29_counts$n_cells),as.integer(expected_counts)),
    is.finite(s29_mean_diff)&&s29_mean_diff<=TOL,
    all(as.logical(s29_rg$Pass)),
    all(as.logical(cross_g$Pass)),
    all(as.logical(s30_g$Pass)),
    nrow(target)==1L&&abs(target$AUC-0.6057142857142858)<=TOL,
    nrow(target)==1L&&abs(target$Permutation_P-0.1569686062787443)<=TOL,
    all(as.logical(t29_g$Pass)),
    nrow(t29)==20L,
    all(as.logical(t30_g$Pass)),
    all(file.exists(expected_files))
  ),
  stringsAsFactors=FALSE
)
readr::write_csv(checks,file.path(AUDIT_DIR,"15_final_gate_summary.csv"))
blocking <- sum(!checks$Pass)
status <- data.frame(Status=ifelse(blocking==0,"READY_FOR_PUBLIC_FREEZE","BLOCKED"),Blocking_failures=blocking,stringsAsFactors=FALSE)
readr::write_csv(status,file.path(AUDIT_DIR,"15_freeze_status.csv"))

cat("\n============================================================\n")
cat("15 ADDITIONAL ROBUSTNESS REPRODUCTION CHECK\n")
cat("S29 lineage: ",sum(as.logical(s29_g$Pass)),"/",nrow(s29_g),"\n",sep="")
cat("S29 max patient-state mean diff: ",format(s29_mean_diff,scientific=TRUE),"\n",sep="")
cat("Cross-cohort producer: ",sum(as.logical(cross_g$Pass)),"/",nrow(cross_g),"\n",sep="")
cat("S30 reporting: ",sum(as.logical(s30_g$Pass)),"/",nrow(s30_g),"\n",sep="")
cat("Table S29 inventory rows: ",nrow(t29),"\n",sep="")
cat("Table S30 gates: ",sum(as.logical(t30_g$Pass)),"/",nrow(t30_g),"\n",sep="")
cat("Artifacts: ",sum(file.exists(expected_files)),"/",length(expected_files),"\n",sep="")
cat("Blocking failures: ",blocking,"\n",sep="")
cat("Overall: ",ifelse(blocking==0,"PASS","FAIL"),"\n",sep="")
cat("Decision: ",ifelse(blocking==0,"READY_FOR_PUBLIC_FREEZE","BLOCKED"),"\n",sep="")
cat("============================================================\n")
if(blocking>0) stop("15 Additional robustness audit failed: ",paste(checks$Section[!checks$Pass],collapse=" | "),call.=FALSE)
