############################################################
## 09_CellChat_LR_make_Supplementary_Tables_S24_S26.R
##
## Public/frozen Supplementary Tables S24-S26 producer.
## Reads only Step05 clean reporting source tables.
############################################################

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(PROJECT_DIR)) PROJECT_DIR <- "D:/ICB_resistance_project"
PROJECT_DIR <- normalizePath(PROJECT_DIR,winslash="/",mustWork=TRUE)

if (!requireNamespace("openxlsx",quietly=TRUE)) {
  stop("Package 'openxlsx' is required.",call.=FALSE)
}

TABLE_DIR <- file.path(
  PROJECT_DIR,
  "results", "tables",
  "CellChat_LR_reporting"
)

REPORT_DIR <- TABLE_DIR

dir.create(
  REPORT_DIR,
  recursive=TRUE,
  showWarnings=FALSE
)

readc <- function(nm) {
  p <- file.path(TABLE_DIR,nm)
  if (!file.exists(p)) stop("Missing canonical table: ",p,call.=FALSE)
  utils::read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
}

nom <- readc("CellChat_candidate_nomination.csv")
sum16 <- readc("CellChat_candidate_summary.csv")
pairs <- readc("CellChat_candidate_pairs.csv")
counts <- readc("CellChat_source_target_nonzero_counts.csv")
k20all <- readc("CosMx_LR_K20_all.csv")
k20sum <- readc("CosMx_LR_K20_summary.csv")
alls <- readc("CosMx_LR_all_kNN.csv")
ksum <- readc("CosMx_LR_pair_summary_by_k.csv")
rob <- readc("CosMx_LR_kNN_robustness.csv")

header_style <- openxlsx::createStyle(
  fontColour="#FFFFFF",fgFill="#4F81BD",
  textDecoration="bold",halign="center",valign="center",
  border="Bottom"
)
wrap_style <- openxlsx::createStyle(wrapText=TRUE,valign="top")

add_df <- function(wb,sheet,df) {
  openxlsx::addWorksheet(wb,sheet)
  openxlsx::writeData(wb,sheet,df,withFilter=TRUE,headerStyle=header_style)
  openxlsx::freezePane(wb,sheet,firstRow=TRUE)
  openxlsx::setColWidths(wb,sheet,cols=seq_len(ncol(df)),widths="auto")
  openxlsx::addStyle(
    wb,sheet,wrap_style,
    rows=1:(nrow(df)+1),cols=seq_len(ncol(df)),
    gridExpand=TRUE,stack=TRUE
  )
}

make_readme <- function(table_id,purpose) {
  data.frame(
    Item=c("Table","Purpose","Interpretation boundary","Reporting source"),
    Value=c(
      table_id,
      purpose,
      "Candidate nomination/spatial hypothesis support only; not receptor activation, functional signaling, causal mechanism, therapeutic efficacy, or response validation.",
      "results/canonical/CellChat_LR_v1.0/"
    ),
    stringsAsFactors=FALSE
  )
}

# S24
wb24 <- openxlsx::createWorkbook()
add_df(wb24,"README",make_readme(
  "Supplementary Table S24",
  "CellChat candidate ligand–receptor nomination in GSE244983."
))
add_df(wb24,"Candidate_nomination",nom)
add_df(wb24,"Candidate_summary",sum16)
add_df(wb24,"Candidate_pairs",pairs)
add_df(wb24,"Source_target_counts",counts)

out24 <- file.path(REPORT_DIR,"Supplementary_Table_S24_CellChat_candidate_LR.xlsx")
openxlsx::saveWorkbook(wb24,out24,overwrite=TRUE)

# S25
wb25 <- openxlsx::createWorkbook()
add_df(wb25,"README",make_readme(
  "Supplementary Table S25",
  "Primary CosMx spatial proximity-constrained analysis of CellChat-nominated axes at k=20 with 999 permutations."
))
add_df(wb25,"K20_all_tests",k20all)
add_df(wb25,"K20_pair_summary",k20sum)

method25 <- data.frame(
  parameter=c("k","n_perm","high definition","low-positive rule","FDR","test universe"),
  value=c("20","999","top 25% among expressing cells","<50 positive cells => all positive cells are high","Benjamini-Hochberg","12 axes x 5 contrasts = 60 tests"),
  stringsAsFactors=FALSE
)
add_df(wb25,"Method_contract",method25)

out25 <- file.path(REPORT_DIR,"Supplementary_Table_S25_CosMx_LR_K20.xlsx")
openxlsx::saveWorkbook(wb25,out25,overwrite=TRUE)

# S26
wb26 <- openxlsx::createWorkbook()
add_df(wb26,"README",make_readme(
  "Supplementary Table S26",
  "kNN sensitivity audit for CosMx spatial ligand–receptor proximity at k=10,20,30."
))
add_df(wb26,"All_kNN_tests",alls)
add_df(wb26,"Pair_summary_by_k",ksum)
add_df(wb26,"Robustness_summary",rob)

method26 <- data.frame(
  parameter=c("k settings","primary k","n_perm per k","high definition","FDR","interpretation"),
  value=c("10;20;30","20","999","top 25% among expressing cells; <50 positives => all positive","BH separately within each k-specific universe","spatial-neighborhood sensitivity audit only"),
  stringsAsFactors=FALSE
)
add_df(wb26,"Method_contract",method26)

out26 <- file.path(REPORT_DIR,"Supplementary_Table_S26_CosMx_LR_kNN_sensitivity.xlsx")
openxlsx::saveWorkbook(wb26,out26,overwrite=TRUE)

cat("\nCELLCHAT/LR TABLES S24-S26 COMPLETED\n")
cat(out24,"\n",out25,"\n",out26,"\n",sep="")
