############################################################
## 03_make_Table_S29_statistical_analysis_summary.R
## Reporting-only Excel export of the frozen Table S29 inventory.
############################################################
options(stringsAsFactors=FALSE)
required <- c("readr","openxlsx")
miss <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)>0) stop("Missing package(s): ",paste(miss,collapse=", "),call.=FALSE)
PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR",unset="D:/ICB_resistance_project")
PROJECT_DIR <- normalizePath(PROJECT_DIR,winslash="/",mustWork=TRUE)
SRC_DIR <- file.path(PROJECT_DIR,"results","tables","additional_robustness","statistical_testing_inventory")
OUT_DIR <- file.path(PROJECT_DIR,"results","tables")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
inv <- readr::read_csv(file.path(SRC_DIR,"Table_S29_statistical_testing_inventory.csv"),show_col_types=FALSE) |> as.data.frame(check.names=FALSE)
readme <- readr::read_csv(file.path(SRC_DIR,"Table_S29_README.csv"),show_col_types=FALSE) |> as.data.frame(check.names=FALSE)
gate <- readr::read_csv(file.path(SRC_DIR,"Table_S29_inventory_gate.csv"),show_col_types=FALSE) |> as.data.frame(check.names=FALSE)
if(!all(as.logical(gate$Pass))) stop("Table S29 source inventory gate is not PASS.",call.=FALSE)

wb <- openxlsx::createWorkbook(creator="15_Additional_robustness")
openxlsx::addWorksheet(wb,"README")
openxlsx::writeData(wb,"README",readme)
openxlsx::addWorksheet(wb,"Table S29")
openxlsx::writeData(wb,"Table S29",inv)
header <- openxlsx::createStyle(textDecoration="bold",fgFill="#D9EAF7",border="Bottom",halign="center",valign="center",wrapText=TRUE)
body <- openxlsx::createStyle(valign="top",wrapText=TRUE)
openxlsx::addStyle(wb,"README",header,rows=1,cols=1:ncol(readme),gridExpand=TRUE)
openxlsx::addStyle(wb,"Table S29",header,rows=1,cols=1:ncol(inv),gridExpand=TRUE)
openxlsx::addStyle(wb,"Table S29",body,rows=2:(nrow(inv)+1),cols=1:ncol(inv),gridExpand=TRUE)
openxlsx::freezePane(wb,"Table S29",firstRow=TRUE,firstCol=TRUE)
openxlsx::setColWidths(wb,"Table S29",cols=1:ncol(inv),widths=c(7,25,22,20,38,35,30,28,34,35,34,34,28,32))
out <- file.path(OUT_DIR,"Supplementary Table S29. Statistical analysis and multiple-testing summary.xlsx")
openxlsx::saveWorkbook(wb,out,overwrite=TRUE)
cat("Saved Table S29: ",normalizePath(out,winslash="/",mustWork=FALSE),"\n",sep="")
