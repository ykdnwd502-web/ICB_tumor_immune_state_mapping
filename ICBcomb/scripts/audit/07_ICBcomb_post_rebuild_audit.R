## ICBcomb post-rebuild blocking audit
options(stringsAsFactors=FALSE)
ROOT <- Sys.getenv("ICBCOMB_ROOT"); if(!nzchar(ROOT)) ROOT <- Sys.getenv("ICB_PUBLIC_ROOT"); if(!nzchar(ROOT)) ROOT <- getwd()
ROOT <- normalizePath(ROOT,winslash="/",mustWork=TRUE)

T <- file.path(ROOT,"results","tables","ICBcomb")
F <- file.path(ROOT,"results","figures")
R <- file.path(ROOT,"results","reports")

q <- read.csv(file.path(T,"ICBcomb_query_counts.csv"),check.names=FALSE)
p <- read.csv(file.path(T,"ICBcomb_positive_prioritization_rows.csv"),check.names=FALSE)
s <- read.csv(file.path(T,"ICBcomb_state_strategy_summary.csv"),check.names=FALSE)
c <- read.csv(file.path(T,"ICBcomb_manuscript_claim_support.csv"),check.names=FALSE)

g <- c(
  nrow(q)==3L && identical(as.integer(q$n_query_genes),c(98L,73L,81L)),
  nrow(p)==42L && all(as.numeric(p$NES)>0),
  nrow(s)==10L,
  nrow(c)==6L && all(as.logical(c$pass)),
  all(file.exists(file.path(F,c(
    "Figure9_ICBcomb_perturbation_prioritization.jpg",
    "Supplementary_Figure_S27_ICBcomb_query_strategy_mapping.jpg",
    "Supplementary_Figure_S28_ICBcomb_positive_evidence.jpg"
  )))),
  all(file.exists(file.path(R,c(
    "Supplementary_Table_S27A_ICBcomb_positive_rows.xlsx",
    "Supplementary_Table_S27B_ICBcomb_strategy_summary.xlsx"
  ))))
)
cat("ICBcomb post-rebuild gates:",sum(g),"/",length(g),"\n")
if(!all(g)) stop("ICBcomb post-rebuild audit failed.")
cat("Decision: READY_FOR_MODULE_FREEZE\n")
