############################################################
## 08_CellChat_LR_make_S26.R
##
## Public/frozen reporting figure S26.
## Fixed reporting-source inputs only; no fuzzy auto-search.
##
## STYLE ONLY:
##   typography harmonized to the frozen S20 reference.
##   Analytical inputs/results are unchanged.
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(PROJECT_DIR)) PROJECT_DIR <- "D:/ICB_resistance_project"
PROJECT_DIR <- normalizePath(PROJECT_DIR,winslash="/",mustWork=TRUE)

TABLE_DIR <- file.path(
  PROJECT_DIR,
  "results", "tables",
  "CellChat_LR_reporting"
)

FIG_DIR <- file.path(
  PROJECT_DIR,
  "results", "figures",
  "CellChat_LR_reporting"
)
dir.create(FIG_DIR,recursive=TRUE,showWarnings=FALSE)

ksum_file <- file.path(TABLE_DIR,"CosMx_LR_pair_summary_by_k.csv")
rob_file <- file.path(TABLE_DIR,"CosMx_LR_kNN_robustness.csv")

if (!file.exists(ksum_file)) stop("Missing canonical kNN summary.",call.=FALSE)
if (!file.exists(rob_file)) stop("Missing canonical robustness summary.",call.=FALSE)

ks <- utils::read.csv(ksum_file,stringsAsFactors=FALSE,check.names=FALSE)
rb <- utils::read.csv(rob_file,stringsAsFactors=FALSE,check.names=FALSE)

if (!setequal(unique(as.integer(ks$k)),c(10L,20L,30L))) stop("S26 requires k=10/20/30.",call.=FALSE)
if (nrow(ks)!=36L) stop("S26 kNN pair-summary source must contain 36 rows.",call.=FALSE)
if (nrow(rb)!=12L) stop("S26 robustness source must contain 12 axes.",call.=FALSE)

ks <- ks %>%
  mutate(
    k=factor(k,levels=c(10,20,30)),
    max_log2_spatial_enrichment=as.numeric(max_log2_spatial_enrichment),
    min_FDR=as.numeric(min_FDR),
    sig=ifelse(is.finite(min_FDR)&min_FDR<0.05,"*","")
  )

pair_order <- rb$pair_family[
  order(
    pmax(rb$k10_max_log2,rb$k20_max_log2,rb$k30_max_log2,na.rm=TRUE),
    decreasing=FALSE
  )
]
pair_order <- unique(as.character(pair_order))

ks$pair_family <- factor(ks$pair_family,levels=pair_order)

pA <- ggplot(ks,aes(x=k,y=pair_family,fill=max_log2_spatial_enrichment)) +
  geom_tile(color="white",linewidth=0.65) +
  geom_text(aes(label=sig),size=4.5,fontface="bold") +
  scale_fill_gradient2(
    low="#2C7BB6",mid="white",high="#D7191C",midpoint=0,
    name="Maximum log2\nspatial enrichment"
  ) +
  labs(
    title="A  kNN sensitivity of candidate ligand–receptor proximity",
    subtitle="Asterisks indicate minimum pair-level FDR < 0.05 within the corresponding k-specific test universe",
    x="k nearest neighbors",y=NULL
  ) +
  theme_bw(base_size=12.5, base_family="sans") +
  theme(
    plot.title=element_text(face="bold",size=15),
    plot.subtitle=element_text(size=10.5),
    axis.title=element_text(size=12.5,face="bold"),
    axis.text=element_text(size=10.75),
    axis.text.y=element_text(size=10.5),
    legend.title=element_text(size=12,face="bold"),
    legend.text=element_text(size=10.5),
    panel.grid=element_blank()
  )

rb <- rb %>%
  mutate(
    robustness_class=case_when(
      significant_positive_kNN_n==3 ~ "3/3 k settings",
      significant_positive_kNN_n==2 ~ "2/3 k settings",
      significant_positive_kNN_n==1 ~ "1/3 k settings",
      TRUE ~ "0/3 k settings"
    ),
    pair_family=factor(pair_family,levels=pair_order)
  )

pB <- ggplot(
  rb,
  aes(
    x=significant_positive_kNN_n,
    y=pair_family,
    size=mean_max_neighbor_fraction_percent,
    shape=robustness_class
  )
) +
  geom_point(alpha=0.9) +
  scale_x_continuous(breaks=0:3,limits=c(-0.15,3.15)) +
  scale_size_continuous(range=c(2.5,8),name="Mean maximum\nneighbor fraction (%)") +
  labs(
    title="B  Robustness classification across kNN settings",
    subtitle="x-axis counts k settings with positive maximum enrichment and FDR < 0.05",
    x="Number of significant-positive kNN settings",
    y=NULL,
    shape="Robustness"
  ) +
  theme_bw(base_size=12.5, base_family="sans") +
  theme(
    plot.title=element_text(face="bold",size=15),
    plot.subtitle=element_text(size=10.5),
    axis.title=element_text(size=12.5,face="bold"),
    axis.text=element_text(size=10.75),
    axis.text.y=element_text(size=10.5),
    legend.title=element_text(size=12,face="bold"),
    legend.text=element_text(size=10.5),
    panel.grid.minor=element_blank(),
    legend.position="right"
  )

fig <- pA / pB +
  plot_annotation(
    title="kNN sensitivity summary for CosMx spatial ligand–receptor proximity analysis",
    subtitle="Spatial-neighborhood sensitivity audit only; not functional ligand–receptor validation",
    theme=theme(
      plot.title=element_text(size=16,face="bold",family="sans"),
      plot.subtitle=element_text(size=11.5,family="sans")
    )
  )

prefix <- "Supplementary_Figure_S26_CosMx_LR_kNN_sensitivity"
out_png <- file.path(FIG_DIR,paste0(prefix,".png"))
out_jpg <- file.path(FIG_DIR,paste0(prefix,".jpg"))
out_pdf <- file.path(FIG_DIR,paste0(prefix,".pdf"))

ggsave(out_png,fig,width=10.8,height=9.5,units="in",dpi=300,bg="white",limitsize=FALSE)
ggsave(out_jpg,fig,width=10.8,height=9.5,units="in",dpi=300,bg="white",limitsize=FALSE)
ggsave(out_pdf,fig,width=10.8,height=9.5,units="in",bg="white",limitsize=FALSE)

audit <- data.frame(
  input_mode="fixed_public_k10_k20_k30",
  k_summary_file=normalizePath(ksum_file,winslash="/",mustWork=TRUE),
  robustness_file=normalizePath(rob_file,winslash="/",mustWork=TRUE),
  k_settings="10;20;30",
  pair_rows=nrow(ks),
  candidate_axes=nrow(rb),
  fuzzy_autosearch=FALSE,
  interpretation="spatial-neighborhood sensitivity audit only",
  stringsAsFactors=FALSE
)
utils::write.csv(audit,file.path(FIG_DIR,paste0(prefix,"_audit.csv")),row.names=FALSE)

cat("\nCELLCHAT/LR S26 COMPLETED\n")
cat("kNN summary rows: ",nrow(ks),"\n",sep="")
