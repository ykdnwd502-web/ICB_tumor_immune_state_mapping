############################################################
## 07_CellChat_LR_make_S25.R
##
## Public/frozen reporting figure S25.
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

k20_file <- file.path(TABLE_DIR,"CosMx_LR_K20_all.csv")
if (!file.exists(k20_file)) stop("Missing canonical K20 input. Run 21J2A first.",call.=FALSE)

d <- utils::read.csv(k20_file,stringsAsFactors=FALSE,check.names=FALSE)
req <- c("pair_family","contrast","status","log2_spatial_enrichment",
         "neighbor_fraction_percent","FDR")
if (!all(req %in% names(d))) stop("K20 canonical source missing required columns.",call.=FALSE)
if (!("k" %in% names(d)) || !all(d$k==20)) stop("S25 source must be k=20 only.",call.=FALSE)
if (nrow(d)!=60L) stop("S25 locked test universe must contain 60 rows.",call.=FALSE)

d <- d %>%
  mutate(
    log2_spatial_enrichment=as.numeric(log2_spatial_enrichment),
    neighbor_fraction_percent=as.numeric(neighbor_fraction_percent),
    FDR=as.numeric(FDR),
    significant=ifelse(
      status=="ok" & is.finite(FDR) & FDR<0.05,
      "FDR < 0.05","FDR ≥ 0.05 / not evaluable"
    )
  )

pair_order <- d %>%
  group_by(pair_family) %>%
  summarise(mx=max(log2_spatial_enrichment,na.rm=TRUE),.groups="drop") %>%
  arrange(mx) %>%
  pull(pair_family)

contrast_order <- c(
  "Dediff/Stromal-high → Myeloid–Treg-high",
  "Myeloid–Treg-high → Dediff/Stromal-high",
  "Strict Dediff-only → Myeloid-only",
  "Strict Myeloid-only → Dediff-only",
  "Within Both-high niche"
)

d$pair_family <- factor(d$pair_family,levels=pair_order)
d$contrast <- factor(d$contrast,levels=contrast_order)

pA <- ggplot(d,aes(x=contrast,y=pair_family)) +
  geom_point(
    aes(size=neighbor_fraction_percent,fill=log2_spatial_enrichment,color=significant),
    shape=21,stroke=0.75,alpha=0.95
  ) +
  scale_fill_gradient2(
    low="#2C7BB6",mid="white",high="#D7191C",midpoint=0,
    name="log2 spatial\nenrichment"
  ) +
  scale_color_manual(
    values=c("FDR < 0.05"="black","FDR ≥ 0.05 / not evaluable"="grey75"),
    name=NULL
  ) +
  scale_size_continuous(range=c(1.5,7),name="Ligand-high /\nreceptor-high (%)") +
  labs(
    title="A  Primary CosMx spatial proximity analysis (k = 20)",
    subtitle="999 permutations; BH-FDR within the 60-test k=20 universe",
    x=NULL,y=NULL
  ) +
  theme_bw(base_size=12.5, base_family="sans") +
  theme(
    plot.title=element_text(face="bold",size=15),
    plot.subtitle=element_text(size=10.5),
    axis.text.x=element_text(angle=35,hjust=1,size=10.25),
    axis.text.y=element_text(size=10.25),
    legend.title=element_text(size=12,face="bold"),
    legend.text=element_text(size=10.5),
    panel.grid.minor=element_blank(),
    legend.position="right"
  )

rank_df <- d %>%
  mutate(
    pair_label=paste0(as.character(pair_family)," | ",as.character(contrast))
  ) %>%
  arrange(FDR,desc(log2_spatial_enrichment),pair_label) %>%
  mutate(rank=row_number())

pB <- ggplot(rank_df,aes(x=rank,y=log2_spatial_enrichment)) +
  geom_hline(yintercept=0,linetype=2,color="grey50") +
  geom_point(aes(fill=log2_spatial_enrichment,color=significant),shape=21,size=2.5,stroke=0.7) +
  scale_fill_gradient2(low="#2C7BB6",mid="white",high="#D7191C",midpoint=0,guide="none") +
  scale_color_manual(
    values=c("FDR < 0.05"="black","FDR ≥ 0.05 / not evaluable"="grey75"),
    name=NULL
  ) +
  labs(
    title="B  Ranked axis–contrast combinations",
    subtitle="Ranking is lexicographic: lower FDR, then higher log2 spatial enrichment",
    x="Rank among 60 k=20 tests",
    y="log2 spatial enrichment"
  ) +
  theme_bw(base_size=12.5, base_family="sans") +
  theme(
    plot.title=element_text(face="bold",size=15),
    plot.subtitle=element_text(size=10.5),
    axis.title=element_text(size=12.5,face="bold"),
    axis.text=element_text(size=10.75),
    legend.title=element_text(size=12,face="bold"),
    legend.text=element_text(size=10.5),
    panel.grid.minor=element_blank(),
    legend.position="bottom"
  )

fig <- pA / pB +
  plot_annotation(
    title="CosMx spatial proximity-constrained analysis of CellChat-nominated candidate axes",
    subtitle="Spatially grounded hypothesis support only; not receptor activation or functional ligand–receptor validation",
    theme=theme(
      plot.title=element_text(size=16,face="bold",family="sans"),
      plot.subtitle=element_text(size=11.5,family="sans")
    )
  )

prefix <- "Supplementary_Figure_S25_CosMx_spatial_proximity_candidate_LR_axes"
out_png <- file.path(FIG_DIR,paste0(prefix,".png"))
out_jpg <- file.path(FIG_DIR,paste0(prefix,".jpg"))
out_pdf <- file.path(FIG_DIR,paste0(prefix,".pdf"))

ggsave(out_png,fig,width=11.5,height=10.5,units="in",dpi=300,bg="white",limitsize=FALSE)
ggsave(out_jpg,fig,width=11.5,height=10.5,units="in",dpi=300,bg="white",limitsize=FALSE)
ggsave(out_pdf,fig,width=11.5,height=10.5,units="in",bg="white",limitsize=FALSE)

utils::write.csv(
  rank_df,
  file.path(FIG_DIR,"Supplementary_Figure_S25_ranked_axis_contrast_source.csv"),
  row.names=FALSE
)

audit <- data.frame(
  input_mode="fixed_public_K20",
  k20_file=normalizePath(k20_file,winslash="/",mustWork=TRUE),
  k=20,
  n_tests=nrow(d),
  n_contrasts=length(unique(d$contrast)),
  n_pairs=length(unique(d$pair_family)),
  n_perm=999,
  FDR_method="BH within k=20 test universe",
  fuzzy_autosearch=FALSE,
  interpretation="spatial hypothesis support only; not functional signaling",
  stringsAsFactors=FALSE
)
utils::write.csv(audit,file.path(FIG_DIR,paste0(prefix,"_audit.csv")),row.names=FALSE)

cat("\nCELLCHAT/LR S25 COMPLETED\n")
cat("K20 tests: ",nrow(d),"\n",sep="")
