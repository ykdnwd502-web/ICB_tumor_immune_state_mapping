############################################################
## 06_CellChat_LR_make_S24.R
##
## Public/frozen reporting figure S24.
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
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash="/", mustWork=TRUE)

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
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

nom_file <- file.path(TABLE_DIR, "CellChat_candidate_nomination.csv")
sum_file <- file.path(TABLE_DIR, "CellChat_candidate_summary.csv")
count_file <- file.path(TABLE_DIR, "CellChat_source_target_nonzero_counts.csv")

for (p in c(nom_file,sum_file,count_file)) {
  if (!file.exists(p)) stop("Missing canonical input: ", p, call.=FALSE)
}

nom <- utils::read.csv(nom_file, stringsAsFactors=FALSE, check.names=FALSE)
sm <- utils::read.csv(sum_file, stringsAsFactors=FALSE, check.names=FALSE)
cnt <- utils::read.csv(count_file, stringsAsFactors=FALSE, check.names=FALSE)

axis_order <- sm$pair_family[order(sm$max_prob, decreasing=TRUE, na.last=TRUE)]
axis_order <- unique(as.character(axis_order))

# Panel A: all tracked nonzero source-target candidate rows; no arbitrary top-N truncation.
pa_df <- nom %>%
  mutate(
    pair_family = factor(pair_family, levels=rev(axis_order)),
    source_target = factor(source_target, levels=rev(sort(unique(as.character(source_target))))),
    prob = as.numeric(prob),
    pval = as.numeric(pval)
  )

pA <- ggplot(pa_df, aes(x=pair_family, y=source_target)) +
  geom_point(
    aes(size=prob, fill=prob),
    shape=21, color="grey25", stroke=0.25, alpha=0.92
  ) +
  scale_fill_gradient(low="white", high="#2C7FB8", name="CellChat\nscore") +
  scale_size_continuous(range=c(1.5,6.0), name="CellChat\nscore") +
  labs(
    title="A  Candidate ligand–receptor scores across source–target groups",
    subtitle="All tracked nonzero CellChat nomination rows; no top-N display truncation",
    x=NULL, y=NULL
  ) +
  theme_bw(base_size=12.5, base_family="sans") +
  theme(
    plot.title=element_text(face="bold", size=15),
    plot.subtitle=element_text(size=10.5),
    axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=9.75),
    axis.text.y=element_text(size=9.25),
    legend.title=element_text(size=12, face="bold"),
    legend.text=element_text(size=10.5),
    panel.grid.minor=element_blank(),
    legend.position="right"
  )

pb_df <- sm %>%
  mutate(
    pair_family=factor(pair_family, levels=rev(axis_order)),
    max_prob=as.numeric(max_prob),
    nominated=ifelse(nominated_by_cellchat, "Nominated", "Not nominated")
  )

pB <- ggplot(pb_df, aes(x=max_prob, y=pair_family)) +
  geom_segment(aes(x=0, xend=max_prob, y=pair_family, yend=pair_family), linewidth=0.55, color="grey70") +
  geom_point(aes(shape=nominated), size=3.1) +
  scale_shape_manual(values=c("Nominated"=16, "Not nominated"=1), name=NULL) +
  labs(
    title="B  Maximum CellChat-inferred score by tracked axis",
    x="Maximum CellChat-inferred score", y=NULL
  ) +
  theme_bw(base_size=12.5, base_family="sans") +
  theme(
    plot.title=element_text(face="bold", size=15),
    axis.text=element_text(size=10.5),
    axis.text.y=element_text(size=10.25),
    legend.title=element_text(size=12, face="bold"),
    legend.text=element_text(size=10.5),
    panel.grid.minor=element_blank(),
    legend.position="bottom"
  )

# Panel C: count of all nonzero inferred CellChat interactions across source/target groups.
pc_df <- cnt %>%
  mutate(
    source=factor(source, levels=sort(unique(as.character(source)))),
    target=factor(target, levels=rev(sort(unique(as.character(target))))),
    nonzero_interaction_count=as.numeric(nonzero_interaction_count)
  )

pC <- ggplot(pc_df, aes(x=source, y=target, fill=nonzero_interaction_count)) +
  geom_tile(color="white", linewidth=0.35) +
  geom_text(aes(label=nonzero_interaction_count), size=3.4) +
  scale_fill_gradient(low="white", high="#7A0177", name="Nonzero\ninteractions") +
  labs(
    title="C  Aggregated nonzero inferred interaction counts",
    subtitle="Counts are across the full source–target CellChat output",
    x="Source group", y="Target group"
  ) +
  theme_bw(base_size=12.5, base_family="sans") +
  theme(
    plot.title=element_text(face="bold", size=15),
    plot.subtitle=element_text(size=10.5),
    axis.title=element_text(size=12.5, face="bold"),
    axis.text.x=element_text(angle=45, hjust=1, size=9.75),
    axis.text.y=element_text(size=10.25),
    legend.title=element_text(size=12, face="bold"),
    legend.text=element_text(size=10.5),
    panel.grid=element_blank()
  )

fig <- pA / (pB | pC) +
  plot_annotation(
    title="CellChat-inferred candidate ligand–receptor axis nomination",
    subtitle="GSE244983; candidate nomination only, not evidence of causal intercellular communication or functional ligand–receptor signaling",
    theme=theme(
      plot.title=element_text(size=16, face="bold", family="sans"),
      plot.subtitle=element_text(size=11.5, family="sans")
    )
  )

prefix <- "Supplementary_Figure_S24_CellChat_candidate_LR_axes"
out_png <- file.path(FIG_DIR, paste0(prefix, ".png"))
out_jpg <- file.path(FIG_DIR, paste0(prefix, ".jpg"))
out_pdf <- file.path(FIG_DIR, paste0(prefix, ".pdf"))

ggsave(out_png, fig, width=12.5, height=10.5, units="in", dpi=300, bg="white", limitsize=FALSE)
ggsave(out_jpg, fig, width=12.5, height=10.5, units="in", dpi=300, bg="white", limitsize=FALSE)
ggsave(out_pdf, fig, width=12.5, height=10.5, units="in", bg="white", limitsize=FALSE)

audit <- data.frame(
  input_mode="fixed_public_reporting_inputs",
  nomination_file=normalizePath(nom_file,winslash="/",mustWork=TRUE),
  summary_file=normalizePath(sum_file,winslash="/",mustWork=TRUE),
  count_file=normalizePath(count_file,winslash="/",mustWork=TRUE),
  nomination_rows=nrow(nom),
  candidate_axes=nrow(sm),
  source_target_count_rows=nrow(cnt),
  fuzzy_autosearch=FALSE,
  interpretation="candidate nomination only; not functional signaling",
  stringsAsFactors=FALSE
)
utils::write.csv(audit, file.path(FIG_DIR,paste0(prefix,"_audit.csv")), row.names=FALSE)

cat("\nCELLCHAT/LR S24 COMPLETED\n")
cat("Candidate axes: ", nrow(sm), "\n", sep="")
