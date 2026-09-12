# =============================================================================
# Supplementary_Figure_S10_build_final_v2.3.3_FINAL_layout_polish.R
#
# Based on:
#   Supplementary_Figure_S10_build_final_v2.3.2.1_FINAL_title_minus15.R
#
# Accepted modifications:
#   1. Panel D program order restored to biological program hierarchy
#   2. Panel C title unchanged:
#        "Average tumor-intrinsic program scores by malignant-cell subcluster"
#   3. Panel A gene labels enlarged
#   4. Panel B facet titles enlarged
#   5. Legend typography refined
#
# Frozen:
#   - canonical CSV sources
#   - program identity
#   - colors
#   - panel structure
#   - statistics
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if(!nzchar(project_dir)){
  project_dir <- "D:/ICB_resistance_project"
}

source_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "GSE244983",
  "malignant_subclustering"
)

figure_dir <- file.path(
  project_dir,
  "results",
  "figures",
  "all_supplementary_figures"
)

dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)

s10a <- read.csv(file.path(source_dir,
                           "GSE244983_SuppFigureS10A_marker_dotplot_source.csv"),
                 check.names=FALSE)

s10b <- read.csv(file.path(source_dir,
                           "GSE244983_SuppFigureS10B_program_z_distribution_source.csv"),
                 check.names=FALSE)

s10c <- read.csv(file.path(source_dir,
                           "GSE244983_SuppFigureS10C_mean_program_z_source.csv"),
                 check.names=FALSE)

s10d <- read.csv(file.path(source_dir,
                           "GSE244983_SuppFigureS10D_program_correlation_source.csv"),
                 check.names=FALSE)

program_display_levels <- c(
  "Tumor plasticity/dedifferentiation",
  "Stromal/ECM remodeling",
  "Melanocytic differentiation",
  "Antigen presentation",
  "IFN response"
)

s10b$ProgramDisplay <- factor(
  s10b$ProgramDisplay,
  levels=program_display_levels
)

s10c$ProgramDisplay <- factor(
  s10c$ProgramDisplay,
  levels=program_display_levels
)

s10d$ProgramRowDisplay <- factor(
  s10d$ProgramRowDisplay,
  levels=rev(program_display_levels)
)

s10d$ProgramColumnDisplay <- factor(
  s10d$ProgramColumnDisplay,
  levels=program_display_levels
)

theme_s10 <- theme_bw(base_size=13.8) +
  theme(
    plot.title=element_text(
      size=14.6625,
      face="bold",
      hjust=0.5,
      margin=margin(b=10)
    ),
    axis.title=element_text(size=13.8),
    axis.text=element_text(size=11.5),
    strip.text=element_text(size=12,face="bold"),
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    panel.grid=element_blank(),
    plot.margin=margin(12,14,12,14)
  )

pA <- ggplot(s10a,aes(
  x=Gene,
  y=MalignantSubcluster))+
  geom_point(aes(
    size=PercentExpressed,
    color=AverageExpressionScaledDisplay))+
  scale_color_gradient2(
    low="#4575B4",
    mid="white",
    high="#D73027",
    midpoint=0,
    name="Average\nExpression")+
  scale_size(
    range=c(1,7),
    name="Percent\nExpressed")+
  labs(
    title="Malignant-cell marker profile across malignant-cell subclusters",
    x=NULL,
    y="Malignant-cell subcluster")+
  theme_s10+
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1,
      size=10)
  )

pB <- ggplot(s10b,aes(
  x=MalignantSubcluster,
  y=ProgramZ,
  fill=MalignantSubcluster))+
  geom_violin(
    trim=FALSE,
    color="grey30")+
  geom_boxplot(
    width=0.08,
    fill="white",
    outlier.shape=NA)+
  facet_wrap(
    ~ProgramDisplay,
    ncol=3,
    scales="free_y")+
  labs(
    title="Tumor-intrinsic program scores across malignant-cell subclusters",
    x="Malignant-cell subcluster",
    y="Program z-score")+
  theme_s10+
  theme(
    legend.position="none",
    axis.text.x=element_text(
      angle=45,
      hjust=1,
      size=9.2)
  )

pC <- ggplot(s10c,aes(
  x=ProgramDisplay,
  y=MalignantSubcluster,
  fill=MeanProgramZ))+
  geom_tile(color="white")+
  geom_text(
    aes(label=sprintf("%.2f",MeanProgramZ)),
    size=3.45)+
  scale_fill_gradient2(
    low="#4575B4",
    mid="white",
    high="#D73027",
    midpoint=0,
    name="Average\nz-score")+
  labs(
    title="Average tumor-intrinsic program scores by malignant-cell subcluster",
    x=NULL,
    y="Malignant-cell subcluster")+
  theme_s10+
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1,
      size=9.2)
  )

pD <- ggplot(s10d,aes(
  x=ProgramColumnDisplay,
  y=ProgramRowDisplay,
  fill=SpearmanRho))+
  geom_tile(color="white")+
  geom_text(
    aes(label=sprintf("%.2f",SpearmanRho)),
    size=3.45)+
  scale_fill_gradient2(
    low="#4575B4",
    mid="white",
    high="#D73027",
    midpoint=0,
    limits=c(-1,1),
    name="Spearman\nrho")+
  labs(
    title="Correlation among tumor-intrinsic malignant-cell programs",
    x=NULL,
    y=NULL)+
  theme_s10+
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1,
      size=9.2)
  )

final <- (pA/pB)/(pC|pD)+
  plot_annotation(
    tag_levels="A",
    theme=theme(
      plot.tag=element_text(
        size=21,
        face="bold")
    )
  )

outfile <- file.path(
  figure_dir,
  "Supplementary Figure S10. Supplementary malignant-cell marker and program analyses in GSE244983"
)

ggsave(
  paste0(outfile,".jpg"),
  final,
  width=19,
  height=16.2,
  units="in",
  dpi=300,
  quality=100)

ggsave(
  paste0(outfile,".png"),
  final,
  width=19,
  height=16.2,
  units="in",
  dpi=300)

cat("Supplementary Figure S10 v2.3.3 FINAL layout polish completed\n")
