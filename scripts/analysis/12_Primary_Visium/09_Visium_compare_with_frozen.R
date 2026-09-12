############################################################
## 09_Visium_compare_with_frozen.R
## FINAL CHECK ONLY — no analysis.
## Scientific/core CSV outputs are compared with FROZEN_BACKUP.
## Figure files are checked for presence because reporting style was
## intentionally standardized and byte identity is not expected.
############################################################

options(stringsAsFactors = FALSE)
PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR", unset = "D:/ICB_resistance_project")
FROZEN_ROOT <- Sys.getenv("ICB_FROZEN_ROOT", unset = "D:/ICB_resistance_project_FROZEN_BACKUP")
for (x in c(PROJECT_DIR,FROZEN_ROOT)) if (!dir.exists(x)) stop("Required root missing: ",x,call.=FALSE)
PROJECT_DIR <- normalizePath(PROJECT_DIR,winslash="/",mustWork=TRUE)
FROZEN_ROOT <- normalizePath(FROZEN_ROOT,winslash="/",mustWork=TRUE)
AUDIT_DIR <- file.path(PROJECT_DIR,"results","audit","12_Primary_Visium_sequential_reproduction_check")
dir.create(AUDIT_DIR,recursive=TRUE,showWarnings=FALSE)

compare_csv <- function(cur, fro, group, tol=1e-10) {
  if (!file.exists(cur) || !file.exists(fro)) return(data.frame(group=group,file=basename(cur),current_exists=file.exists(cur),frozen_exists=file.exists(fro),rows_exact=FALSE,columns_exact=FALSE,max_numeric_abs_diff=Inf,nonnumeric_exact=FALSE,pass=FALSE,stringsAsFactors=FALSE))
  a <- utils::read.csv(cur,check.names=FALSE,stringsAsFactors=FALSE)
  b <- utils::read.csv(fro,check.names=FALSE,stringsAsFactors=FALSE)
  if (nrow(a)!=nrow(b) || !identical(names(a),names(b))) return(data.frame(group=group,file=basename(cur),current_exists=TRUE,frozen_exists=TRUE,rows_exact=nrow(a)==nrow(b),columns_exact=identical(names(a),names(b)),max_numeric_abs_diff=Inf,nonnumeric_exact=FALSE,pass=FALSE,stringsAsFactors=FALSE))
  md <- 0; ne <- TRUE
  for (nm in names(a)) {
    aa<-a[[nm]]; bb<-b[[nm]]
    if (is.numeric(aa) && is.numeric(bb)) {
      if (!identical(is.na(aa),is.na(bb))) {md<-Inf; next}
      ok<-!is.na(aa)
      if (any(ok)) {
        dd<-abs(aa[ok]-bb[ok]); fi<-is.finite(dd)
        if(any(fi)) md<-max(md,max(dd[fi]))
        if(any(!fi) && !identical(aa[ok][!fi],bb[ok][!fi])) md<-Inf
      }
    } else {
      aa<-as.character(aa); bb<-as.character(bb)
      if(!identical(is.na(aa),is.na(bb))) ne<-FALSE else {ok<-!is.na(aa); if(!identical(aa[ok],bb[ok])) ne<-FALSE}
    }
  }
  data.frame(group=group,file=basename(cur),current_exists=TRUE,frozen_exists=TRUE,rows_exact=TRUE,columns_exact=TRUE,max_numeric_abs_diff=md,nonnumeric_exact=ne,pass=is.finite(md)&&md<=tol&&ne,stringsAsFactors=FALSE)
}

file_pairs <- list(
  Step01_state_scoring=data.frame(
    current=c(
      "results/tables/spatial_melanoma_validation/Step11A_spatial_final_four_ICB_state_scores.csv",
      "results/tables/spatial_melanoma_validation/Step11A_spatial_dominant_state_counts.csv",
      "results/tables/spatial_melanoma_validation/Step11A_spatial_final_state_signature_mapping.csv",
      "results/tables/spatial_melanoma_validation/Step11A_spatial_gene_set_presence_summary.csv"
    ),
    frozen_legacy=c(
      "results/tables/spatial_melanoma_validation/Step11A_spatial_final_four_ICB_state_scores_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11A_spatial_dominant_state_counts_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11A_spatial_final_state_signature_mapping_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11A_spatial_gene_set_presence_summary_CLEAN.csv"
    ),
    stringsAsFactors=FALSE
  ),
  Step02_spatial_QC_Moran=data.frame(
    current=c(
      "results/tables/spatial_melanoma_validation/Step11_spatial_QC_state_Spearman_correlation.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_raw_vs_QCresidual_state_correlation.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_Morans_I_raw_and_QCresidual.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_high_state_neighbor_enrichment_top25_pooled.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_dual_high_spot_level_OR_top25_pooled.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_dual_high_threshold_sensitivity_OR.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_coordinate_source_and_merge_audit.csv"
    ),
    frozen_legacy=c(
      "results/tables/spatial_melanoma_validation/Step11_spatial_QC_state_Spearman_correlation_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_raw_vs_QCresidual_state_correlation_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_Morans_I_raw_and_QCresidual_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_high_state_neighbor_enrichment_top25_pooled_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_dual_high_spot_level_OR_top25_pooled_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_dual_high_threshold_sensitivity_OR_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step11_spatial_coordinate_source_and_merge_audit_CLEAN.csv"
    ),
    stringsAsFactors=FALSE
  ),
  Step03_dual_high_niche=data.frame(
    current=c(
      "results/tables/spatial_melanoma_validation/Step12_dual_high_cutoff_audit.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_group_counts.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_spot_level_OR_top25_pooled.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_candidate_gene_screen_Both_vs_Neither.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_top_candidate_genes_Figure8C.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_Figure8D_four_group_candidate_gene_expression_summary.csv"
    ),
    frozen_legacy=c(
      "results/tables/spatial_melanoma_validation/Step12_dual_high_cutoff_audit_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_group_counts_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_spot_level_OR_top25_pooled_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_candidate_gene_screen_Both_vs_Neither_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_top_candidate_genes_Figure8C_CLEAN.csv",
      "results/tables/spatial_melanoma_validation/Step12_dual_high_Figure8D_four_group_candidate_gene_expression_summary_CLEAN.csv"
    ),
    stringsAsFactors=FALSE
  ),
  Step04_threshold_independent=data.frame(
    current=c(
      "results/tables/revision_spatial_covariation/Visium_continuous_state_correlation_and_partial_correlation.csv",
      "results/tables/revision_spatial_covariation/Visium_bivariate_spatial_Moran_kNN_sensitivity.csv",
      "results/tables/revision_spatial_covariation/Visium_threshold_grid_dual_high_Fisher_OR_sensitivity.csv",
      "results/tables/revision_spatial_covariation/Visium_local_bivariate_colocalization_scores_k6_primary.csv"
    ),
    frozen_legacy=c(
      "results/tables/revision_spatial_covariation/Visium_continuous_state_correlation_and_partial_correlation_FIXED.csv",
      "results/tables/revision_spatial_covariation/Visium_bivariate_spatial_Moran_kNN_sensitivity_FIXED.csv",
      "results/tables/revision_spatial_covariation/Visium_threshold_grid_dual_high_Fisher_OR_sensitivity_FIXED.csv",
      "results/tables/revision_spatial_covariation/Visium_local_bivariate_colocalization_scores_k6_primary_FIXED.csv"
    ),
    stringsAsFactors=FALSE
  ),
  Step05_composition=data.frame(
    current=c(
      "results/tables/revision_visium_composition/18E_Visium_spot_state_QC_composition_table.csv",
      "results/tables/revision_visium_composition/18E_composition_BothHigh_vs_NeitherHigh_Wilcoxon.csv",
      "results/tables/revision_visium_composition/18E_Visium_state_composition_Spearman_correlation.csv",
      "results/tables/revision_visium_composition/18E_Visium_spot_state_composition_residualized_scores.csv",
      "results/tables/revision_visium_composition/18E_raw_vs_residualized_dual_high_Fisher_OR.csv",
      "results/tables/revision_visium_composition/18E_raw_vs_residualized_state_score_Spearman_correlation.csv"
    ),
    frozen_legacy=c(
      "results/tables/revision_visium_composition/18E_Visium_spot_state_QC_composition_table.csv",
      "results/tables/revision_visium_composition/18E_composition_BothHigh_vs_NeitherHigh_Wilcoxon.csv",
      "results/tables/revision_visium_composition/18E_Visium_state_composition_Spearman_correlation.csv",
      "results/tables/revision_visium_composition/18E_Visium_spot_state_composition_residualized_scores.csv",
      "results/tables/revision_visium_composition/18E_raw_vs_residualized_dual_high_Fisher_OR.csv",
      "results/tables/revision_visium_composition/18E_raw_vs_residualized_state_score_Spearman_correlation.csv"
    ),
    stringsAsFactors=FALSE
  ),
  Step06_TLS_context=data.frame(
    current=c(
      "results/tables/revision_humoral_TLS_like/18F_humoral_TLS_like_module_scores.csv",
      "results/tables/revision_humoral_TLS_like/18F_Visium_spot_state_composition_TLS_like_table.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_module_summary_by_raw_dual_high_category.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_module_correlations_with_states_and_composition.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_module_univariate_spatial_Moran.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_bivariate_spatial_colocalization_Moran.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_operational_cohigh_enrichment_BothHigh_vs_NeitherHigh.csv"
    ),
    frozen_legacy=c(
      "results/tables/revision_humoral_TLS_like/18F_humoral_TLS_like_module_scores.csv",
      "results/tables/revision_humoral_TLS_like/18F_Visium_spot_state_composition_TLS_like_table.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_module_summary_by_raw_dual_high_category.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_module_correlations_with_states_and_composition.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_module_univariate_spatial_Moran.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_bivariate_spatial_colocalization_Moran.csv",
      "results/tables/revision_humoral_TLS_like/18F_TLS_like_operational_cohigh_enrichment_BothHigh_vs_NeitherHigh.csv"
    ),
    stringsAsFactors=FALSE
  ),
  Step07_primary_canonical=data.frame(
    current=c(
      "results/tables/primary_visium_spatial/panel_sources/S15/21_S15A_dominant_state_spatial_map.csv",
      "results/tables/primary_visium_spatial/panel_sources/S15/21_S15B_corrected_univariate_Moran_raw_QC_k6.csv",
      "results/tables/primary_visium_spatial/panel_sources/S15/21_S15C_corrected_high_state_neighbor_OR_DESCRIPTIVE.csv",
      "results/tables/primary_visium_spatial/panel_sources/S16/21_S16A_continuous_primary_pair_summary.csv",
      "results/tables/primary_visium_spatial/panel_sources/S16/21_S16B_symmetric_bivariate_Moran_k4_k6_k8_k12.csv",
      "results/tables/primary_visium_spatial/panel_sources/S16/21_S16C_threshold_grid_dual_high_OR_DESCRIPTIVE.csv",
      "results/tables/primary_visium_spatial/panel_sources/S17/21_S17B_state_composition_Spearman.csv",
      "results/tables/primary_visium_spatial/panel_sources/S18/21_S18D_bivariate_Moran_and_permutation_k6.csv",
      "results/tables/primary_visium_spatial/panel_sources/S19/21_S19C_humoral_TLS_module_spatial_autocorrelation_k_sensitivity.csv",
      "results/tables/primary_visium_spatial/panel_sources/S19/21_S19D_operational_humoral_TLS_cohigh_features.csv"
    ),
    frozen_legacy=c(
      "results/tables/primary_visium_spatial/panel_sources/S15/21_S15A_dominant_state_spatial_map.csv",
      "results/tables/primary_visium_spatial/panel_sources/S15/21_S15B_corrected_univariate_Moran_raw_QC_k6.csv",
      "results/tables/primary_visium_spatial/panel_sources/S15/21_S15C_corrected_high_state_neighbor_OR_DESCRIPTIVE.csv",
      "results/tables/primary_visium_spatial/panel_sources/S16/21_S16A_continuous_primary_pair_summary.csv",
      "results/tables/primary_visium_spatial/panel_sources/S16/21_S16B_symmetric_bivariate_Moran_k4_k6_k8_k12.csv",
      "results/tables/primary_visium_spatial/panel_sources/S16/21_S16C_threshold_grid_dual_high_OR_DESCRIPTIVE.csv",
      "results/tables/primary_visium_spatial/panel_sources/S17/21_S17B_state_composition_Spearman.csv",
      "results/tables/primary_visium_spatial/panel_sources/S18/21_S18D_bivariate_Moran_and_permutation_k6.csv",
      "results/tables/primary_visium_spatial/panel_sources/S19/21_S19C_humoral_TLS_module_spatial_autocorrelation_k_sensitivity.csv",
      "results/tables/primary_visium_spatial/panel_sources/S19/21_S19D_operational_humoral_TLS_cohigh_features.csv"
    ),
    stringsAsFactors=FALSE
  )
)
rows<-list(); i<-1L
for(g in names(file_pairs)) {
  z <- file_pairs[[g]]
  for(k in seq_len(nrow(z))) {
    rows[[i]] <- compare_csv(file.path(PROJECT_DIR,z$current[k]), file.path(FROZEN_ROOT,z$frozen_legacy[k]), g)
    i <- i + 1L
  }
}
cmp<-do.call(rbind,rows)
utils::write.csv(cmp,file.path(AUDIT_DIR,"01_CORE_TABLES_vs_FROZEN.csv"),row.names=FALSE)
summary<-do.call(rbind,lapply(split(cmp,cmp$group),function(z)data.frame(group=z$group[1],passed=sum(z$pass),total=nrow(z),pass=all(z$pass),stringsAsFactors=FALSE)))
rownames(summary)<-NULL
utils::write.csv(summary,file.path(AUDIT_DIR,"02_GROUP_SUMMARY.csv"),row.names=FALSE)

# Scientific contract.
state_scores<-utils::read.csv(file.path(PROJECT_DIR,"results/tables/spatial_melanoma_validation/Step11A_spatial_final_four_ICB_state_scores.csv"),check.names=FALSE)
s15a<-utils::read.csv(file.path(PROJECT_DIR,"results/tables/primary_visium_spatial/panel_sources/S15/21_S15A_dominant_state_spatial_map.csv"),check.names=FALSE)
s19b<-utils::read.csv(file.path(PROJECT_DIR,"results/tables/primary_visium_spatial/panel_sources/S19/21_S19B_humoral_TLS_module_spatial_scores_long.csv"),check.names=FALSE)
contract<-data.frame(check=c("state_score_spots","S15_spots","S19_unique_spots","four_state_columns_present"),observed=c(nrow(state_scores),nrow(s15a),length(unique(s19b$Spot)),all(c("Immune_defective_Cold","Myeloid_Treg_Immunosuppressive","Tumor_dedifferentiation_Stromal_remodeling","Melanocytic_Differentiation") %in% names(state_scores))),expected=c("3458","3458","3458","TRUE"),pass=c(nrow(state_scores)==3458L,nrow(s15a)==3458L,length(unique(s19b$Spot))==3458L,all(c("Immune_defective_Cold","Myeloid_Treg_Immunosuppressive","Tumor_dedifferentiation_Stromal_remodeling","Melanocytic_Differentiation") %in% names(state_scores))),stringsAsFactors=FALSE)
utils::write.csv(contract,file.path(AUDIT_DIR,"03_REPRODUCTION_CONTRACT.csv"),row.names=FALSE)

fig_current<-c(
 "results/figures/spatial_melanoma_validation/Figure8A_dual_high_spatial_colocalization.png",
 "results/figures/spatial_melanoma_validation/Figure8B_dual_high_enrichment_OR.png",
 "results/figures/spatial_melanoma_validation/Figure8C_candidate_mechanism_genes_Both_vs_Neither.png",
 "results/figures/spatial_melanoma_validation/Figure8D_four_group_candidate_gene_dotplot.png",
 paste0("results/figures/primary_visium_spatial/Supplementary_Figure_S",15:19,".pdf"),
 paste0("results/figures/primary_visium_spatial/Supplementary_Figure_S",15:19,".jpg")
)
fig_frozen_legacy<-c(
 "results/figures/spatial_melanoma_validation/Figure8A_dual_high_spatial_colocalization_STEP11_CLEAN.png",
 "results/figures/spatial_melanoma_validation/Figure8B_dual_high_enrichment_OR_CLEAN.png",
 "results/figures/spatial_melanoma_validation/Figure8C_candidate_mechanism_genes_Both_vs_Neither_CLEAN.png",
 "results/figures/spatial_melanoma_validation/Figure8D_four_group_candidate_gene_dotplot_CLEAN.png",
 paste0("results/figures/primary_visium_spatial/Supplementary_Figure_S",15:19,"_CANONICAL.pdf"),
 paste0("results/figures/primary_visium_spatial/Supplementary_Figure_S",15:19,"_CANONICAL.jpg")
)
figs<-data.frame(
  current_relative_path=fig_current,
  frozen_legacy_relative_path=fig_frozen_legacy,
  current_exists=file.exists(file.path(PROJECT_DIR,fig_current)),
  frozen_exists=file.exists(file.path(FROZEN_ROOT,fig_frozen_legacy)),
  stringsAsFactors=FALSE
)
# Frozen presence is informative; style-standardized current files need only current presence.
figs$pass<-figs$current_exists
utils::write.csv(figs,file.path(AUDIT_DIR,"04_FIGURE_PRESENCE.csv"),row.names=FALSE)

blocking<-sum(!cmp$pass)+sum(!contract$pass)+sum(!figs$pass)
cat("\\n============================================================\\n")
cat("12_Primary_Visium SEQUENTIAL REPRODUCTION CHECK\\n")
cat("============================================================\\n")
for(j in seq_len(nrow(summary))) cat(summary$group[j],": ",summary$passed[j],"/",summary$total[j],ifelse(summary$pass[j]," PASS"," FAIL"),"\\n",sep="")
cat("Reproduction contract: ",sum(contract$pass),"/",nrow(contract)," PASS\\n",sep="")
cat("Public figure files: ",sum(figs$pass),"/",nrow(figs)," PRESENT\\n",sep="")
cat("Blocking failures: ",blocking,"\\n",sep="")
cat("Overall: ",ifelse(blocking==0L,"PASS","FAIL"),"\\n",sep="")
cat("Audit: ",AUDIT_DIR,"\\n",sep="")
cat("============================================================\\n")
if(blocking!=0L){if(any(!cmp$pass))print(cmp[!cmp$pass,,drop=FALSE],row.names=FALSE);if(any(!contract$pass))print(contract[!contract$pass,,drop=FALSE],row.names=FALSE);stop("Primary Visium reproduction comparison failed; inspect audit CSVs.",call.=FALSE)}
