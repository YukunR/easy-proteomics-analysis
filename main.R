# ===========================================
# ==== Proteomics Data Analysis Pipeline ====
# ===========================================
# This script provides a complete proteomics analysis workflow with 
# automatic checkpoint and resume functionality

# ===========================
# ========= Configs =========
# ===========================

# Usually, you only need to modify this section 
# unless you know exactly what you're doing

# ==== Normalization ====
base_dir <- "./res/"
protein_expr_file <- "./data/origin_data.txt"
sample_info_file <- "./data/sample_info.txt"

na_threshold <- 0.6  # Proteins with an NA ratio greater than this value will be discarded
normalization_method <- "global"  # should be one of "global" and "within_group"
imputation_method <- "knn"  # should be one of "knn" and "perseus"

# If you have preferred colors uncomment and run the line below.
# Otherwise, just leave the line below commented out
# custom_colors <- c("#FF5733", "#33C3FF", "#75FF33")

# ==== Differential expression analysis ====
# Ensure that the group names below match the group names used in sample_info
comparisons <- list(
  list(control = "HC", treatment = "NC", name = "HC_vs_NC"),
  list(control = "HD", treatment = "HC", name = "HD_vs_HC")
)

# ==== Enrichment analysis ====
go_background_file <- "./data/all_uniprot_go_background.csv"
kegg_background_file <- "./data/pathfromKegg_mmu.txt"

# ==================================
# ========= Usage Examples =========
# ==================================

# # Run complete analysis pipeline
# result <- run_proteomics_analysis()
# 
# # Check project status
# check_project_status()
# 
# # If interrupted, run the same command to resume
# result <- run_proteomics_analysis()
# 
# # Force rerun specific steps
# rerun_step("normalization")  # Rerun normalization step
# rerun_step(c("pca", "differential_analysis_KI_vs_Old"))  # Rerun multiple steps
# 
# # Clean project and restart
# clean_project()

# ================================
# ========= Run Pipeline =========
# ================================
source("analysis_steps.R")
result <- run_proteomics_analysis()
