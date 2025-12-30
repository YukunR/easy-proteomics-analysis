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
imputation_method <- "auto"  # should be one of "auto", "knn" or "perseus"

# If you have preferred colors uncomment and run the line below.
# Otherwise, just leave the line below commented out
# custom_colors <- c("#FF5733", "#33C3FF", "#75FF33")

# ==== Threshold Settings ====
# FC threshold mode:
#   - "auto": automatically calculate using coverage analysis
#     (global_fc_threshold and comparison fc_threshold will be IGNORED with warning)
#   - "global": use global_fc_threshold for all comparisons
#   - "per_comparison": use fc_threshold in each comparison, or ask interactively if not specified
fc_threshold_mode <- "auto"
global_fc_threshold <- 2  # Used when fc_threshold_mode = "global"

# P-value threshold mode (default: 0.05):
#   - "global": use global_p_threshold for all comparisons
#   - "per_comparison": use p_threshold in each comparison, or ask interactively if not specified
p_threshold_mode <- "global"
global_p_threshold <- 0.05  # Used when p_threshold_mode = "global"

# ==== Differential expression analysis ====
# Ensure that the group names below match the group names used in sample_info
comparisons <- list(
  list(control = "HC", treatment = "NC", name = "HC_vs_NC"),
  list(control = "HD", treatment = "HC", name = "HD_vs_HC", fc_threshold = 2),
  list(control = c("HC", "HD"), treatment = "NC", name = "NC_vs_Rest"),
  list(control = "HD", treatment = c("HC", "NC"), name = "Rest_vs_HD"),
  list(control = "HC", treatment = c("HD", "NC"), name = "Rest_vs_HC")
  # This causes an error because NC appears in both the control and treatment groups.
  # list(control = c"HC", treatment = c("HD", "NC"), name = "Rest_vs_HC")
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
