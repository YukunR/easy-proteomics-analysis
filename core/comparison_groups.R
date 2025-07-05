# ==== Example usage ====
# Create all comparisons
# all_comparisons <- create_all_pairwise_comparisons(
#   expression_data = processed_data,
#   sample_info = sample_info,
#   control_groups = c("Control"),
#   exclude_groups = NULL
# )

# Control vs all other groups
# control_comparisons <- create_control_vs_all_comparisons(
#   expression_data = processed_data,
#   sample_info = sample_info,
#   control_group = "Control"
# )

# Control vs specific groups
# control_comparisons <- create_control_vs_all_comparisons(
#   expression_data = processed_data,
#   sample_info = sample_info,
#   control_group = "Control",
#   treatment_groups = c("Treatment_A", "High_dose")  # 只比较这些组
# )

library(dplyr)
library(tidyr)

#' Generate comparison groups for differential analysis
#' 
#' @description
#' Creates comparison groups by automatically selecting samples based on user-defined
#' control and treatment group combinations
#' 
#' @param expression_data Expression data matrix (with or without ID column)
#' @param sample_info Sample information data frame
#' @param comparisons List of comparisons to perform. Each element should be a named list with:
#'   \itemize{
#'     \item control: Name of control group
#'     \item treatment: Name of treatment group  
#'     \item name: Optional custom name for comparison (default: "treatment_vs_control")
#'   }
#' @param id_col ID column name in expression data, default "Accession"
#' @param sample_col Sample column name in sample_info, default "Sample"
#' @param group_col Group column name in sample_info, default "Group"
#' @param include_id_col Whether to include ID column in output, default TRUE
#' 
#' @return Named list of comparison groups, each containing:
#'   \item{expression_data}{Filtered expression data for the comparison}
#'   \item{sample_info}{Filtered sample info for the comparison}
#'   \item{control_name}{Name of control group}
#'   \item{treatment_name}{Name of treatment group}
#'   \item{n_control}{Number of control samples}
#'   \item{n_treatment}{Number of treatment samples}
#' 
#' @export
#' 
#' @examples
#' # Define comparisons
#' my_comparisons <- list(
#'   list(control = "Control", treatment = "Treatment_A"),
#'   list(control = "Control", treatment = "Treatment_B", name = "TreatB_vs_Ctrl"),
#'   list(control = "Treatment_A", treatment = "Treatment_B")
#' )
#' 
#' # Generate comparison groups
#' comparison_groups <- create_comparison_groups(
#'   expression_data = processed_data,
#'   sample_info = sample_info,
#'   comparisons = my_comparisons
#' )
create_comparison_groups <- function(expression_data,
                                     sample_info,
                                     comparisons,
                                     id_col = "Accession",
                                     sample_col = "Sample",
                                     group_col = "Group",
                                     include_id_col = TRUE) {
  
  # Input validation
  if (!is.list(comparisons)) {
    stop("'comparisons' must be a list")
  }
  
  if (!is.data.frame(sample_info)) {
    stop("'sample_info' must be a data frame")
  }
  
  # Check required columns
  required_cols <- c(sample_col, group_col)
  missing_cols <- setdiff(required_cols, colnames(sample_info))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in sample_info:", paste(missing_cols, collapse = ", ")))
  }
  
  # Prepare expression data
  if (id_col %in% colnames(expression_data)) {
    expr_matrix <- expression_data
    has_id_col <- TRUE
  } else {
    expr_matrix <- expression_data
    has_id_col <- FALSE
  }
  
  # Get available groups
  available_groups <- unique(sample_info[[group_col]])
  cat("Available groups:", paste(available_groups, collapse = ", "), "\n\n")
  
  # Initialize results list
  comparison_results <- list()
  
  # Process each comparison
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    
    # Validate comparison structure
    if (!("control" %in% names(comp)) || !("treatment" %in% names(comp))) {
      stop(paste("Comparison", i, "must have 'control' and 'treatment' fields"))
    }
    
    control_group <- comp$control
    treatment_group <- comp$treatment
    
    # Check if groups exist
    if (!control_group %in% available_groups) {
      stop(paste("Control group '", control_group, "' not found in data"))
    }
    
    if (!treatment_group %in% available_groups) {
      stop(paste("Treatment group '", treatment_group, "' not found in data"))
    }
    
    # Generate comparison name
    if ("name" %in% names(comp) && !is.null(comp$name)) {
      comp_name <- comp$name
    } else {
      comp_name <- paste0(treatment_group, "_vs_", control_group)
    }
    
    # Get samples for each group
    control_samples <- sample_info[[sample_col]][sample_info[[group_col]] == control_group]
    treatment_samples <- sample_info[[sample_col]][sample_info[[group_col]] == treatment_group]
    
    # Combine samples for this comparison
    selected_samples <- c(control_samples, treatment_samples)
    
    # Check if samples exist in expression data
    if (has_id_col) {
      available_samples <- setdiff(colnames(expr_matrix), id_col)
    } else {
      available_samples <- colnames(expr_matrix)
    }
    
    valid_samples <- intersect(selected_samples, available_samples)
    missing_samples <- setdiff(selected_samples, available_samples)
    
    if (length(missing_samples) > 0) {
      warning(paste("Missing samples in expression data for comparison", comp_name, ":", 
                    paste(missing_samples, collapse = ", ")))
    }
    
    if (length(valid_samples) < 2) {
      warning(paste("Not enough valid samples for comparison", comp_name, ". Skipping."))
      next
    }
    
    # Filter expression data
    if (has_id_col && include_id_col) {
      filtered_expr <- expr_matrix[, c(id_col, valid_samples), drop = FALSE]
    } else {
      filtered_expr <- expr_matrix[, valid_samples, drop = FALSE]
    }
    
    # Filter sample info
    filtered_sample_info <- sample_info[sample_info[[sample_col]] %in% valid_samples, , drop = FALSE]
    
    # Count samples
    n_control <- sum(filtered_sample_info[[group_col]] == control_group)
    n_treatment <- sum(filtered_sample_info[[group_col]] == treatment_group)
    
    # Store results
    comparison_results[[comp_name]] <- list(
      expression_data = filtered_expr,
      sample_info = filtered_sample_info,
      control_name = control_group,
      treatment_name = treatment_group,
      n_control = n_control,
      n_treatment = n_treatment
    )
    
    # Print summary
    cat("Comparison:", comp_name, "\n")
    cat("  Control (", control_group, "):", n_control, "samples\n")
    cat("  Treatment (", treatment_group, "):", n_treatment, "samples\n")
    cat("  Total samples:", length(valid_samples), "\n")
    cat("  Data dimensions:", nrow(filtered_expr), "proteins ×", ncol(filtered_expr) - (has_id_col && include_id_col), "samples\n\n")
  }
  
  cat("Generated", length(comparison_results), "comparison groups successfully.\n")
  
  return(comparison_results)
}

#' Generate all pairwise comparisons between groups
#' 
#' @description
#' Automatically generates all possible pairwise comparisons between groups
#' 
#' @param expression_data Expression data matrix
#' @param sample_info Sample information data frame
#' @param exclude_groups Groups to exclude from comparisons (optional)
#' @param control_groups Groups to always use as control when present (optional)
#' @param id_col ID column name, default "Accession"
#' @param sample_col Sample column name, default "Sample"
#' @param group_col Group column name, default "Group"
#' @param include_id_col Whether to include ID column in output, default TRUE
#' 
#' @return Named list of all pairwise comparison groups
#' 
#' @export
create_all_pairwise_comparisons <- function(expression_data,
                                            sample_info,
                                            exclude_groups = NULL,
                                            control_groups = NULL,
                                            id_col = "Accession",
                                            sample_col = "Sample", 
                                            group_col = "Group",
                                            include_id_col = TRUE) {
  
  # Get all unique groups
  all_groups <- unique(sample_info[[group_col]])
  
  # Remove excluded groups
  if (!is.null(exclude_groups)) {
    all_groups <- setdiff(all_groups, exclude_groups)
  }
  
  if (length(all_groups) < 2) {
    stop("Need at least 2 groups for pairwise comparisons")
  }
  
  cat("Generating pairwise comparisons for groups:", paste(all_groups, collapse = ", "), "\n\n")
  
  # Generate all pairwise combinations
  comparisons <- list()
  
  for (i in 1:(length(all_groups) - 1)) {
    for (j in (i + 1):length(all_groups)) {
      group1 <- all_groups[i]
      group2 <- all_groups[j]
      
      # Determine control and treatment
      if (!is.null(control_groups)) {
        if (group1 %in% control_groups && !(group2 %in% control_groups)) {
          control <- group1
          treatment <- group2
        } else if (group2 %in% control_groups && !(group1 %in% control_groups)) {
          control <- group2
          treatment <- group1
        } else {
          # Both or neither are in control_groups, use alphabetical order
          if (group1 < group2) {
            control <- group1
            treatment <- group2
          } else {
            control <- group2
            treatment <- group1
          }
        }
      } else {
        # Use alphabetical order
        if (group1 < group2) {
          control <- group1
          treatment <- group2
        } else {
          control <- group2
          treatment <- group1
        }
      }
      
      comparisons[[length(comparisons) + 1]] <- list(
        control = control,
        treatment = treatment,
        name = paste0(treatment, "_vs_", control)
      )
    }
  }
  
  # Use create_comparison_groups to generate the actual comparisons
  return(create_comparison_groups(
    expression_data = expression_data,
    sample_info = sample_info,
    comparisons = comparisons,
    id_col = id_col,
    sample_col = sample_col,
    group_col = group_col,
    include_id_col = include_id_col
  ))
}

#' Create comparison groups with specific control
#' 
#' @description
#' Creates comparisons where one group is always the control against all others
#' 
#' @param expression_data Expression data matrix
#' @param sample_info Sample information data frame
#' @param control_group Name of the control group
#' @param treatment_groups Names of treatment groups (optional, uses all others if NULL)
#' @param id_col ID column name, default "Accession"
#' @param sample_col Sample column name, default "Sample"
#' @param group_col Group column name, default "Group"
#' @param include_id_col Whether to include ID column in output, default TRUE
#' 
#' @return Named list of comparison groups
#' 
#' @export
create_control_vs_all_comparisons <- function(expression_data,
                                              sample_info,
                                              control_group,
                                              treatment_groups = NULL,
                                              id_col = "Accession",
                                              sample_col = "Sample",
                                              group_col = "Group",
                                              include_id_col = TRUE) {
  
  # Get all unique groups
  all_groups <- unique(sample_info[[group_col]])
  
  # Check if control group exists
  if (!control_group %in% all_groups) {
    stop(paste("Control group '", control_group, "' not found in data"))
  }
  
  # Determine treatment groups
  if (is.null(treatment_groups)) {
    treatment_groups <- setdiff(all_groups, control_group)
  } else {
    # Validate treatment groups
    missing_treatments <- setdiff(treatment_groups, all_groups)
    if (length(missing_treatments) > 0) {
      stop(paste("Treatment groups not found:", paste(missing_treatments, collapse = ", ")))
    }
  }
  
  if (length(treatment_groups) == 0) {
    stop("No treatment groups specified or available")
  }
  
  cat("Creating comparisons with", control_group, "as control against:", 
      paste(treatment_groups, collapse = ", "), "\n\n")
  
  # Generate comparisons
  comparisons <- lapply(treatment_groups, function(treat) {
    list(
      control = control_group,
      treatment = treat,
      name = paste0(treat, "_vs_", control_group)
    )
  })
  
  # Use create_comparison_groups to generate the actual comparisons
  return(create_comparison_groups(
    expression_data = expression_data,
    sample_info = sample_info,
    comparisons = comparisons,
    id_col = id_col,
    sample_col = sample_col,
    group_col = group_col,
    include_id_col = include_id_col
  ))
}

#' Print summary of comparison groups
#' 
#' @description
#' Prints a summary table of all comparison groups
#' 
#' @param comparison_groups List of comparison groups (output from create_comparison_groups)
#' 
#' @export
print_comparison_summary <- function(comparison_groups) {
  
  if (length(comparison_groups) == 0) {
    cat("No comparison groups found.\n")
    return(invisible())
  }
  
  # Create summary data frame
  summary_df <- data.frame(
    Comparison = names(comparison_groups),
    Control = sapply(comparison_groups, function(x) x$control_name),
    Treatment = sapply(comparison_groups, function(x) x$treatment_name),
    N_Control = sapply(comparison_groups, function(x) x$n_control),
    N_Treatment = sapply(comparison_groups, function(x) x$n_treatment),
    Total_Samples = sapply(comparison_groups, function(x) x$n_control + x$n_treatment),
    N_Proteins = sapply(comparison_groups, function(x) nrow(x$expression_data)),
    stringsAsFactors = FALSE
  )
  
  cat("=== Comparison Groups Summary ===\n\n")
  print(summary_df)
  cat("\nTotal comparisons:", nrow(summary_df), "\n")
  
  return(invisible(summary_df))
}