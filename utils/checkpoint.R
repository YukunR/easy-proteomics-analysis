# ===============================================
# ==== Checkpoint System for R Data Analysis ====
# ===============================================
# A robust checkpoint system that allows resuming analysis from interruption points
# with automatic file integrity checking and error recovery

library(jsonlite)
library(digest)

# Create workspace manager
create_workspace <- function(project_name, base_dir = "./res/") {
  workspace <- list(
    project_name = project_name,
    base_dir = base_dir,
    checkpoint_file = file.path(base_dir, "checkpoint.json"),
    log_file = file.path(base_dir, "analysis.log")
  )
  
  # Create base directory
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
    cat("Created directory:", base_dir, "\n")
  }
  
  return(workspace)
}

# Load checkpoint status
load_checkpoint <- function(workspace) {
  if (file.exists(workspace$checkpoint_file)) {
    checkpoint <- fromJSON(workspace$checkpoint_file)
    
    # Handle potential interrupted state from previous session
    if (!is.null(checkpoint$current_step)) {
      log_message(workspace, paste("Detected interrupted step from previous session:", checkpoint$current_step), "WARN")
      log_message(workspace, "This step will be restarted", "INFO")
      
      # Remove interrupted step from completed steps to ensure it runs again
      checkpoint$completed_steps <- setdiff(checkpoint$completed_steps, checkpoint$current_step)
      checkpoint$current_step <- NULL
      checkpoint$step_start_time <- NULL
      
      # Save the cleaned state
      write(toJSON(checkpoint, pretty = TRUE), workspace$checkpoint_file)
    }
    
    return(checkpoint)
  } else {
    return(list(
      completed_steps = character(0),
      current_step = NULL,
      last_update = Sys.time(),
      file_hashes = list(),
      comparison_groups_completed = character(0),
      step_start_time = NULL
    ))
  }
}

# Save checkpoint status
save_checkpoint <- function(workspace, checkpoint) {
  checkpoint$last_update <- Sys.time()
  write(toJSON(checkpoint, pretty = TRUE), workspace$checkpoint_file)
}

# Log messages with timestamp
log_message <- function(workspace, message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s", timestamp, level, message)
  cat(log_entry, "\n")
  write(log_entry, workspace$log_file, append = TRUE)
}

# Calculate file hash for integrity checking
calculate_file_hash <- function(filepath) {
  if (file.exists(filepath)) {
    return(digest(filepath, algo = "md5", file = TRUE))
  }
  return(NULL)
}

# Verify file integrity
verify_file_integrity <- function(filepath, expected_hash = NULL) {
  if (!file.exists(filepath)) {
    return(FALSE)
  }
  
  if (!is.null(expected_hash)) {
    current_hash <- calculate_file_hash(filepath)
    return(current_hash == expected_hash)
  }
  
  # Basic check: file size > 0
  return(file.size(filepath) > 0)
}

# Clean up step-related files
cleanup_step_files <- function(workspace, step_name, file_patterns = NULL) {
  log_message(workspace, paste("Cleaning files for step:", step_name))
  
  if (!is.null(file_patterns)) {
    for (pattern in file_patterns) {
      files_to_remove <- list.files(workspace$base_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
      for (file in files_to_remove) {
        if (file.exists(file)) {
          file.remove(file)
          log_message(workspace, paste("Removed file:", basename(file)))
        }
      }
    }
  }
}

# Main function to execute steps with checkpoint support
execute_step <- function(workspace, step_name, step_function, 
                         output_files = NULL, cleanup_patterns = NULL,
                         dependencies = NULL, force_rerun = FALSE, ...) {
  
  checkpoint <- load_checkpoint(workspace)
  
  # Check for interrupted step and reset if necessary
  if (!is.null(checkpoint$current_step) && checkpoint$current_step == step_name) {
    log_message(workspace, paste("Step", step_name, "was interrupted previously, resetting"), "WARN")
    checkpoint$current_step <- NULL
    checkpoint$completed_steps <- setdiff(checkpoint$completed_steps, step_name)
    save_checkpoint(workspace, checkpoint)
  }
  
  # Check if dependent steps are completed
  if (!is.null(dependencies)) {
    missing_deps <- setdiff(dependencies, checkpoint$completed_steps)
    if (length(missing_deps) > 0) {
      stop(paste("Missing dependent steps:", paste(missing_deps, collapse = ", ")))
    }
  }
  
  # Check if step is already completed and not forced to rerun
  if (step_name %in% checkpoint$completed_steps && !force_rerun) {
    # Verify output file integrity
    all_files_ok <- TRUE
    if (!is.null(output_files)) {
      for (file in output_files) {
        full_path <- file.path(workspace$base_dir, file)
        expected_hash <- checkpoint$file_hashes[[file]]
        if (!verify_file_integrity(full_path, expected_hash)) {
          log_message(workspace, paste("File", file, "corrupted or missing, need to regenerate"), "WARN")
          all_files_ok <- FALSE
          break
        }
      }
    }
    
    if (all_files_ok) {
      log_message(workspace, paste("Step", step_name, "already completed, skipping"))
      return(invisible(TRUE))
    } else {
      log_message(workspace, paste("Step", step_name, "output files have issues, re-executing"), "WARN")
      # Remove from completed list
      checkpoint$completed_steps <- setdiff(checkpoint$completed_steps, step_name)
    }
  }
  
  log_message(workspace, paste("Starting execution of step:", step_name))
  checkpoint$current_step <- step_name
  checkpoint$step_start_time <- Sys.time()
  save_checkpoint(workspace, checkpoint)
  
  # Clean up potentially incomplete files
  if (!is.null(cleanup_patterns)) {
    cleanup_step_files(workspace, step_name, cleanup_patterns)
  }
  
  # Execute step with enhanced error handling
  tryCatch({
    result <- step_function(workspace, ...)
    
    # Verify output files and record hashes
    if (!is.null(output_files)) {
      for (file in output_files) {
        full_path <- file.path(workspace$base_dir, file)
        if (file.exists(full_path)) {
          checkpoint$file_hashes[[file]] <- calculate_file_hash(full_path)
          log_message(workspace, paste("Generated file:", file))
        } else {
          log_message(workspace, paste("Warning: Expected output file not generated:", file), "WARN")
        }
      }
    }
    
    # Mark step as completed
    checkpoint$completed_steps <- unique(c(checkpoint$completed_steps, step_name))
    checkpoint$current_step <- NULL
    checkpoint$step_start_time <- NULL
    save_checkpoint(workspace, checkpoint)
    
    log_message(workspace, paste("Step", step_name, "completed successfully"))
    return(result)
    
  }, error = function(e) {
    log_message(workspace, paste("Step", step_name, "execution failed:", e$message), "ERROR")
    if (!is.null(cleanup_patterns)) {
      cleanup_step_files(workspace, step_name, cleanup_patterns)
    }
    # Clear current step on error
    checkpoint$current_step <- NULL
    checkpoint$step_start_time <- NULL
    save_checkpoint(workspace, checkpoint)
    stop(e)
  }, interrupt = function(int) {
    log_message(workspace, paste("Step", step_name, "was interrupted by user"), "WARN")
    if (!is.null(cleanup_patterns)) {
      cleanup_step_files(workspace, step_name, cleanup_patterns)
    }
    # Clear current step on interrupt
    checkpoint$current_step <- NULL
    checkpoint$step_start_time <- NULL
    save_checkpoint(workspace, checkpoint)
    stop("Analysis interrupted by user")
  })
}

# Check project status
check_project_status <- function(project_name = "proteomics_project", base_dir = "./res/") {
  workspace <- create_workspace(project_name, base_dir)
  
  if (!file.exists(workspace$checkpoint_file)) {
    cat("Project", project_name, "does not exist or has not started\n")
    return(invisible(NULL))
  }
  
  checkpoint <- fromJSON(workspace$checkpoint_file)
  
  cat("=== Project Status ===\n")
  cat("Project:", project_name, "\n")
  cat("Base directory:", base_dir, "\n")
  cat("Completed steps:", ifelse(length(checkpoint$completed_steps) == 0, "None", paste(checkpoint$completed_steps, collapse = ", ")), "\n")
  cat("Current step:", ifelse(is.null(checkpoint$current_step), "None", checkpoint$current_step), "\n")
  cat("Last update:", checkpoint$last_update, "\n")
  
  # Check for interrupted steps
  if (!is.null(checkpoint$current_step)) {
    cat("*** WARNING: Found interrupted step '", checkpoint$current_step, "' ***\n", sep = "")
    cat("This step will be restarted on next run.\n")
  }
  
  # Show file integrity status
  if (length(checkpoint$file_hashes) > 0) {
    cat("\nFile integrity check:\n")
    for (file in names(checkpoint$file_hashes)) {
      full_path <- file.path(base_dir, file)
      if (file.exists(full_path)) {
        current_hash <- calculate_file_hash(full_path)
        if (current_hash == checkpoint$file_hashes[[file]]) {
          cat("  ✓", file, "- OK\n")
        } else {
          cat("  ✗", file, "- CORRUPTED (will be regenerated)\n")
        }
      } else {
        cat("  ✗", file, "- MISSING (will be regenerated)\n")
      }
    }
  }
  
  cat("======================\n\n")
  
  # Show recent logs
  if (file.exists(workspace$log_file)) {
    cat("Recent logs (last 10 lines):\n")
    log_lines <- readLines(workspace$log_file)
    recent_logs <- tail(log_lines, 10)
    cat(paste(recent_logs, collapse = "\n"), "\n")
  }
  
  return(invisible(checkpoint))
}

# Reset interrupted step
reset_interrupted_step <- function(project_name = "proteomics_project", base_dir = "./res/") {
  workspace <- create_workspace(project_name, base_dir)
  
  if (!file.exists(workspace$checkpoint_file)) {
    cat("No checkpoint file found for project:", project_name, "\n")
    return(invisible(NULL))
  }
  
  checkpoint <- fromJSON(workspace$checkpoint_file)
  
  if (!is.null(checkpoint$current_step)) {
    interrupted_step <- checkpoint$current_step
    cat("Resetting interrupted step:", interrupted_step, "\n")
    
    # Remove from completed steps if it was there
    checkpoint$completed_steps <- setdiff(checkpoint$completed_steps, interrupted_step)
    checkpoint$current_step <- NULL
    checkpoint$step_start_time <- NULL
    
    save_checkpoint(workspace, checkpoint)
    log_message(workspace, paste("Manual reset of interrupted step:", interrupted_step))
    
    cat("Step", interrupted_step, "has been reset and will run from the beginning.\n")
  } else {
    cat("No interrupted step found.\n")
  }
}

# Verify all checkpoint files
verify_checkpoint_files <- function(project_name = "proteomics_project", base_dir = "./res/") {
  workspace <- create_workspace(project_name, base_dir)
  
  if (!file.exists(workspace$checkpoint_file)) {
    cat("No checkpoint file found for project:", project_name, "\n")
    return(invisible(FALSE))
  }
  
  checkpoint <- fromJSON(workspace$checkpoint_file)
  
  cat("=== File Verification ===\n")
  all_ok <- TRUE
  
  for (file in names(checkpoint$file_hashes)) {
    full_path <- file.path(base_dir, file)
    expected_hash <- checkpoint$file_hashes[[file]]
    
    if (!file.exists(full_path)) {
      cat("✗ MISSING:", file, "\n")
      all_ok <- FALSE
    } else {
      current_hash <- calculate_file_hash(full_path)
      if (current_hash == expected_hash) {
        cat("✓ OK:", file, "\n")
      } else {
        cat("✗ CORRUPTED:", file, "(hash mismatch)\n")
        all_ok <- FALSE
      }
    }
  }
  
  cat("========================\n")
  
  if (all_ok) {
    cat("All files are intact and valid.\n")
  } else {
    cat("Some files are missing or corrupted. They will be regenerated on next run.\n")
  }
  
  return(invisible(all_ok))
}

# Force rerun specific steps
force_rerun_steps <- function(step_names, project_name = "proteomics_project", base_dir = "./res/") {
  workspace <- create_workspace(project_name, base_dir)
  checkpoint <- load_checkpoint(workspace)
  checkpoint$completed_steps <- setdiff(checkpoint$completed_steps, step_names)
  save_checkpoint(workspace, checkpoint)
  log_message(workspace, paste("Forced rerun of steps:", paste(step_names, collapse = ", ")))
}

# Clean project (restart from beginning)
clean_project <- function(project_name = "proteomics_project", base_dir = "./res/") {
  workspace <- create_workspace(project_name, base_dir)
  if (file.exists(workspace$checkpoint_file)) {
    file.remove(workspace$checkpoint_file)
  }
  if (file.exists(workspace$log_file)) {
    file.remove(workspace$log_file)
  }
  cat("Project", project_name, "cleaned, next run will start from beginning\n")
}