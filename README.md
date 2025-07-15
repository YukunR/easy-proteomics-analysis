# Proteomics Data Analysis Pipeline

An automated proteomics data analysis tool designed for bioinformatics beginners, featuring checkpoint/resume functionality and covering the complete workflow from data normalization to differential expression analysis.

## Key Features

- **Checkpoint/Resume** - Resume analysis from interruption points without starting over
- **One-click Analysis** - Complete analysis with a single command after configuration
- **Complete Workflow** - Covers normalization, PCA, batch removal, differential expression, enrichment analysis, and GSEA
- **Batch Effect Removal** - Optional ComBat-based batch correction with interactive guidance
- **Smart Recovery** - Automatically detects interrupted states and continues from correct position
- **Detailed Logging** - Records execution status and timing for each step
- **Troubleshooting Tools** - Built-in diagnostic tools to help resolve issues

## Analysis Workflow

1. **Data Normalization and Imputation** - Missing value handling, normalization, log2 transformation
2. **Principal Component Analysis (PCA)** - Sample quality control and visualization
3. **Batch Effect Removal (Optional)** - Interactive batch assessment and correction using ComBat
4. **Comparison Group Preparation** - Automatic generation of comparison groups based on configuration
5. **Differential Expression Analysis** - t-tests, volcano plots, coverage analysis
6. **Functional Enrichment Analysis** - GO and KEGG pathway enrichment
7. **Gene Set Enrichment Analysis (GSEA)** - Pathway-level functional analysis

## Installation and Setup

### Step 1: Install R and RStudio

#### Download and Install R

1. Go to [https://cran.r-project.org/](https://cran.r-project.org/)
2. Choose your operating system (Windows/Mac/Linux)
3. Download and install the latest version of R (4.4.0 or higher recommended)

#### Download and Install RStudio (Recommended)

1. Go to [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)
2. Download RStudio Desktop (free version)
3. Install RStudio after R installation is complete

### Step 2: Set Up the Project Environment

#### Open RStudio and Install renv

```r
# Install renv package manager
install.packages("renv")
```

#### Clone or Download This Project

1. Download the project files to your computer
2. Open RStudio
3. Go to File → Open Project and select the project folder
4. Or set working directory: `setwd("path/to/your/project")`

#### Restore Package Environment

```r
# Restore all required packages
renv::restore()
```

## Quick Start

### Step 1: Prepare Your Data Files

Ensure you have the following files in the `data/` folder:

#### Sample Information File (See `data/sample_info.txt` for an example.)

Tab-separated file with sample names and groups:

```
Sample	Group
NC_1	NC
NC_2	NC
NC_3	NC
HC_1	HC
HC_2	HC
HC_3	HC
HD_1	HD
HD_2	HD
HD_3	HD
```

#### Protein Expression Data (See `data/origin_data.txt` for an example.)

Tab-separated file with protein expression values (first columns should be protein annotations, followed by sample columns)

#### Batch Information (Optional)

If you know batch information beforehand, prepare a CSV file:

```csv
Sample,Batch
NC_1,batch1
NC_2,batch1
NC_3,batch1
HC_1,batch2
HC_2,batch2
HC_3,batch2
HD_1,batch3
HD_2,batch3
HD_3,batch3
```

### Step 2: Modify Configuration Script

Modify the `Configs` section in `main.R`, then run the script with `Ctrl + Shift + Enter`.

### Step 3: Monitor Progress (Optional)

```r
# Check current status
check_status()

# View configuration
print_config()
```

### Step 4: Batch Effect Assessment

After PCA analysis, the pipeline will:

1. Display the PCA plot (PC1 vs PC2)
2. Ask if you want to perform batch removal
3. If yes, guide you through batch assignment with multiple options:
   - Enter batch for each sample individually
   - Enter all batches as comma-separated values
   - Specify that samples are already in batch order
   - Load batch information from a file

## Common Operations

### Resume After Interruption

If analysis is interrupted (power failure, manual interruption, etc.), simply rerun the same command:

```r
result <- run_proteomics_analysis()
```

### Rerun Specific Steps

```r
# Rerun normalization step
rerun_step("normalization")

# Rerun batch removal step
rerun_step("batch_removal")

# Rerun multiple steps
rerun_step(c("pca", "batch_removal", "differential_analysis_KI_vs_Old"))

# Rerun all differential analysis for all comparison groups
rerun_step(c("differential_analysis_KI_vs_Old",
            "differential_analysis_Old_vs_Young",
            "differential_analysis_KI_vs_Young"))
```

### Start Fresh

```r
# Safe cleanup (with confirmation prompt)
clean_project_safe()

# Then rerun analysis
result <- run_proteomics_analysis()
```

## Output Files Description

### Normalization Results (`[base_dir]/norm_results/`)

- `protein_abundance_data.csv` - Normalized protein abundance data
- `normalization_workflow_comparison.pdf` - Normalization workflow visualization
- `color_scheme.csv` - Color scheme for plots

### PCA Analysis Results (`[base_dir]/pca_results/`)

- `pca_biplot_PC1_PC2.pdf` - Main PCA biplot (used for batch assessment)
- `pca_biplot_PC1_PC3.pdf` - Alternative PC combination
- `pca_biplot_PC2_PC3.pdf` - Alternative PC combination
- `pca_screeplot.pdf` - Variance explained by each PC
- `pca_variance_explained.pdf` - Bar plot of variance explained
- `pca_loadings.pdf` - Top contributing variables
- `sample_correlation_heatmap.pdf` - Sample correlation heatmap
- `sample_dendrogram.pdf` - Sample clustering dendrogram

### Batch Removal Results (`[base_dir]/batch_removal_results/`)

Only created if batch removal is performed:

- `batch_assignments.csv` - Record of sample-batch assignments
- `batch_correction_pca_by_batch.pdf` - PCA comparison colored by batch
- `batch_correction_pca_by_group.pdf` - PCA comparison colored by biological group
- `batch_correction_pc1_boxplot.pdf` - PC1 distribution by batch
- `batch_correction_summary.pdf` - Combined visualization of all batch effects
- `pca_after_correction/` - Complete PCA analysis on corrected data

### Differential Expression Results (`[base_dir]/dea_results/[comparison_name]/`)

Each comparison group has its own folder containing:

- `t_test_result.csv` - Complete t-test results
- `regulated_data.csv` - Differentially expressed protein summary
- `volcano_plot.pdf` - Volcano plot
- `enrichment_results/` - Enrichment analysis results
  - `analysis_up.csv` - Upregulated protein enrichment results
  - `analysis_down.csv` - Downregulated protein enrichment results
  - `GObarplot_*.tiff` - GO enrichment bar plots
  - `KEGGdotplot_*.tiff` - KEGG enrichment dot plots
- `gsea_results/` - GSEA analysis results
  - `gsea_results.csv` - GSEA results in CSV format
  - `gsea_ridge_plot.pdf` - GSEA results visualized as a ridge plot
  - `gsea_summary_plot.pdf` - GSEA summary shown as a dot plot
  - `pathway_plots/` - Individual enrichment plots (gseaplot) for each pathway

## Configuration Parameters

### Comparison Groups Configuration

```r
comparisons <- list(
  list(control = "control_group_name", treatment = "treatment_group_name", name = "comparison_name"),
  # Add more comparisons as needed
)
```

### Parameter Descriptions

- `na_threshold`: Missing value filter threshold; proteins with missing rates above this value will be removed
- `normalization_method`:
  - `"global"`: Global median normalization
  - `"within_group"`: Within-group normalization
- `imputation_method`:
  - `"knn"`: K-nearest neighbors imputation
  - `"perseus"`: Perseus software method

## Batch Effect Removal

### When to Use Batch Removal

Consider batch removal if:

- Samples cluster by processing batch rather than biological groups in PCA
- Technical batches (e.g., different run dates, operators) are known
- Systematic differences exist between sample processing groups

### Batch Assignment Methods

1. **Interactive Entry**: Enter batch for each sample when prompted
2. **Comma-separated List**: Provide all batch assignments at once
3. **Sequential Batches**: Specify number of batches if samples are already ordered
4. **File Upload**: Load pre-prepared batch assignments from CSV

### Reference Batch Selection

- Optionally specify a reference batch for ComBat
- If not specified, ComBat will use default parameters
- Reference batch should be the most stable or control batch

## Troubleshooting

### Common Issues

**Q: Analysis restarted from beginning after interruption?**

```r
# Check project status
check_status()

# Use troubleshooting tool
troubleshoot()

# If interrupted step detected, reset status
reset_interrupted()
```

**Q: A step keeps failing?**

```r
# View detailed logs
check_status()  # Shows last 10 log lines

# Force rerun that step
rerun_step("step_name")
```

**Q: Want to change parameters and reanalyze?**

```r
# After modifying configuration parameters, clean project
clean_project_safe()

# Rerun analysis
result <- run_proteomics_analysis()
```

**Q: Batch removal failed with error?**

Common causes:

- Batches are confounded with biological groups
- Some batches have only one sample
- Missing values in critical samples

Solutions:

- Ensure each batch has multiple samples
- Check that biological groups span multiple batches
- Review sample quality before batch removal

### Diagnostic Tools

```r
# Comprehensive troubleshooting
troubleshoot()

# Check project status
check_status()

# Check batch removal status
check_batch_status()

# Verify project structure
verify_files()

# Reset interrupted state
reset_interrupted()

# View current configuration
print_config()
```

### Log Files

Detailed execution logs are saved in `[base_dir]/analysis.log`, containing start time, completion status, and error messages for each step.

## Understanding Results

### Differential Expression Results

- `log2_fold_change`: Log2 fold change; positive values indicate upregulation, negative values indicate downregulation
- `p_value`: Statistical test p-value
- `regulation_group`: Regulation direction ("up", "down", "ns")

### Enrichment Analysis Results

- `p.adjust`: Multiple testing corrected p-value
- `Count`: Number of genes enriched in this pathway
- `GeneRatio`: Ratio of enriched genes to input genes

### Batch Correction Assessment

Review the following plots to assess batch correction effectiveness:

- **PCA by Batch**: Samples should be less clustered by batch after correction
- **PCA by Group**: Biological groups should remain well-separated
- **PC1 Boxplot**: Batch-related variation in PC1 should be reduced

## Best Practices for Beginners

### Before Starting

1. Ensure your sample information file matches exactly with sample names in expression data
2. Check that group names in comparisons match those in sample_info.txt
3. Verify all required data files are in the correct locations
4. Consider batch effects if samples were processed on different days/batches

### During Analysis

1. Don't close R/RStudio while analysis is running
2. Monitor the console output for any error messages
3. If you see errors, use `troubleshoot()` for guidance
4. Carefully review PCA plots before deciding on batch removal

### After Completion

1. Check the `[base_dir]/analysis.log` file for any warnings
2. Verify that all expected output files were generated
3. Review the normalization plots before interpreting differential expression results
4. If batch removal was performed, check the comparison plots

## Getting Help

If you encounter issues:

1. Run `troubleshoot()` for diagnostic information
2. Check `[base_dir]/analysis.log` for detailed error messages
3. Verify data file formats are correct
4. Ensure all required R packages are properly installed
5. For batch removal issues, check that batches aren't confounded with groups

## Package Management with renv

This project uses `renv` for reproducible package management:

```r
# Check package status
renv::status()

# Update packages
renv::update()

# Create snapshot of current packages
renv::snapshot()
```

## License

This project is intended for academic research purposes. Please cite appropriately when using.

---

**Note**: This tool is designed for bioinformatics beginners to simplify the proteomics data analysis workflow. For more complex analyses or customized functionality, please consult with bioinformatics experts.
