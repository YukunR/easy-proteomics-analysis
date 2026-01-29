[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue)](https://github.com/YukunR/easy-proteomics-analysis)
[![Gitee](https://img.shields.io/badge/Gitee-Repository-red)](https://gitee.com/yukun-r/easy-proteomics-analysis)

# Easy Proteomics Analysis

## 目录 / Table of Contents

**中文**

- [简介](#简介)
- [主要特性](#主要特性)
- [安装](#安装)
- [快速开始](#快速开始)
- [理解结果](#理解结果)
- [详细文档](#详细文档)
- [获取帮助](#获取帮助)
- [许可证](#许可证)

**English**

- [Introduction](#introduction)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Understanding Results](#understanding-results)
- [Documentation](#documentation)
- [Getting Help](#getting-help)
- [License](#license)

---

## 简介

Easy Proteomics Analysis 是一个专为蛋白质组学差异分析设计的 **R 语言工具**。它提供了从数据预处理到统计分析再到可视化的完整自动化工作流程，特别适合需要处理大规模蛋白质组学数据但缺乏编程经验的研究人员。

本工具基于配置文件运行，您只需修改 `main.R` 中的参数，无需编写任何 R 代码。同时提供了智能检查点系统，确保分析过程可以随时中断和恢复。

## 主要特性

- **🔄 检查点与恢复** - 智能检查点系统，分析可随时中断并从上次位置恢复，永不丢失进度
- **📊 完整分析流程** - 一站式解决方案：数据加载 → 缺失值插补 → 数据标准化 → PCA 分析 → 差异分析 → 火山图 → 富集分析
- **🎯 灵活阈值设定** - 自动计算最佳 FC 阈值，支持全局或逐比较设置，交互式调整
- **👥 零编程门槛** - 基于配置文件，只需修改 main.R 参数，无需编写 R 代码

## 安装

### 系统要求

- **R 语言**：R >= 4.4.0
- **RStudio**：推荐使用（可选）
- **操作系统**：Windows、macOS 或 Linux

### 安装步骤

#### 1. 安装 R 和 RStudio

**安装 R 语言：**

- 访问 R 官方网站：https://cran.r-project.org/
- 下载并安装适合您操作系统的最新版本（>= 4.4.0）

**安装 RStudio（推荐）：**

- 访问 RStudio 官方网站：https://posit.co/download/rstudio-desktop/
- 下载并安装 RStudio Desktop 免费版

> 💡 **新手提示**：RStudio 是一个强大的 R 集成开发环境，强烈推荐使用。它让 R 编程变得更加简单直观。

#### 2. 获取项目代码

**方法 1：使用 Git（推荐）**

```bash
# 从 Gitee 克隆
git clone https://gitee.com/yukun-r/easy-proteomics-analysis.git

# 或从 GitHub 克隆
git clone https://github.com/YukunR/easy-proteomics-analysis.git
```

**方法 2：直接下载**

- 访问 [Gitee](https://gitee.com/yukun-r/easy-proteomics-analysis) 或 [GitHub](https://github.com/YukunR/easy-proteomics-analysis) 项目页面
- 点击 "下载" 或 "Download ZIP"
- 解压到您选择的目录

#### 3. 打开项目

1. 启动 RStudio
2. 点击菜单 `File` → `Open Project...`
3. 浏览到项目目录，选择 `QuickProtAna.Rproj` 文件
4. 点击 "Open"

#### 4. 恢复 R 包环境

项目使用 `renv` 管理 R 包依赖。首次使用时，在 R 控制台运行：

```r
# 恢复项目依赖的 R 包
renv::restore()
```

这个过程可能需要 10-30 分钟，取决于您的网络速度。

**常见问题：**

- **下载速度慢？** 配置国内镜像：

  ```r
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  renv::restore()
  ```

- **包安装失败？** 查看 [详细文档](DOCUMENTATION.md#故障排除) 中的解决方案

## 快速开始

### 1. 准备数据文件

您需要准备两个文件：

**① 蛋白质表达数据文件** ([示例文件](data/origin_data.txt))

- 制表符分隔的文本文件
- 第一列：蛋白质 ID
- 其他列：各样本的蛋白质丰度值

```
Protein	Sample1	Sample2	Sample3	Sample4
P12345	1250.5	1180.3	850.2	920.9
Q67890	520.2	510.1	1030.5	1120.1
```

**② 样本信息文件** ([示例文件](data/sample_info.txt))

- 制表符分隔的文本文件
- 第一列：样本名称（必须与表达数据文件的列名一致）
- 第二列：分组信息

```
Sample	Group
Sample1	Control
Sample2	Control
Sample3	Treatment
Sample4	Treatment
```

> 💡 **提示**：点击上面的"示例文件"链接查看实际的数据格式

### 2. 配置参数

打开 `main.R` 文件，修改以下关键参数：

```r
# 数据文件路径
protein_expr_file <- "./data/origin_data.txt"
sample_info_file <- "./data/sample_info.txt"

# 定义比较组
comparisons <- list(
  list(control = "Control", treatment = "Treatment", name = "Treatment_vs_Control")
)

# 其他参数通常使用默认值即可
```

### 3. 运行分析

在 RStudio 中：

1. 打开 `main.R` 文件
2. 点击 `Source` 按钮（或按 `Ctrl+Shift+S`）
3. 等待分析完成

或在 R 控制台中运行：

```r
source("main.R")
```

### 4. 查看结果

所有结果保存在 `res/` 目录中：

- 火山图、PCA 图等可视化结果（PNG 格式）
- 差异分析统计结果（CSV 格式）
- 标准化和插补后的数据

## 理解结果

### 输出文件结构

```
res/
├── normalized_data.csv          # 标准化后的数据
├── pca_plot.png                 # PCA 图
└── comparison_Treatment_vs_Control/  # 比较结果目录
    ├── volcano_plot.png         # 火山图
    ├── differential_results.csv # 完整差异分析结果
    ├── significant_up.csv       # 显著上调蛋白质
    └── significant_down.csv     # 显著下调蛋白质
```

### 关键结果文件

**火山图 (volcano_plot.png)**

- 可视化展示差异表达蛋白质
- 红色点：显著上调
- 蓝色点：显著下调
- 灰色点：无显著变化

**差异分析结果 (differential_results.csv)**

- 包含所有蛋白质的统计信息
- 主要列：
  - `Protein`: 蛋白质 ID
  - `log2FC`: log2 倍数变化
  - `pvalue`: P 值
  - `adj_pvalue`: 调整后 P 值（FDR）
  - `significant`: 是否显著
  - `regulation`: Up/Down/Not significant

### 结果解读

**显著性判断标准：**

- 倍数变化（FC）> 阈值（默认自动计算）
- 调整后 P 值 < 0.05

**如何解读火山图：**

- X 轴：log2(倍数变化) - 表示表达变化的方向和幅度
- Y 轴：-log10(P 值) - 表示统计显著性
- 越靠近图的上方，统计显著性越高
- 越远离中心，表达变化越大

> 📖 更多详细说明请参阅 [详细文档](DOCUMENTATION.md#理解结果)

## 详细文档

如需了解更多信息，请参阅 [DOCUMENTATION.md](DOCUMENTATION.md)，包括：

- 详细的 R 和 RStudio 安装指南
- renv 环境恢复详细说明
- 数据文件格式要求
- 所有配置选项详解
- 完整的分析流程说明
- 高级功能使用方法
- 故障排除指南
- 常见问题解答

## 获取帮助

如果您遇到问题或有功能建议，请在以下平台提交 Issue：

- **Gitee**: https://gitee.com/yukun-r/easy-proteomics-analysis/issues
- **GitHub**: https://github.com/YukunR/easy-proteomics-analysis/issues

**提交 Issue 时，请包含：**

- 问题的详细描述
- 错误信息（如有）
- 您的 R 版本（运行 `R.version.string` 查看）
- 您的操作系统

## 许可证

本项目采用 MIT 许可证。详见 [LICENSE](LICENSE) 文件。

---

## Introduction

Easy Proteomics Analysis is an **R-based tool** designed for proteomics differential analysis. It provides a complete automated workflow from data preprocessing to statistical analysis and visualization, particularly suitable for researchers who need to process large-scale proteomics data but lack programming experience.

This tool runs based on configuration files. You only need to modify parameters in `main.R` without writing any R code. It also provides a smart checkpoint system that ensures the analysis can be interrupted and resumed at any time.

## Key Features

- **🔄 Checkpoint & Resume** - Smart checkpoint system allows analysis to be interrupted and resumed from the last position, never losing progress
- **📊 Complete Workflow** - One-stop solution: Data loading → Imputation → Normalization → PCA → Differential analysis → Volcano plots → Enrichment analysis
- **🎯 Flexible Thresholds** - Automatic FC threshold calculation, supports global or per-comparison settings, interactive adjustment
- **👥 No Coding Required** - Configuration-based, only need to modify main.R parameters, no R coding required

## Installation

### Requirements

- **R Language**: R >= 4.4.0
- **RStudio**: Recommended (optional)
- **Operating System**: Windows, macOS, or Linux

### Installation Steps

#### 1. Install R and RStudio

**Install R:**

- Visit R official website: https://cran.r-project.org/
- Download and install the latest version for your OS (>= 4.4.0)

**Install RStudio (Recommended):**

- Visit RStudio official website: https://posit.co/download/rstudio-desktop/
- Download and install RStudio Desktop free version

> 💡 **Beginner Tip**: RStudio is a powerful R integrated development environment, highly recommended. It makes R programming much simpler and more intuitive.

#### 2. Get Project Code

**Method 1: Using Git (Recommended)**

```bash
# Clone from Gitee
git clone https://gitee.com/yukun-r/easy-proteomics-analysis.git

# Or clone from GitHub
git clone https://github.com/YukunR/easy-proteomics-analysis.git
```

**Method 2: Direct Download**

- Visit [Gitee](https://gitee.com/yukun-r/easy-proteomics-analysis) or [GitHub](https://github.com/YukunR/easy-proteomics-analysis) project page
- Click "Download" or "Download ZIP"
- Extract to your chosen directory

#### 3. Open Project

1. Launch RStudio
2. Click menu `File` → `Open Project...`
3. Navigate to project directory and select `QuickProtAna.Rproj` file
4. Click "Open"

#### 4. Restore R Package Environment

The project uses `renv` to manage R package dependencies. On first use, run in R console:

```r
# Restore project R package dependencies
renv::restore()
```

This process may take 10-30 minutes depending on your network speed.

**Common Issues:**

- **Slow download?** Configure domestic mirror:

  ```r
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  renv::restore()
  ```

- **Package installation fails?** See solutions in [Detailed Documentation](DOCUMENTATION.md#troubleshooting)

## Quick Start

### 1. Prepare Data Files

You need to prepare two files:

**① Protein Expression Data File** ([Example file](data/origin_data.txt))

- Tab-separated text file
- First column: Protein IDs
- Other columns: Protein abundance values for each sample

```
Protein	Sample1	Sample2	Sample3	Sample4
P12345	1250.5	1180.3	850.2	920.9
Q67890	520.2	510.1	1030.5	1120.1
```

**② Sample Information File** ([Example file](data/sample_info.txt))

- Tab-separated text file
- First column: Sample names (must match column names in expression data file)
- Second column: Group information

```
Sample	Group
Sample1	Control
Sample2	Control
Sample3	Treatment
Sample4	Treatment
```

> 💡 **Tip**: Click the "Example file" links above to view actual data formats

### 2. Configure Parameters

Open `main.R` file and modify the following key parameters:

```r
# Data file paths
protein_expr_file <- "./data/origin_data.txt"
sample_info_file <- "./data/sample_info.txt"

# Define comparison groups
comparisons <- list(
  list(control = "Control", treatment = "Treatment", name = "Treatment_vs_Control")
)

# Other parameters can usually use default values
```

### 3. Run Analysis

In RStudio:

1. Open `main.R` file
2. Click `Source` button (or press `Ctrl+Shift+S`)
3. Wait for analysis to complete

Or run in R console:

```r
source("main.R")
```

### 4. View Results

All results are saved in `res/` directory:

- Visualization results like volcano plots, PCA plots (PNG format)
- Differential analysis statistical results (CSV format)
- Normalized and imputed data

## Understanding Results

### Output File Structure

```
res/
├── normalized_data.csv          # Normalized data
├── pca_plot.png                 # PCA plot
└── comparison_Treatment_vs_Control/  # Comparison results directory
    ├── volcano_plot.png         # Volcano plot
    ├── differential_results.csv # Complete differential analysis results
    ├── significant_up.csv       # Significantly upregulated proteins
    └── significant_down.csv     # Significantly downregulated proteins
```

### Key Result Files

**Volcano Plot (volcano_plot.png)**

- Visualizes differentially expressed proteins
- Red dots: Significantly upregulated
- Blue dots: Significantly downregulated
- Gray dots: No significant change

**Differential Analysis Results (differential_results.csv)**

- Contains statistical information for all proteins
- Main columns:
  - `Protein`: Protein ID
  - `log2FC`: log2 fold change
  - `pvalue`: P-value
  - `adj_pvalue`: Adjusted P-value (FDR)
  - `significant`: Whether significant
  - `regulation`: Up/Down/Not significant

### Result Interpretation

**Significance Criteria:**

- Fold Change (FC) > threshold (automatically calculated by default)
- Adjusted P-value < 0.05

**How to Interpret Volcano Plot:**

- X-axis: log2(Fold Change) - Indicates direction and magnitude of expression change
- Y-axis: -log10(P-value) - Indicates statistical significance
- Higher on the plot = higher statistical significance
- Further from center = larger expression change

> 📖 For more details, see [Detailed Documentation](DOCUMENTATION.md#understanding-results)

## Documentation

For more information, please refer to [DOCUMENTATION.md](DOCUMENTATION.md), including:

- Detailed R and RStudio installation guide
- renv environment restoration detailed instructions
- Data file format requirements
- All configuration options explained
- Complete analysis workflow description
- Advanced features usage
- Troubleshooting guide
- FAQ

## Getting Help

If you encounter issues or have feature suggestions, please submit an Issue on:

- **Gitee**: https://gitee.com/yukun-r/easy-proteomics-analysis/issues
- **GitHub**: https://github.com/YukunR/easy-proteomics-analysis/issues

**When submitting an Issue, please include:**

- Detailed description of the problem
- Error messages (if any)
- Your R version (run `R.version.string` to check)
- Your operating system

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) file for details.
