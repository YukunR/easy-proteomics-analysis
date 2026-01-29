# Easy Proteomics Analysis - 详细文档 / Detailed Documentation

## 目录 / Table of Contents

**中文版 / Chinese Version**

1. [项目简介](#项目简介)
2. [安装指南](#安装指南)
3. [数据准备](#数据准备)
4. [配置选项](#配置选项)
5. [分析流程](#分析流程)
6. [理解结果](#理解结果)
7. [高级功能](#高级功能)
8. [故障排除](#故障排除)
9. [常见问题](#常见问题)

**English Version**

1. [Introduction](#introduction)
2. [Installation Guide](#installation-guide)
3. [Data Preparation](#data-preparation)
4. [Configuration Options](#configuration-options)
5. [Analysis Workflow](#analysis-workflow)
6. [Understanding Results](#understanding-results)
7. [Advanced Usage](#advanced-usage)
8. [Troubleshooting](#troubleshooting)
9. [FAQ](#faq)

---

## 项目简介

Easy Proteomics Analysis 是一个专为蛋白质组学差异分析设计的 R 语言工具。它提供了从数据预处理到统计分析再到可视化的完整工作流程，特别适合需要处理大规模蛋白质组学数据但缺乏编程经验的研究人员。

**核心功能：**

- 自动化的差异蛋白质分析流程（数据加载 → 插补 → 标准化 → PCA → 差异分析 → 火山图 → 富集分析）
- 智能检查点系统，支持中断后恢复
- 灵活的阈值设定（自动计算最佳 FC 阈值或手动设置）
- 基于配置文件，只需修改 main.R 参数，无需编写 R 代码
- 完整的统计分析和可视化（火山图、PCA、热图等）
- 支持批次效应去除
- 使用 renv 管理 R 包依赖，确保环境可重现

## 安装指南

### 系统要求

- **R 语言**：R >= 4.4.0
- **RStudio**：推荐使用最新版本（可选但强烈推荐）
- **操作系统**：Windows、macOS 或 Linux
- **内存**：至少 4GB RAM（推荐 8GB 以上）
- **磁盘空间**：足够的磁盘空间用于存储结果文件和 R 包

### 安装步骤

#### 1. 安装 R 语言

**Windows 用户：**

1. 访问 R 官方网站：https://cran.r-project.org/
2. 点击 "Download R for Windows"
3. 点击 "base"
4. 下载最新版本的 R 安装程序（例如：R-4.4.0-win.exe）
5. 运行安装程序，按照提示完成安装
6. 建议使用默认安装路径

**macOS 用户：**

1. 访问 R 官方网站：https://cran.r-project.org/
2. 点击 "Download R for macOS"
3. 根据您的 macOS 版本选择合适的安装包
4. 下载并安装 R
5. 如果使用 Apple Silicon (M1/M2/M3)，请选择 ARM64 版本

**Linux 用户：**

1. 使用包管理器安装 R
2. Ubuntu/Debian: `sudo apt-get install r-base r-base-dev`
3. Fedora/CentOS: `sudo yum install R`
4. 或访问 https://cran.r-project.org/ 查看详细说明

#### 2. 安装 RStudio（推荐）

RStudio 是一个强大的 R 集成开发环境，强烈推荐使用。

1. 访问 RStudio 官方网站：https://posit.co/download/rstudio-desktop/
2. 下载适合您操作系统的 RStudio Desktop 免费版
3. 运行安装程序，按照提示完成安装
4. 安装完成后，启动 RStudio

#### 3. 获取项目代码

**方法 1：使用 Git（推荐）**

```bash
# 从 Gitee 克隆
git clone https://gitee.com/yukun-r/easy-proteomics-analysis.git

# 或从 GitHub 克隆
git clone https://github.com/YukunR/easy-proteomics-analysis.git
```

**方法 2：直接下载**

1. 访问 Gitee 或 GitHub 项目页面
2. 点击 "下载" 或 "Download ZIP"
3. 解压下载的 ZIP 文件到您选择的目录

#### 4. 打开项目

1. 启动 RStudio
2. 点击菜单 `File` → `Open Project...`
3. 浏览到项目目录，选择 `QuickProtAna.Rproj` 文件
4. 点击 "Open" 打开项目

#### 5. 使用 renv 恢复 R 包环境

本项目使用 `renv` 来管理 R 包依赖，确保所有用户使用相同版本的包。

**首次使用时：**

1. 打开项目后，RStudio 会自动加载 renv
2. 在 R 控制台中运行以下命令：

```r
# 恢复项目依赖的 R 包
renv::restore()
```

3. renv 会自动下载并安装所有需要的 R 包
4. 这个过程可能需要 10-30 分钟，取决于您的网络速度
5. 安装完成后，您会看到类似 "\* The library has been restored." 的消息

**常见问题：**

- **问题 1：renv::restore() 下载速度很慢**
  - 解决方案：配置国内镜像源

  ```r
  # 设置清华大学 CRAN 镜像
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

  # 然后再运行
  renv::restore()
  ```

- **问题 2：某些包安装失败**
  - 解决方案：检查错误信息，可能需要安装系统依赖
  - Windows 用户：确保安装了 Rtools
  - macOS 用户：确保安装了 Xcode Command Line Tools
  - Linux 用户：可能需要安装开发库（如 libxml2-dev, libcurl4-openssl-dev）

- **问题 3：Bioconductor 包安装失败**
  - 解决方案：手动安装 BiocManager
  ```r
  install.packages("BiocManager")
  BiocManager::install(version = "3.20")
  renv::restore()
  ```

#### 6. 验证安装

运行以下命令验证环境是否正确配置：

```r
# 检查 R 版本
R.version.string

# 检查关键包是否已安装
library(dplyr)
library(ggplot2)

# 如果没有错误，说明安装成功
```

## 数据准备

### 输入文件格式

本项目需要两个主要的输入文件：

#### 1. 蛋白质表达数据文件 (origin_data.txt)

这是一个制表符分隔的文本文件，包含蛋白质表达矩阵。

**文件格式：**

- **第一列**：蛋白质 ID（例如：Uniprot ID、基因名称等）
- **后续列**：各样本的蛋白质丰度值
- **列名**：样本名称（必须与 sample_info.txt 中的样本名称一致）
- **分隔符**：制表符（Tab）

**示例：**

```
Protein	Sample1	Sample2	Sample3	Sample4	Sample5	Sample6
P12345	1250.5	1180.3	1320.1	850.2	920.9	880.3
Q67890	520.2	510.1	498.8	1030.5	1120.1	1050.2
P11111	NA	850.3	820.5	450.2	480.1	470.3
Q22222	2100.5	2050.3	2180.1	1850.2	1920.9	1880.3
```

**数据要求：**

- 数值应为正数（原始丰度值或已 log 转换的值）
- 缺失值可以用 `NA`、空白或 `0` 表示
- 避免包含非数值字符（除了蛋白质 ID 列和列名）
- 建议至少有 3 个生物学重复

**查看示例文件：**
项目中提供了示例数据文件，位于 [data/origin_data.txt](data/origin_data.txt)

#### 2. 样本信息文件 (sample_info.txt)

这是一个制表符分隔的文本文件，包含样本的分组信息。

**文件格式：**

- **第一列**：样本名称（必须与 origin_data.txt 中的列名完全一致）
- **第二列**：分组信息（例如：Control、Treatment、HC、NC 等）
- **列名**：第一列为 `Sample`，第二列为 `Group`
- **分隔符**：制表符（Tab）

**示例：**

```
Sample	Group
Sample1	HC
Sample2	HC
Sample3	HC
Sample4	NC
Sample5	NC
Sample6	NC
```

**数据要求：**

- 样本名称必须与 origin_data.txt 中的列名完全匹配
- 分组名称将用于后续的差异分析比较
- 每个组建议至少有 3 个样本

**查看示例文件：**
项目中提供了示例数据文件，位于 [data/sample_info.txt](data/sample_info.txt)

### 数据质量要求

- **完整性**：确保数据已经过基本的质量控制
- **缺失值**：程序会自动处理缺失值，但缺失率过高的蛋白质会被过滤（默认阈值：60%）
- **数据范围**：确保数值在合理范围内，异常值可能影响分析结果
- **样本数量**：每组至少 3 个生物学重复，推荐 5 个以上

### 准备您自己的数据

1. **导出蛋白质表达数据**
   - 从您的蛋白质组学分析软件（如 MaxQuant、Proteome Discoverer 等）导出数据
   - 确保格式为制表符分隔的文本文件
   - 第一列为蛋白质 ID，其他列为样本丰度值

2. **创建样本信息文件**
   - 使用 Excel 或文本编辑器创建
   - 确保样本名称与表达数据文件中的列名完全一致
   - 保存为制表符分隔的文本文件

3. **放置文件**
   - 将两个文件放在项目的 `data/` 目录下
   - 或者放在任意位置，然后在 `main.R` 中指定完整路径

## 配置选项

### 配置文件说明

本项目的所有配置都在 `main.R` 文件中。您只需要修改这个文件中的参数，无需编写任何 R 代码。

### 主要配置项详解

#### 1. 基本文件路径配置

```r
base_dir <- "./res/"                          # 结果输出目录
protein_expr_file <- "./data/origin_data.txt" # 蛋白质表达数据文件
sample_info_file <- "./data/sample_info.txt"  # 样本信息文件
```

**说明：**

- `base_dir`：所有分析结果将保存在此目录下
- `protein_expr_file`：蛋白质表达矩阵文件的路径
- `sample_info_file`：样本分组信息文件的路径

#### 2. 数据预处理配置

```r
na_threshold <- 0.6  # 缺失值阈值
```

**说明：**

- 缺失值比例超过此阈值的蛋白质将被过滤掉
- 默认值 0.6 表示如果某个蛋白质在 60% 以上的样本中缺失，则该蛋白质会被移除
- 建议范围：0.5 - 0.7

#### 3. 标准化方法配置

```r
normalization_method <- "global"  # 标准化方法
use_common_proteins_for_norm <- TRUE  # 是否仅使用共同鉴定的蛋白质进行标准化
```

**标准化方法选项：**

- `"global"`：全局标准化（推荐）- 使用所有样本的中位数进行标准化
- `"within_group"`：组内标准化 - 在每个组内分别进行标准化

**use_common_proteins_for_norm：**

- `TRUE`：仅使用在所有样本中都被鉴定到的蛋白质计算标准化因子（推荐）
- `FALSE`：使用所有蛋白质计算标准化因子

#### 4. 缺失值插补方法

```r
imputation_method <- "auto"  # 插补方法
```

**插补方法选项：**

- `"auto"`：自动选择（推荐）- 根据数据特征自动选择最佳方法
- `"knn"`：K-最近邻插补 - 基于相似蛋白质的表达模式进行插补
- `"perseus"`：Perseus 风格插补 - 从正态分布中随机抽取低值进行插补

#### 5. 倍数变化（FC）阈值配置

```r
fc_threshold_mode <- "auto"      # FC 阈值模式
global_fc_threshold <- 2         # 全局 FC 阈值
```

**FC 阈值模式选项：**

- `"auto"`：自动计算（推荐）- 使用覆盖度分析自动确定最佳 FC 阈值
- `"global"`：全局阈值 - 所有比较使用相同的 `global_fc_threshold`
- `"per_comparison"`：逐比较设置 - 在每个比较中单独指定 FC 阈值

**注意：**

- 当使用 `"auto"` 模式时，`global_fc_threshold` 会被忽略
- FC 阈值表示倍数变化，例如 FC=2 表示表达量变化 2 倍

#### 6. P 值阈值配置

```r
p_threshold_mode <- "global"     # P 值阈值模式
global_p_threshold <- 0.05       # 全局 P 值阈值
```

**P 值阈值模式选项：**

- `"global"`：全局阈值（推荐）- 所有比较使用相同的 P 值阈值
- `"per_comparison"`：逐比较设置 - 在每个比较中单独指定 P 值阈值

#### 7. 差异表达分析比较组配置

```r
comparisons <- list(
  list(control = "HC", treatment = "NC", name = "HC_vs_NC"),
  list(control = "HD", treatment = "HC", name = "HD_vs_HC", fc_threshold = 2),
  list(control = c("HC", "HD"), treatment = "NC", name = "NC_vs_Rest")
)
```

**配置说明：**

- 每个比较是一个 list，包含以下字段：
  - `control`：对照组名称（必须与 sample_info.txt 中的 Group 名称一致）
  - `treatment`：处理组名称
  - `name`：比较的名称（用于结果文件命名）
  - `fc_threshold`（可选）：该比较特定的 FC 阈值（仅在 fc_threshold_mode = "per_comparison" 时有效）
  - `p_threshold`（可选）：该比较特定的 P 值阈值（仅在 p_threshold_mode = "per_comparison" 时有效）

**支持的比较类型：**

- **单组 vs 单组**：`list(control = "HC", treatment = "NC", name = "HC_vs_NC")`
- **多组 vs 单组**：`list(control = c("HC", "HD"), treatment = "NC", name = "NC_vs_Rest")`
- **单组 vs 多组**：`list(control = "HC", treatment = c("HD", "NC"), name = "Rest_vs_HC")`

**注意事项：**

- 组名必须与 sample_info.txt 中的 Group 列完全一致
- 不能在 control 和 treatment 中同时包含相同的组
- 建议使用有意义的 name，便于识别结果文件

#### 8. 富集分析配置（可选）

```r
go_background_file <- "./data/all_uniprot_go_background.csv"
kegg_background_file <- "./data/pathfromKegg_mmu.txt"
```

**说明：**

- 如果需要进行 GO 和 KEGG 富集分析，需要提供背景数据库文件
- 如果不进行富集分析，可以注释掉这些行或留空

#### 9. 自定义颜色（可选）

```r
# custom_colors <- c("#FF5733", "#33C3FF", "#75FF33")
```

**说明：**

- 如果想自定义图表颜色，取消注释并修改颜色代码
- 颜色代码使用十六进制格式（HEX）
- 如果不设置，程序会使用默认配色方案

### 配置示例

**示例 1：基本配置（推荐新手使用）**

```r
base_dir <- "./res/"
protein_expr_file <- "./data/origin_data.txt"
sample_info_file <- "./data/sample_info.txt"

na_threshold <- 0.6
normalization_method <- "global"
use_common_proteins_for_norm <- TRUE
imputation_method <- "auto"

fc_threshold_mode <- "auto"
p_threshold_mode <- "global"
global_p_threshold <- 0.05

comparisons <- list(
  list(control = "Control", treatment = "Treatment", name = "Treatment_vs_Control")
)
```

**示例 2：多组比较配置**

```r
comparisons <- list(
  list(control = "WT", treatment = "KO", name = "KO_vs_WT"),
  list(control = "WT", treatment = "HET", name = "HET_vs_WT"),
  list(control = "HET", treatment = "KO", name = "KO_vs_HET")
)
```

**示例 3：自定义阈值配置**

```r
fc_threshold_mode <- "per_comparison"
p_threshold_mode <- "per_comparison"

comparisons <- list(
  list(control = "HC", treatment = "NC", name = "HC_vs_NC",
       fc_threshold = 1.5, p_threshold = 0.05),
  list(control = "HD", treatment = "HC", name = "HD_vs_HC",
       fc_threshold = 2.0, p_threshold = 0.01)
)
```

## 分析流程

### 完整流程概述

分析流程包含以下步骤：

1. **数据加载** - 读取蛋白质表达数据和样本信息
2. **数据过滤** - 根据缺失值阈值过滤蛋白质
3. **缺失值插补** - 使用选定的方法处理缺失值
4. **批次效应去除**（可选）- 去除技术批次效应
5. **数据标准化** - 标准化样本间的系统差异
6. **PCA 分析** - 主成分分析，查看样本整体分布
7. **差异表达分析** - 对每个比较组进行统计检验
8. **阈值确定** - 自动或手动确定 FC 和 P 值阈值
9. **火山图绘制** - 可视化差异表达蛋白质
10. **富集分析**（可选）- GO 和 KEGG 富集分析
11. **结果导出** - 保存所有分析结果

### 运行分析

#### 基本运行方式

1. **在 RStudio 中打开项目**
   - 双击 `QuickProtAna.Rproj` 文件
   - 或在 RStudio 中选择 `File` → `Open Project...`

2. **配置参数**
   - 打开 `main.R` 文件
   - 根据您的数据修改配置参数（参见"配置选项"章节）

3. **运行分析**
   - 在 RStudio 中打开 `main.R`
   - 点击 `Source` 按钮运行整个脚本
   - 或者在 R 控制台中运行：

   ```r
   source("main.R")
   ```

4. **查看进度**
   - 程序会在控制台输出当前进度
   - 每个步骤完成后会显示相应信息

5. **查看结果**
   - 所有结果保存在 `res/` 目录下
   - 包括图表（PNG 格式）和数据文件（CSV 格式）

#### 使用命令行运行

您也可以在命令行中运行分析：

```r
# 运行完整分析流程
result <- run_proteomics_analysis()

# 查看项目状态
check_project_status()
```

### 检查点系统

程序在每个主要步骤后自动保存检查点。如果分析被中断，您可以从上次完成的步骤继续，而无需重新运行整个流程。

#### 检查点文件

检查点文件保存在 `res/checkpoints/` 目录中：

- `checkpoint_load_data.rds` - 数据加载后
- `checkpoint_imputation.rds` - 插补后的数据
- `checkpoint_batch_removal.rds` - 批次去除后的数据（如果启用）
- `checkpoint_normalization.rds` - 标准化后的数据
- `checkpoint_pca.rds` - PCA 分析结果
- `checkpoint_differential_*.rds` - 每个比较的差异分析结果

#### 恢复分析

如果分析被中断（例如：电脑关机、R 崩溃等）：

1. **重新打开项目**
   - 在 RStudio 中打开项目

2. **运行相同的命令**

   ```r
   source("main.R")
   # 或
   result <- run_proteomics_analysis()
   ```

3. **程序会自动检测检查点**
   - 程序会询问是否从检查点恢复
   - 输入 `y` 或 `yes` 继续之前的分析
   - 输入 `n` 或 `no` 从头开始

#### 强制重新运行特定步骤

如果您想重新运行某个特定步骤：

```r
# 重新运行标准化步骤
rerun_step("normalization")

# 重新运行多个步骤
rerun_step(c("pca", "differential_analysis"))

# 重新运行特定比较的差异分析
rerun_step("differential_analysis_HC_vs_NC")
```

#### 清理项目并重新开始

如果您想完全重新开始分析：

```r
# 清理所有检查点和结果
clean_project()

# 然后重新运行分析
result <- run_proteomics_analysis()
```

**注意：** `clean_project()` 会删除所有结果文件和检查点，请谨慎使用！

### 分析时间估计

分析所需时间取决于数据集大小和计算机性能：

- **小数据集**（< 2000 个蛋白质，< 20 个样本）：5-15 分钟
- **中等数据集**（2000-5000 个蛋白质，20-50 个样本）：15-30 分钟
- **大数据集**（> 5000 个蛋白质，> 50 个样本）：30-60 分钟

**注意：** 首次运行时，如果需要安装 R 包，可能需要额外的时间。

### 分析过程中的交互

在某些情况下，程序可能会要求您的输入：

1. **检查点恢复确认**
   - 程序检测到检查点时会询问是否恢复
   - 输入 `y` 继续，`n` 重新开始

2. **阈值设置**（如果使用 per_comparison 模式且未指定阈值）
   - 程序会询问每个比较的 FC 和 P 值阈值
   - 输入数值或按 Enter 使用默认值

3. **覆盖度分析结果**（如果使用 auto FC 阈值模式）
   - 程序会显示推荐的 FC 阈值
   - 您可以接受推荐值或输入自定义值

## 理解结果

### 输出文件说明

分析完成后，结果将保存在 `res/` 目录中。目录结构如下：

```
res/
├── checkpoints/              # 检查点文件
├── normalized_data.csv       # 标准化后的数据
├── imputed_data.csv         # 插补后的数据
├── pca_plot.png             # PCA 图
├── comparison_HC_vs_NC/     # 比较结果目录（每个比较一个目录）
│   ├── volcano_plot.png     # 火山图
│   ├── differential_results.csv  # 差异分析结果
│   ├── significant_up.csv   # 显著上调蛋白质
│   ├── significant_down.csv # 显著下调蛋白质
│   └── enrichment/          # 富集分析结果（如果启用）
│       ├── go_enrichment.csv
│       └── kegg_enrichment.csv
└── ...
```

### 主要输出文件详解

#### 1. 标准化和插补数据

**normalized_data.csv**

- 标准化后的蛋白质表达矩阵
- 可用于后续的自定义分析
- 格式与输入文件相同，但数值已标准化

**imputed_data.csv**

- 缺失值插补后的数据
- 所有缺失值已被填充

#### 2. PCA 分析结果

**pca_plot.png**

- 主成分分析图
- 展示样本在 PC1 和 PC2 上的分布
- 不同组用不同颜色标记
- 用于评估样本间的整体差异和分组效果

**如何解读：**

- 同组样本应该聚集在一起
- 不同组样本应该分开
- 如果组内样本分散或组间重叠，可能表示：
  - 生物学差异较小
  - 技术变异较大
  - 需要批次效应去除

#### 3. 差异分析结果

每个比较会生成一个独立的目录，包含以下文件：

**differential_results.csv**

- 完整的差异分析结果表
- 包含所有蛋白质的统计信息

**列说明：**

- `Protein`：蛋白质 ID
- `log2FC`：log2 倍数变化（Log2 Fold Change）
- `FC`：倍数变化（Fold Change）
- `pvalue`：原始 P 值
- `adj_pvalue`：调整后的 P 值（FDR，使用 Benjamini-Hochberg 方法）
- `significant`：是否显著（TRUE/FALSE）
- `regulation`：调控方向（Up/Down/Not significant）
- `mean_control`：对照组平均表达量
- `mean_treatment`：处理组平均表达量

**significant_up.csv**

- 显著上调的蛋白质列表
- 仅包含满足 FC 和 P 值阈值的上调蛋白质

**significant_down.csv**

- 显著下调的蛋白质列表
- 仅包含满足 FC 和 P 值阈值的下调蛋白质

#### 4. 火山图

**volcano_plot.png**

- 火山图可视化差异表达蛋白质
- X 轴：log2(Fold Change) - 表示表达变化的方向和幅度
- Y 轴：-log10(P-value) - 表示统计显著性

**图例说明：**

- **红色点**：显著上调的蛋白质（FC > 阈值 且 P < 阈值）
- **蓝色点**：显著下调的蛋白质（FC < 1/阈值 且 P < 阈值）
- **灰色点**：无显著变化的蛋白质
- **虚线**：阈值线
  - 垂直虚线：FC 阈值
  - 水平虚线：P 值阈值

**如何解读：**

- 越靠近图的上方，统计显著性越高
- 越远离中心（X=0），表达变化越大
- 右上角：显著上调且变化大的蛋白质
- 左上角：显著下调且变化大的蛋白质

#### 5. 富集分析结果（如果启用）

**go_enrichment.csv**

- GO（Gene Ontology）富集分析结果
- 包含生物学过程（BP）、细胞组分（CC）、分子功能（MF）

**kegg_enrichment.csv**

- KEGG 通路富集分析结果
- 展示显著富集的代谢通路和信号通路

**列说明：**

- `Term`：GO 术语或 KEGG 通路名称
- `Count`：该术语中的蛋白质数量
- `Pvalue`：富集的 P 值
- `adj_Pvalue`：调整后的 P 值
- `Genes`：该术语中的蛋白质列表

### 结果解读指南

#### 倍数变化（Fold Change）解读

- **FC > 2**：表达量上调 2 倍以上
- **FC = 1**：表达量无变化
- **FC < 0.5**：表达量下调 2 倍以上（即 1/2）

**log2FC 与 FC 的关系：**

- log2FC = 1 → FC = 2（上调 2 倍）
- log2FC = 0 → FC = 1（无变化）
- log2FC = -1 → FC = 0.5（下调 2 倍）

#### P 值解读

- **P < 0.05**：统计显著（常用阈值）
- **P < 0.01**：高度显著
- **P < 0.001**：极显著

**调整后 P 值（FDR）：**

- 考虑了多重检验校正
- 更保守，假阳性率更低
- 推荐使用调整后 P 值进行筛选

#### 显著性判断标准

蛋白质被认为显著差异表达需要同时满足：

1. **FC 阈值**：|log2FC| > log2(FC_threshold)
2. **P 值阈值**：adj_pvalue < p_threshold

例如，如果 FC 阈值 = 2，P 值阈值 = 0.05：

- 显著上调：FC > 2 且 adj_pvalue < 0.05
- 显著下调：FC < 0.5 且 adj_pvalue < 0.05

### 结果验证建议

1. **检查 PCA 图**
   - 确认样本分组合理
   - 识别可能的异常样本

2. **查看差异蛋白质数量**
   - 如果显著差异蛋白质过少（< 10），可能需要：
     - 降低阈值
     - 检查数据质量
     - 增加样本量
   - 如果显著差异蛋白质过多（> 50%），可能需要：
     - 提高阈值
     - 检查是否存在批次效应

3. **验证关键蛋白质**
   - 选择几个关键的差异蛋白质
   - 使用 Western Blot 或其他方法验证
   - 确认生物学意义

4. **富集分析解读**
   - 关注显著富集的通路
   - 结合生物学背景解释结果
   - 识别潜在的生物学机制

## 高级功能

### 批次效应去除

如果您的数据来自不同的实验批次（例如：不同时间点的实验、不同仪器等），可能存在批次效应。批次效应会掩盖真实的生物学差异，需要在分析前去除。

#### 如何判断是否需要批次去除？

1. **查看 PCA 图**
   - 如果样本按批次而非生物学分组聚类，说明存在批次效应
   - 例如：同一批次的不同组样本聚在一起，而不是同一组的样本聚在一起

2. **已知的批次信息**
   - 如果数据来自不同的实验批次，建议进行批次去除

#### 启用批次去除

本项目支持使用 `limma` 包的 `removeBatchEffect` 函数进行批次去除。

**配置方法：**

1. 在 `sample_info.txt` 中添加批次信息列：

```
Sample	Group	Batch
Sample1	HC	Batch1
Sample2	HC	Batch1
Sample3	HC	Batch2
Sample4	NC	Batch1
Sample5	NC	Batch2
Sample6	NC	Batch2
```

2. 在 `main.R` 中启用批次去除功能（需要修改代码）：

```r
# 在 analysis_steps.R 中找到批次去除相关代码
# 取消注释并配置批次列名
```

**注意事项：**

- 批次去除应在标准化之后、差异分析之前进行
- 确保每个批次中都有来自不同组的样本
- 批次去除可能会降低统计功效，仅在确实存在批次效应时使用

### 自定义阈值调整

#### 方法 1：使用自动 FC 阈值（推荐）

```r
fc_threshold_mode <- "auto"
```

程序会使用覆盖度分析（Coverage Analysis）自动确定最佳 FC 阈值：

- 分析不同 FC 阈值下的蛋白质覆盖度
- 选择能够最大化覆盖度同时保持合理筛选的阈值
- 程序会显示推荐的阈值，您可以接受或修改

#### 方法 2：使用全局阈值

```r
fc_threshold_mode <- "global"
global_fc_threshold <- 2  # 所有比较使用 FC = 2

p_threshold_mode <- "global"
global_p_threshold <- 0.05  # 所有比较使用 P = 0.05
```

#### 方法 3：为每个比较设置不同阈值

```r
fc_threshold_mode <- "per_comparison"
p_threshold_mode <- "per_comparison"

comparisons <- list(
  list(control = "HC", treatment = "NC", name = "HC_vs_NC",
       fc_threshold = 1.5, p_threshold = 0.05),
  list(control = "HD", treatment = "HC", name = "HD_vs_HC",
       fc_threshold = 2.0, p_threshold = 0.01)
)
```

#### 方法 4：分析后交互式调整

如果您想在看到初步结果后调整阈值：

```r
# 重新运行特定比较的差异分析，使用新阈值
rerun_step("differential_analysis_HC_vs_NC")
# 程序会询问新的阈值
```

### 多组比较

#### 合并多个组进行比较

您可以将多个组合并为一个组进行比较：

```r
comparisons <- list(
  # 将 HC 和 HD 合并为对照组，与 NC 比较
  list(control = c("HC", "HD"), treatment = "NC", name = "NC_vs_Rest"),

  # 将 HD 和 NC 合并为处理组，与 HC 比较
  list(control = "HC", treatment = c("HD", "NC"), name = "Rest_vs_HC")
)
```

**注意：**

- 合并组时，程序会计算组内所有样本的平均值
- 确保合并的组在生物学上有意义

#### 多重比较

您可以定义多个比较，程序会依次进行分析：

```r
comparisons <- list(
  list(control = "WT", treatment = "KO", name = "KO_vs_WT"),
  list(control = "WT", treatment = "HET", name = "HET_vs_WT"),
  list(control = "HET", treatment = "KO", name = "KO_vs_HET")
)
```

每个比较会生成独立的结果目录。

### 自定义可视化

#### 修改图表颜色

在 `main.R` 中设置自定义颜色：

```r
custom_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
```

颜色将按顺序应用于不同的组。

#### 修改图表样式

如果您熟悉 R 和 ggplot2，可以修改 `core/volcano_plot.R` 和 `core/pca.R` 中的绘图代码来自定义图表样式。

**常见修改：**

- 调整点的大小
- 修改字体大小
- 更改图例位置
- 添加标签

### 导出数据用于其他分析

所有中间数据都以 CSV 格式保存，您可以：

1. **在 Excel 中打开**
   - 直接双击 CSV 文件
   - 进行进一步的筛选和排序

2. **在 R 中读取**

   ```r
   # 读取标准化数据
   norm_data <- read.csv("res/normalized_data.csv", row.names = 1)

   # 读取差异分析结果
   diff_results <- read.csv("res/comparison_HC_vs_NC/differential_results.csv")
   ```

3. **在 Python 中读取**

   ```python
   import pandas as pd

   # 读取标准化数据
   norm_data = pd.read_csv("res/normalized_data.csv", index_col=0)

   # 读取差异分析结果
   diff_results = pd.read_csv("res/comparison_HC_vs_NC/differential_results.csv")
   ```

4. **导入其他工具**
   - GraphPad Prism：用于绘制出版级图表
   - Cytoscape：用于蛋白质互作网络分析
   - STRING：用于蛋白质功能分析

### 处理大数据集

如果您的数据集很大（> 10000 个蛋白质或 > 100 个样本）：

1. **增加内存限制**

   ```r
   # 在 main.R 开头添加
   options(java.parameters = "-Xmx8g")  # 分配 8GB 内存
   ```

2. **使用更快的插补方法**

   ```r
   imputation_method <- "perseus"  # 比 knn 更快
   ```

3. **减少不必要的可视化**
   - 注释掉不需要的绘图代码
   - 仅生成关键图表

4. **分批处理**
   - 将数据分成多个子集
   - 分别进行分析
   - 最后合并结果

### 重现性和版本控制

#### 记录分析参数

程序会自动保存分析参数到 `res/analysis_config.txt`，包括：

- 使用的 R 版本
- 所有配置参数
- 分析时间
- 使用的包版本

#### 使用 renv 锁定包版本

```r
# 更新 renv.lock 文件（如果您安装了新包）
renv::snapshot()

# 在其他机器上恢复相同环境
renv::restore()
```

#### Git 版本控制

建议使用 Git 管理您的分析：

```bash
# 初始化 Git 仓库
git init

# 添加文件
git add main.R data/

# 提交更改
git commit -m "Initial analysis setup"
```

**注意：** 不要将大数据文件和结果文件提交到 Git，使用 `.gitignore` 排除它们。

## 故障排除

### R 和 RStudio 相关问题

#### 1. R 版本过低

**问题：** 提示需要 R >= 4.4.0

**解决方案：**

1. 访问 https://cran.r-project.org/
2. 下载并安装最新版本的 R
3. 重新启动 RStudio

#### 2. RStudio 无法找到 R

**问题：** RStudio 启动时提示找不到 R

**解决方案（Windows）：**

1. 打开 RStudio
2. 点击 `Tools` → `Global Options`
3. 选择 `General` → `R version`
4. 点击 `Change...` 并选择正确的 R 安装路径

**解决方案（macOS/Linux）：**

```bash
# 检查 R 是否在 PATH 中
which R

# 如果没有，添加到 PATH
export PATH="/usr/local/bin:$PATH"
```

#### 3. 无法打开 .Rproj 文件

**问题：** 双击 .Rproj 文件没有反应

**解决方案：**

1. 右键点击 `QuickProtAna.Rproj`
2. 选择 "打开方式" → "RStudio"
3. 或者在 RStudio 中使用 `File` → `Open Project...`

### renv 相关问题

#### 1. renv::restore() 下载速度很慢

**问题：** 包下载速度非常慢或超时

**解决方案 1：使用国内镜像（推荐）**

```r
# 设置清华大学 CRAN 镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 设置 Bioconductor 镜像
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 然后再运行
renv::restore()
```

**解决方案 2：使用其他镜像**

```r
# 中科大镜像
options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

# 或者使用 RStudio 镜像
options(repos = c(CRAN = "https://cran.rstudio.com/"))
```

#### 2. renv::restore() 报错：包安装失败

**问题：** 某些包无法安装，提示编译错误

**解决方案（Windows）：**

1. 安装 Rtools
   - 访问 https://cran.r-project.org/bin/windows/Rtools/
   - 下载并安装与您的 R 版本匹配的 Rtools
   - 重启 RStudio

2. 验证 Rtools 安装
   ```r
   Sys.which("make")
   # 应该显示 Rtools 的路径
   ```

**解决方案（macOS）：**

1. 安装 Xcode Command Line Tools

   ```bash
   xcode-select --install
   ```

2. 如果使用 Apple Silicon (M1/M2/M3)

   ```bash
   # 安装 Homebrew
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

   # 安装必要的库
   brew install gcc
   ```

**解决方案（Linux）：**

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev

# CentOS/RHEL
sudo yum install gcc gcc-c++ make libcurl-devel openssl-devel libxml2-devel
```

#### 3. Bioconductor 包安装失败

**问题：** Bioconductor 相关包无法安装

**解决方案：**

```r
# 手动安装 BiocManager
install.packages("BiocManager")

# 设置 Bioconductor 版本
BiocManager::install(version = "3.20")

# 然后再运行 renv::restore()
renv::restore()
```

#### 4. renv 缓存问题

**问题：** renv 提示缓存损坏或不一致

**解决方案：**

```r
# 清理 renv 缓存
renv::purge()

# 重新安装所有包
renv::restore()
```

### 数据加载问题

#### 1. 无法读取数据文件

**问题：** 提示找不到文件或无法读取

**解决方案：**

1. 检查文件路径是否正确

   ```r
   # 检查文件是否存在
   file.exists("./data/origin_data.txt")
   ```

2. 检查文件编码
   - 确保文件使用 UTF-8 编码
   - 在文本编辑器中另存为 UTF-8 格式

3. 检查文件格式
   - 确保是制表符分隔（Tab-separated）
   - 不要使用逗号或空格分隔

#### 2. 样本名称不匹配

**问题：** 提示样本名称在两个文件中不一致

**解决方案：**

1. 检查 origin_data.txt 的列名
2. 检查 sample_info.txt 的 Sample 列
3. 确保两者完全一致（包括大小写、空格）

```r
# 查看数据文件的列名
data <- read.delim("./data/origin_data.txt", row.names = 1)
colnames(data)

# 查看样本信息文件的样本名
sample_info <- read.delim("./data/sample_info.txt")
sample_info$Sample
```

#### 3. 数据格式错误

**问题：** 提示数据包含非数值字符

**解决方案：**

1. 检查数据文件中是否有非数值字符（除了 NA）
2. 确保第一列是蛋白质 ID，其他列是数值
3. 缺失值应该用 `NA` 表示，不要用其他字符

### 分析运行问题

#### 1. 内存不足

**问题：** 提示内存不足或 R 崩溃

**解决方案：**

1. 关闭其他占用内存的程序
2. 增加 R 的内存限制

   ```r
   # 在 main.R 开头添加
   memory.limit(size = 16000)  # Windows，单位 MB
   ```

3. 使用更简单的插补方法

   ```r
   imputation_method <- "perseus"  # 比 knn 占用内存少
   ```

4. 减少样本数量或蛋白质数量进行测试

#### 2. 分析中断

**问题：** 分析过程中 R 崩溃或被中断

**解决方案：**

1. 重新打开项目
2. 运行相同的命令
3. 程序会自动从检查点恢复
4. 如果仍然失败，尝试清理项目重新开始
   ```r
   clean_project()
   ```

#### 3. 某个步骤一直卡住

**问题：** 程序在某个步骤长时间无响应

**解决方案：**

1. 检查 R 控制台是否有错误信息
2. 按 `Esc` 键中断当前操作
3. 检查数据是否有问题
4. 尝试使用更简单的方法
   ```r
   imputation_method <- "perseus"  # 更快
   ```

### 结果相关问题

#### 1. 没有显著差异蛋白质

**问题：** 分析完成但没有找到显著差异蛋白质

**可能原因和解决方案：**

1. **阈值太严格**
   - 降低 FC 阈值（例如从 2 改为 1.5）
   - 提高 P 值阈值（例如从 0.01 改为 0.05）

2. **生物学差异太小**
   - 检查 PCA 图，看组间是否有分离
   - 可能需要更多样本或更敏感的检测方法

3. **数据质量问题**
   - 检查原始数据质量
   - 查看缺失值比例
   - 检查是否有异常样本

#### 2. 显著差异蛋白质太多

**问题：** 几乎所有蛋白质都显著差异

**可能原因和解决方案：**

1. **阈值太宽松**
   - 提高 FC 阈值
   - 降低 P 值阈值

2. **存在批次效应**
   - 检查 PCA 图
   - 考虑进行批次去除

3. **样本分组错误**
   - 检查 sample_info.txt 中的分组是否正确

#### 3. 火山图显示异常

**问题：** 火山图形状不正常或点分布异常

**解决方案：**

1. 检查数据是否已经 log 转换
2. 检查标准化是否正确
3. 查看原始数据分布

### 包相关错误

#### 1. 找不到函数

**问题：** 提示 `could not find function "xxx"`

**解决方案：**

```r
# 检查包是否已加载
library(dplyr)
library(ggplot2)

# 如果包未安装
renv::restore()
```

#### 2. 包版本冲突

**问题：** 提示包版本不兼容

**解决方案：**

```r
# 使用 renv 恢复正确的包版本
renv::restore()

# 如果仍有问题，清理并重新安装
renv::purge()
renv::restore()
```

#### 3. 无法加载包

**问题：** `library()` 命令失败

**解决方案：**

```r
# 检查包是否已安装
installed.packages()

# 重新安装问题包
install.packages("包名")

# 或使用 renv
renv::install("包名")
```

### 其他常见问题

#### 1. 中文乱码

**问题：** 结果文件中中文显示为乱码

**解决方案：**

```r
# 设置正确的编码
Sys.setlocale("LC_ALL", "Chinese")

# 或在读取文件时指定编码
read.csv("file.csv", fileEncoding = "UTF-8")
```

#### 2. 图片无法保存

**问题：** 提示无法保存图片文件

**解决方案：**

1. 检查 res/ 目录是否存在
2. 检查是否有写入权限
3. 关闭可能打开的图片文件

#### 3. 路径问题（Windows）

**问题：** Windows 路径中的反斜杠导致错误

**解决方案：**

```r
# 使用正斜杠
protein_expr_file <- "./data/origin_data.txt"

# 或使用双反斜杠
protein_expr_file <- ".\\data\\origin_data.txt"

# 或使用 file.path()
protein_expr_file <- file.path("data", "origin_data.txt")
```

### 获取帮助

如果以上方法都无法解决您的问题：

1. **查看错误信息**
   - 仔细阅读 R 控制台中的错误信息
   - 复制完整的错误信息用于搜索

2. **检查日志文件**
   - 查看 res/ 目录下的日志文件

3. **提交 Issue**
   - Gitee: https://gitee.com/yukun-r/easy-proteomics-analysis/issues
   - GitHub: https://github.com/YukunR/easy-proteomics-analysis/issues
   - 提供详细的错误信息、R 版本、操作系统等信息

## 常见问题

**Q1: 我的数据需要预先进行 log 转换吗？**

A: 不需要。程序会自动检测数据是否已经 log 转换。如果数据范围很大（例如 0-100000），程序会自动进行 log2 转换。如果数据范围较小（例如 0-20），程序会假设已经 log 转换。

**Q2: 可以同时比较多个组吗？**

A: 可以。您可以在 `comparisons` 列表中定义多个比较，程序会依次进行分析。每个比较会生成独立的结果目录。

**Q3: 如何选择合适的插补方法？**

A:

- **auto**（推荐）：程序会根据数据特征自动选择最佳方法
- **knn**：适用于大多数情况，考虑了蛋白质间的相关性，但计算较慢
- **perseus**：快速，适用于缺失值代表低表达的情况

**Q4: 程序支持哪些数据格式？**

A: 支持制表符分隔的文本文件（.txt）。数据应该是蛋白质表达矩阵，第一列为蛋白质 ID，其他列为样本丰度值。

**Q5: 分析需要多长时间？**

A: 取决于数据集大小和计算机性能：

- 小数据集（< 2000 蛋白质）：5-15 分钟
- 中等数据集（2000-5000 蛋白质）：15-30 分钟
- 大数据集（> 5000 蛋白质）：30-60 分钟

首次运行时，如果需要安装 R 包（renv::restore()），可能需要额外 10-30 分钟。

**Q6: 可以修改可视化图片的样式吗？**

A: 可以。您可以：

1. 在 main.R 中设置 `custom_colors` 来修改颜色
2. 修改 `core/volcano_plot.R` 和 `core/pca.R` 中的绘图代码来自定义样式
3. 使用导出的数据在其他工具（如 GraphPad Prism）中重新绘图

**Q7: 如何处理技术重复和生物学重复？**

A:

- **技术重复**：建议先在数据预处理阶段取平均值，然后作为一个样本输入
- **生物学重复**：应该作为独立的样本输入，程序会在统计分析中考虑生物学变异

**Q8: FC 阈值应该设置为多少？**

A:

- 使用 `fc_threshold_mode = "auto"` 让程序自动确定（推荐）
- 常用值：1.5（较宽松）、2（常用）、3（严格）
- 取决于您的研究目的和数据特征

**Q9: 为什么我的 PCA 图显示样本没有分组？**

A: 可能的原因：

1. 组间生物学差异较小
2. 组内变异较大
3. 存在批次效应
4. 数据质量问题

建议检查原始数据质量，考虑增加样本量或进行批次去除。

**Q10: 如何引用这个工具？**

A: 如果您在研究中使用了本工具，请在论文的方法部分说明使用了本分析流程，并引用相关的 R 包（如 limma、ggplot2 等）。

**Q11: 可以分析非蛋白质组学数据吗？**

A: 可以。本工具适用于任何类似的组学数据，如：

- 代谢组学数据
- 转录组学数据（基因表达）
- 脂质组学数据

只要数据格式符合要求（表达矩阵），就可以使用本工具进行分析。

**Q12: 如何更新项目到最新版本？**

A:

```bash
# 使用 Git 更新
cd easy-proteomics-analysis
git pull

# 更新 R 包依赖
renv::restore()
```

**Q13: 分析结果可以用于发表吗？**

A: 可以。本工具使用的统计方法和可视化方法都是标准的、被广泛接受的方法。但建议：

1. 对关键结果进行实验验证（如 Western Blot）
2. 在论文中详细描述分析方法和参数
3. 考虑使用专业绘图软件优化图表质量

**Q14: 遇到问题如何获取帮助？**

A:

1. 查看本文档的"故障排除"章节
2. 在 Gitee 或 GitHub 上提交 Issue
3. 提供详细的错误信息、R 版本、操作系统等信息

---

## Introduction

Easy Proteomics Analysis is an R-based tool designed for proteomics differential analysis. It provides a complete workflow from data preprocessing to statistical analysis and visualization, particularly suitable for researchers who need to process large-scale proteomics data but lack programming experience.

**Core Features:**

- Automated differential protein analysis workflow (Data loading → Imputation → Normalization → PCA → Differential analysis → Volcano plots → Enrichment analysis)
- Smart checkpoint system with resume capability
- Flexible threshold determination (automatic FC threshold calculation or manual setting)
- Configuration-based, only need to modify main.R parameters, no R coding required
- Complete statistical analysis and visualization (volcano plots, PCA, heatmaps, etc.)
- Batch effect removal support
- Uses renv to manage R package dependencies, ensuring reproducible environments

## Installation Guide

### System Requirements

- **R Language**: R >= 4.4.0
- **RStudio**: Latest version recommended (optional but highly recommended)
- **Operating System**: Windows, macOS, or Linux
- **Memory**: At least 4GB RAM (8GB+ recommended)
- **Disk Space**: Sufficient disk space for result files and R packages

### Installation Steps

#### 1. Install R Language

**Windows Users:**

1. Visit R official website: https://cran.r-project.org/
2. Click "Download R for Windows"
3. Click "base"
4. Download the latest R installer (e.g., R-4.4.0-win.exe)
5. Run the installer and follow the prompts
6. Recommended to use default installation path

**macOS Users:**

1. Visit R official website: https://cran.r-project.org/
2. Click "Download R for macOS"
3. Choose the appropriate installer for your macOS version
4. Download and install R
5. If using Apple Silicon (M1/M2/M3), choose the ARM64 version

**Linux Users:**

1. Use package manager to install R
2. Ubuntu/Debian: `sudo apt-get install r-base r-base-dev`
3. Fedora/CentOS: `sudo yum install R`
4. Or visit https://cran.r-project.org/ for detailed instructions

#### 2. Install RStudio (Recommended)

RStudio is a powerful R integrated development environment, highly recommended.

1. Visit RStudio official website: https://posit.co/download/rstudio-desktop/
2. Download RStudio Desktop free version for your operating system
3. Run the installer and follow the prompts
4. After installation, launch RStudio

#### 3. Get Project Code

**Method 1: Using Git (Recommended)**

```bash
# Clone from Gitee
git clone https://gitee.com/yukun-r/easy-proteomics-analysis.git

# Or clone from GitHub
git clone https://github.com/YukunR/easy-proteomics-analysis.git
```

**Method 2: Direct Download**

1. Visit Gitee or GitHub project page
2. Click "Download" or "Download ZIP"
3. Extract the downloaded ZIP file to your chosen directory

#### 4. Open Project

1. Launch RStudio
2. Click menu `File` → `Open Project...`
3. Navigate to project directory and select `QuickProtAna.Rproj` file
4. Click "Open" to open the project

#### 5. Restore R Package Environment Using renv

This project uses `renv` to manage R package dependencies, ensuring all users use the same package versions.

**First Time Use:**

1. After opening the project, RStudio will automatically load renv
2. Run the following command in R console:

```r
# Restore project R package dependencies
renv::restore()
```

3. renv will automatically download and install all required R packages
4. This process may take 10-30 minutes depending on your network speed
5. After completion, you'll see a message like "\* The library has been restored."

**Common Issues:**

- **Issue 1: renv::restore() download speed is very slow**
  - Solution: Configure domestic mirror

  ```r
  # Set Tsinghua University CRAN mirror
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

  # Then run
  renv::restore()
  ```

- **Issue 2: Some packages fail to install**
  - Solution: Check error messages, may need to install system dependencies
  - Windows users: Ensure Rtools is installed
  - macOS users: Ensure Xcode Command Line Tools are installed
  - Linux users: May need to install development libraries (e.g., libxml2-dev, libcurl4-openssl-dev)

- **Issue 3: Bioconductor packages fail to install**
  - Solution: Manually install BiocManager
  ```r
  install.packages("BiocManager")
  BiocManager::install(version = "3.20")
  renv::restore()
  ```

#### 6. Verify Installation

Run the following commands to verify the environment is correctly configured:

```r
# Check R version
R.version.string

# Check if key packages are installed
library(dplyr)
library(ggplot2)

# If no errors, installation is successful
```

## Data Preparation

### Input File Format

This project requires two main input files:

#### 1. Protein Expression Data File (origin_data.txt)

This is a tab-separated text file containing the protein expression matrix.

**File Format:**

- **First column**: Protein ID (e.g., Uniprot ID, gene name, etc.)
- **Subsequent columns**: Protein abundance values for each sample
- **Column names**: Sample names (must match sample names in sample_info.txt)
- **Delimiter**: Tab

**Example:**

```
Protein	Sample1	Sample2	Sample3	Sample4	Sample5	Sample6
P12345	1250.5	1180.3	1320.1	850.2	920.9	880.3
Q67890	520.2	510.1	498.8	1030.5	1120.1	1050.2
P11111	NA	850.3	820.5	450.2	480.1	470.3
Q22222	2100.5	2050.3	2180.1	1850.2	1920.9	1880.3
```

**Data Requirements:**

- Values should be positive numbers (raw abundance or log-transformed values)
- Missing values can be represented as `NA`, blank, or `0`
- Avoid non-numeric characters (except in protein ID column and column names)
- Recommend at least 3 biological replicates

**View Example File:**
An example data file is provided in the project at [data/origin_data.txt](data/origin_data.txt)

#### 2. Sample Information File (sample_info.txt)

This is a tab-separated text file containing sample grouping information.

**File Format:**

- **First column**: Sample names (must exactly match column names in origin_data.txt)
- **Second column**: Group information (e.g., Control, Treatment, HC, NC, etc.)
- **Column names**: First column is `Sample`, second column is `Group`
- **Delimiter**: Tab

**Example:**

```
Sample	Group
Sample1	HC
Sample2	HC
Sample3	HC
Sample4	NC
Sample5	NC
Sample6	NC
```

**Data Requirements:**

- Sample names must exactly match column names in origin_data.txt
- Group names will be used for subsequent differential analysis comparisons
- Each group should have at least 3 samples

**View Example File:**
An example data file is provided in the project at [data/sample_info.txt](data/sample_info.txt)

### Data Quality Requirements

- **Completeness**: Ensure data has undergone basic quality control
- **Missing Values**: Program will automatically handle missing values, but proteins with high missing rates will be filtered (default threshold: 60%)
- **Data Range**: Ensure values are within reasonable range, outliers may affect analysis results
- **Sample Size**: At least 3 biological replicates per group, 5+ recommended

### Preparing Your Own Data

1. **Export Protein Expression Data**
   - Export data from your proteomics analysis software (e.g., MaxQuant, Proteome Discoverer, etc.)
   - Ensure format is tab-separated text file
   - First column should be protein IDs, other columns should be sample abundance values

2. **Create Sample Information File**
   - Create using Excel or text editor
   - Ensure sample names exactly match column names in expression data file
   - Save as tab-separated text file

3. **Place Files**
   - Place both files in the project's `data/` directory
   - Or place anywhere and specify full path in `main.R`

---

## Configuration Options

All configuration parameters are set in the [main.R](main.R) file. You only need to modify the values in the configuration section without writing any R code.

### Main Configuration Items

#### 1. Data File Paths

```r
base_dir <- "./res/"                          # Output directory
protein_expr_file <- "./data/origin_data.txt" # Protein expression data file
sample_info_file <- "./data/sample_info.txt"  # Sample information file
```

- **base_dir**: Directory for saving all analysis results
- **protein_expr_file**: Path to protein expression data file (tab-separated)
- **sample_info_file**: Path to sample information file (tab-separated)

#### 2. Data Preprocessing Parameters

```r
na_threshold <- 0.6                    # NA ratio threshold
normalization_method <- "global"       # Normalization method
use_common_proteins_for_norm <- TRUE   # Use common proteins for normalization
imputation_method <- "auto"            # Imputation method
```

- **na_threshold**: Proteins with NA ratio exceeding this value will be removed (default: 0.6, i.e., 60%)
- **normalization_method**:
  - `"global"`: Global normalization (recommended) - uses all samples
  - `"within_group"`: Within-group normalization - normalizes each group separately
- **use_common_proteins_for_norm**: Whether to use only commonly identified proteins for normalization (TRUE/FALSE)
- **imputation_method**:
  - `"auto"`: Automatically select method based on data characteristics (recommended)
  - `"knn"`: K-Nearest Neighbors imputation
  - `"perseus"`: Perseus-style imputation (for MNAR data)

#### 3. Threshold Settings

```r
fc_threshold_mode <- "auto"        # FC threshold mode
global_fc_threshold <- 2           # Global FC threshold
p_threshold_mode <- "global"       # P-value threshold mode
global_p_threshold <- 0.05         # Global P-value threshold
```

**FC Threshold Mode:**

- `"auto"`: Automatically calculate optimal threshold using coverage analysis (recommended)
- `"global"`: Use `global_fc_threshold` for all comparisons
- `"per_comparison"`: Set threshold for each comparison individually (interactive or specified in comparison list)

**P-value Threshold Mode:**

- `"global"`: Use `global_p_threshold` for all comparisons (default: 0.05)
- `"per_comparison"`: Set threshold for each comparison individually

#### 4. Comparison Groups

```r
comparisons <- list(
  list(control = "HC", treatment = "NC", name = "HC_vs_NC"),
  list(control = "HD", treatment = "HC", name = "HD_vs_HC", fc_threshold = 2),
  list(control = c("HC", "HD"), treatment = "NC", name = "NC_vs_Rest")
)
```

Each comparison is a list containing:

- **control**: Control group name(s) - can be a single group or multiple groups (vector)
- **treatment**: Treatment group name(s) - can be a single group or multiple groups (vector)
- **name**: Comparison name (used for output file naming)
- **fc_threshold** (optional): FC threshold for this specific comparison (only used when `fc_threshold_mode = "per_comparison"`)
- **p_threshold** (optional): P-value threshold for this specific comparison (only used when `p_threshold_mode = "per_comparison"`)

**Important Notes:**

- Group names must match those in the sample information file
- Control and treatment groups cannot have overlapping groups
- Multi-group comparisons are supported (e.g., comparing one group against multiple groups combined)

#### 5. Enrichment Analysis Files

```r
go_background_file <- "./data/all_uniprot_go_background.csv"
kegg_background_file <- "./data/pathfromKegg_mmu.txt"
```

- **go_background_file**: GO annotation background file (CSV format)
- **kegg_background_file**: KEGG pathway background file (text format)

#### 6. Custom Colors (Optional)

```r
# Uncomment to use custom colors
# custom_colors <- c("#FF5733", "#33C3FF", "#75FF33")
```

If you want to customize plot colors, uncomment this line and specify your preferred color codes.

### Configuration Examples

#### Example 1: Basic Two-Group Comparison

```r
protein_expr_file <- "./data/origin_data.txt"
sample_info_file <- "./data/sample_info.txt"

comparisons <- list(
  list(control = "Control", treatment = "Treatment", name = "Treatment_vs_Control")
)

fc_threshold_mode <- "auto"
p_threshold_mode <- "global"
global_p_threshold <- 0.05
```

#### Example 2: Multiple Comparisons with Custom Thresholds

```r
protein_expr_file <- "./data/origin_data.txt"
sample_info_file <- "./data/sample_info.txt"

comparisons <- list(
  list(control = "Control", treatment = "Treatment1", name = "T1_vs_Control", fc_threshold = 1.5),
  list(control = "Control", treatment = "Treatment2", name = "T2_vs_Control", fc_threshold = 2.0),
  list(control = "Treatment1", treatment = "Treatment2", name = "T2_vs_T1", fc_threshold = 1.5)
)

fc_threshold_mode <- "per_comparison"
p_threshold_mode <- "global"
global_p_threshold <- 0.05
```

#### Example 3: Multi-Group Comparison

```r
protein_expr_file <- "./data/origin_data.txt"
sample_info_file <- "./data/sample_info.txt"

comparisons <- list(
  list(control = c("Control", "Treatment1"), treatment = "Treatment2", name = "T2_vs_Others"),
  list(control = "Control", treatment = c("Treatment1", "Treatment2"), name = "Treatments_vs_Control")
)

fc_threshold_mode <- "auto"
p_threshold_mode <- "global"
global_p_threshold <- 0.05
```

---

## Analysis Workflow

### Complete Analysis Process

The Easy Proteomics Analysis tool provides a complete automated workflow that includes the following steps:

1. **Data Loading and Validation** - Read protein expression data and sample information, validate data format
2. **Missing Value Handling** - Filter proteins with excessive missing values, impute remaining missing values
3. **Data Separation** - Separate data into commonly identified proteins and group-specific proteins
4. **Data Normalization** - Normalize data to remove systematic differences between samples
5. **Quality Control Visualization** - Generate PCA plots to assess sample quality and grouping
6. **Differential Expression Analysis** - Perform statistical tests for each comparison group
7. **Threshold Determination** - Automatically calculate or interactively set FC and P-value thresholds
8. **Result Visualization** - Generate volcano plots, heatmaps, and other visualizations
9. **Enrichment Analysis** - Perform GO and KEGG pathway enrichment analysis (if background files provided)
10. **Result Export** - Save all analysis results and processed data
11. **Report Generation** - Generate analysis summary report

### Checkpoint System

The tool features a smart checkpoint system that automatically saves progress after each major step. This ensures:

- **Automatic Resume**: If analysis is interrupted, simply run `source("main.R")` again to resume from the last checkpoint
- **No Data Loss**: All intermediate results are saved, preventing loss of computation time
- **Flexible Rerun**: Can selectively rerun specific steps without affecting others
- **Configuration Tracking**: Automatically detects configuration changes and reruns affected steps

#### Checkpoint Files

Checkpoint files are stored in the `res/.checkpoints/` directory:

- `data_loading.rds` - Loaded and validated data
- `imputation.rds` - Data after missing value imputation
- `separation.rds` - Separated common and specific proteins
- `normalization.rds` - Normalized data
- `pca.rds` - PCA analysis results
- `differential_analysis_[comparison_name].rds` - Differential analysis results for each comparison
- `threshold_selection_[comparison_name].rds` - Threshold selection results
- `visualization_[comparison_name].rds` - Visualization results
- `enrichment_[comparison_name].rds` - Enrichment analysis results

#### Checkpoint Management Commands

```r
# Check project status
check_project_status()

# Rerun specific step
rerun_step("normalization")

# Rerun multiple steps
rerun_step(c("pca", "differential_analysis_HC_vs_NC"))

# Clean all checkpoints and restart
clean_project()
```

### Running the Analysis

#### Method 1: Run in RStudio (Recommended)

1. Open RStudio
2. Open the project file `QuickProtAna.Rproj`
3. Open `main.R` file
4. Modify configuration parameters as needed
5. Click the `Source` button (or press `Ctrl+Shift+S` / `Cmd+Shift+S`)
6. Wait for analysis to complete

#### Method 2: Run in R Console

```r
# Navigate to project directory
setwd("/path/to/EasyProteomicsAnalysis")

# Run analysis
source("main.R")
```

#### Method 3: Run from Command Line

```bash
# Navigate to project directory
cd /path/to/EasyProteomicsAnalysis

# Run with R
Rscript main.R
```

### Analysis Progress Monitoring

During analysis, the program will display:

- Current step being executed
- Progress indicators for time-consuming steps
- Warnings and important information
- Checkpoint save confirmations

Example output:

```
[2026-01-29 10:30:15] Starting proteomics analysis pipeline...
[2026-01-29 10:30:16] Step 1/11: Loading data...
[2026-01-29 10:30:18] ✓ Data loaded successfully (2500 proteins, 12 samples)
[2026-01-29 10:30:18] Step 2/11: Handling missing values...
[2026-01-29 10:30:25] ✓ Imputation completed (method: auto -> knn)
[2026-01-29 10:30:25] Checkpoint saved: imputation.rds
...
```

### Resuming Interrupted Analysis

If analysis is interrupted (power failure, manual stop, error, etc.):

1. Simply run `source("main.R")` again
2. The program will automatically detect existing checkpoints
3. Analysis will resume from the last completed step
4. No need to rerun completed steps

**Note:** If you modify configuration parameters, the program will automatically detect changes and rerun affected steps.

---

## Understanding Results

### Output Directory Structure

After analysis completion, all results are saved in the `res/` directory:

```
res/
├── .checkpoints/                    # Checkpoint files (can be deleted after analysis)
├── normalized_data.csv              # Normalized protein expression data
├── imputed_data.csv                 # Data after missing value imputation
├── common_proteins.csv              # Commonly identified proteins across all samples
├── specific_proteins.csv            # Group-specific proteins
├── pca_plot.png                     # PCA plot
├── pca_data.csv                     # PCA coordinates data
└── comparison_[name]/               # Results for each comparison
    ├── differential_results.csv     # Complete differential analysis results
    ├── significant_proteins.csv     # All significant proteins
    ├── significant_up.csv           # Significantly upregulated proteins
    ├── significant_down.csv         # Significantly downregulated proteins
    ├── volcano_plot.png             # Volcano plot
    ├── heatmap.png                  # Heatmap of significant proteins
    ├── go_enrichment_results.csv    # GO enrichment results (if available)
    ├── kegg_enrichment_results.csv  # KEGG enrichment results (if available)
    └── analysis_summary.txt         # Analysis summary report
```

### Key Result Files

#### 1. Differential Analysis Results (differential_results.csv)

Contains complete statistical analysis results for all proteins:

| Column         | Description                                    |
| -------------- | ---------------------------------------------- |
| Protein        | Protein ID                                     |
| control_mean   | Mean expression in control group               |
| treatment_mean | Mean expression in treatment group             |
| log2FC         | Log2 fold change (treatment vs control)        |
| FC             | Fold change (linear scale)                     |
| pvalue         | Original P-value from statistical test         |
| adj_pvalue     | Adjusted P-value (FDR correction)              |
| significant    | Whether protein is significant (TRUE/FALSE)    |
| regulation     | Regulation direction (Up/Down/Not significant) |

#### 2. Volcano Plot (volcano_plot.png)

Visualizes differential expression results:

- **X-axis**: log2(Fold Change) - indicates magnitude and direction of change
- **Y-axis**: -log10(P-value) - indicates statistical significance
- **Red dots**: Significantly upregulated proteins
- **Blue dots**: Significantly downregulated proteins
- **Gray dots**: Non-significant proteins
- **Dashed lines**: Significance thresholds (FC and P-value cutoffs)

#### 3. PCA Plot (pca_plot.png)

Shows overall sample relationships:

- Each point represents one sample
- Different colors represent different groups
- Closer points indicate more similar expression patterns
- Percentage values show variance explained by each principal component

#### 4. Heatmap (heatmap.png)

Displays expression patterns of significant proteins:

- Rows: Significant proteins
- Columns: Samples
- Color intensity: Expression level (red = high, blue = low)
- Hierarchical clustering shows protein and sample relationships

#### 5. Enrichment Analysis Results

**GO Enrichment (go_enrichment_results.csv):**

- GO term ID and description
- Number of significant proteins in each term
- P-value and adjusted P-value
- Enrichment ratio

**KEGG Enrichment (kegg_enrichment_results.csv):**

- KEGG pathway ID and name
- Proteins involved in each pathway
- Statistical significance
- Pathway coverage

### Result Interpretation Guide

#### Significance Criteria

A protein is considered significantly differentially expressed if:

1. **Fold Change (FC)** exceeds the threshold:
   - Upregulated: FC > threshold (e.g., FC > 2 means 2-fold increase)
   - Downregulated: FC < 1/threshold (e.g., FC < 0.5 means 2-fold decrease)

2. **Adjusted P-value** < 0.05 (default):
   - Indicates statistical significance after multiple testing correction
   - FDR (False Discovery Rate) correction is used

#### Interpreting Fold Change

- **log2FC > 0**: Protein is upregulated in treatment group
- **log2FC < 0**: Protein is downregulated in treatment group
- **log2FC = 1**: 2-fold increase (FC = 2)
- **log2FC = -1**: 2-fold decrease (FC = 0.5)
- **log2FC = 2**: 4-fold increase (FC = 4)

#### Interpreting P-values

- **P-value**: Probability that the observed difference occurred by chance
- **Adjusted P-value (adj_pvalue)**: P-value after FDR correction
  - Controls false discovery rate in multiple testing
  - More stringent than raw P-value
  - Recommended for determining significance

#### Quality Control Checks

**PCA Plot:**

- Samples from the same group should cluster together
- If samples don't cluster by group, consider:
  - Batch effects
  - Sample quality issues
  - Biological heterogeneity

**Number of Significant Proteins:**

- Too few (< 10): May indicate insufficient biological difference or too stringent thresholds
- Too many (> 50% of proteins): May indicate technical issues or too lenient thresholds
- Typical range: 5-20% of total proteins

---

## Advanced Usage

### 1. Batch Effect Removal

If your samples were processed in different batches, batch effects may confound biological differences.

**To enable batch removal:**

1. Add batch information to your sample information file:

```
Sample	Group	Batch
Sample1	Control	Batch1
Sample2	Control	Batch1
Sample3	Treatment	Batch2
Sample4	Treatment	Batch2
```

2. The program will automatically detect the `Batch` column and offer batch removal options

**Note:** Batch removal requires the `sva` package. The program will prompt you to install it if needed.

### 2. Interactive Threshold Adjustment

After automatic threshold calculation, you can interactively adjust thresholds:

- The program will display suggested thresholds based on coverage analysis
- You can accept the suggestion or enter custom values
- Volcano plots will be regenerated with new thresholds
- No need to rerun statistical analysis

**When to adjust thresholds:**

- Too many significant proteins: Increase FC threshold
- Too few significant proteins: Decrease FC threshold
- Focus on highly significant changes: Decrease P-value threshold

### 3. Multi-Group Comparisons

You can compare one group against multiple groups combined:

```r
comparisons <- list(
  # Compare Treatment against combined Control groups
  list(
    control = c("Control1", "Control2"),
    treatment = "Treatment",
    name = "Treatment_vs_Controls"
  ),

  # Compare one group against all others
  list(
    control = c("GroupA", "GroupB"),
    treatment = "GroupC",
    name = "GroupC_vs_Others"
  )
)
```

**Important:** Control and treatment groups cannot overlap.

### 4. Custom Normalization Strategy

Choose normalization method based on your data:

**Global Normalization (default):**

```r
normalization_method <- "global"
use_common_proteins_for_norm <- TRUE
```

- Uses all samples together
- Recommended for most cases
- Assumes most proteins don't change

**Within-Group Normalization:**

```r
normalization_method <- "within_group"
```

- Normalizes each group separately
- Useful when groups have very different protein compositions
- Use with caution as it may mask real biological differences

### 5. Selective Step Rerun

Rerun specific analysis steps without affecting others:

```r
# Rerun normalization only
rerun_step("normalization")

# Rerun PCA and all differential analyses
rerun_step(c("pca", "differential_analysis"))

# Rerun specific comparison
rerun_step("differential_analysis_HC_vs_NC")
```

**Available step names:**

- `data_loading`
- `imputation`
- `separation`
- `normalization`
- `pca`
- `differential_analysis` (all comparisons)
- `differential_analysis_[comparison_name]` (specific comparison)
- `threshold_selection_[comparison_name]`
- `visualization_[comparison_name]`
- `enrichment_[comparison_name]`

---

## Troubleshooting

### Common Issues and Solutions

#### 1. R or RStudio Installation Issues

**Issue:** Cannot install R or RStudio

**Solutions:**

- **Windows**: Ensure you have administrator privileges
- **macOS**: May need to allow installation from "System Preferences > Security & Privacy"
- **Linux**: Use package manager (e.g., `sudo apt-get install r-base`)
- Check system requirements: R >= 4.4.0

#### 2. renv::restore() Fails

**Issue:** Package installation fails during `renv::restore()`

**Solutions:**

**Solution 1: Configure CRAN mirror**

```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
renv::restore()
```

**Solution 2: Install system dependencies (Linux)**

```bash
# Ubuntu/Debian
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev

# CentOS/RHEL
sudo yum install libcurl-devel openssl-devel libxml2-devel
```

**Solution 3: Install packages individually**

```r
# If specific package fails, install it separately
install.packages("package_name")
renv::restore()
```

**Solution 4: Update renv**

```r
install.packages("renv")
renv::restore()
```

#### 3. Bioconductor Package Installation Fails

**Issue:** Error installing Bioconductor packages (e.g., `limma`, `impute`)

**Solutions:**

```r
# Install BiocManager first
install.packages("BiocManager")

# Install specific Bioconductor package
BiocManager::install("limma")
BiocManager::install("impute")

# Then restore renv
renv::restore()
```

#### 4. Data Loading Errors

**Issue:** "Error in read.table: cannot open file" or similar

**Solutions:**

- Check file path is correct (use forward slashes `/` or double backslashes `\\`)
- Ensure file exists in specified location
- Verify file is not open in Excel or another program
- Check file permissions (read access required)
- Use absolute path if relative path doesn't work:
  ```r
  protein_expr_file <- "C:/Users/YourName/Documents/data/origin_data.txt"
  ```

#### 5. Sample Name Mismatch

**Issue:** "Sample names in expression data don't match sample info file"

**Solutions:**

- Open both files and verify sample names are identical
- Check for extra spaces, different cases, or special characters
- Ensure sample info file includes ALL samples from expression data
- Sample names are case-sensitive

#### 6. Group Name Not Found

**Issue:** "Group 'GroupName' not found in sample information"

**Solutions:**

- Check group names in `main.R` match exactly with sample info file
- Group names are case-sensitive
- Check for typos or extra spaces
- Verify sample info file has a "Group" column

#### 7. Memory Issues

**Issue:** "Cannot allocate vector of size..." or R crashes

**Solutions:**

- Close other programs to free up memory
- Restart R session: `Session > Restart R` in RStudio
- Increase memory limit (Windows):
  ```r
  memory.limit(size = 16000)  # Set to 16GB
  ```
- Use simpler imputation method:
  ```r
  imputation_method <- "perseus"  # Instead of "knn"
  ```
- Process fewer comparisons at once

#### 8. Plotting Errors

**Issue:** Plots not generated or display incorrectly

**Solutions:**

- Update graphics device:
  ```r
  dev.off()  # Close all graphics devices
  ```
- Reinstall ggplot2:
  ```r
  install.packages("ggplot2")
  ```
- Check if output directory has write permissions
- Ensure sufficient disk space

#### 9. Checkpoint Issues

**Issue:** "Checkpoint file corrupted" or resume fails

**Solutions:**

- Clean all checkpoints and restart:
  ```r
  clean_project()
  source("main.R")
  ```
- Manually delete checkpoint directory:
  - Delete `res/.checkpoints/` folder
  - Run analysis again

#### 10. No Significant Proteins Found

**Issue:** Analysis completes but no significant proteins

**Possible Causes:**

- Biological difference is small
- Thresholds are too stringent
- Sample size is too small
- High variability within groups

**Solutions:**

- Check PCA plot - do groups separate?
- Try less stringent thresholds:
  ```r
  fc_threshold_mode <- "global"
  global_fc_threshold <- 1.5  # Instead of 2
  global_p_threshold <- 0.1   # Instead of 0.05
  ```
- Increase sample size if possible
- Check data quality and normalization

#### 11. Too Many Significant Proteins

**Issue:** Most proteins are significant (> 50%)

**Possible Causes:**

- Thresholds too lenient
- Technical artifacts or batch effects
- Normalization issues

**Solutions:**

- Use more stringent thresholds:
  ```r
  global_fc_threshold <- 3    # Instead of 2
  global_p_threshold <- 0.01  # Instead of 0.05
  ```
- Check for batch effects in PCA plot
- Verify data quality and preprocessing

#### 12. Enrichment Analysis Fails

**Issue:** GO or KEGG enrichment produces no results

**Solutions:**

- Ensure background files are provided and paths are correct
- Check if protein IDs in your data match those in background files
- Verify background file format is correct
- May need to convert protein IDs (e.g., gene names to UniProt IDs)

#### 13. RStudio Freezes or Hangs

**Issue:** RStudio becomes unresponsive during analysis

**Solutions:**

- Wait patiently - some steps (imputation, enrichment) can take several minutes
- Check R console for progress messages
- If truly frozen, restart RStudio and resume from checkpoint
- Consider running from command line instead:
  ```bash
  Rscript main.R
  ```

#### 14. Permission Denied Errors

**Issue:** "Permission denied" when saving results

**Solutions:**

- Close any files in `res/` directory that might be open
- Check folder permissions - ensure write access
- Run RStudio as administrator (Windows)
- Change output directory to a location with write permissions

#### 15. Package Version Conflicts

**Issue:** "Package 'X' version Y is required but Z is installed"

**Solutions:**

```r
# Remove package and reinstall
remove.packages("package_name")
renv::restore()

# Or update specific package
install.packages("package_name")
```

#### 16. Character Encoding Issues

**Issue:** Strange characters in protein names or sample names

**Solutions:**

- Save data files with UTF-8 encoding
- In Excel: Save As > Tools > Web Options > Encoding > UTF-8
- Or specify encoding when reading:
  ```r
  # Modify in analysis_steps.R if needed
  read.table(file, encoding = "UTF-8")
  ```

#### 17. Comparison Groups Overlap Error

**Issue:** "Control and treatment groups cannot overlap"

**Solution:**

- Ensure no group appears in both control and treatment
- Check your comparison definitions in `main.R`
- Example of INVALID comparison:
  ```r
  list(control = c("A", "B"), treatment = c("B", "C"))  # B appears in both!
  ```

#### 18. Path with Spaces Issues

**Issue:** Errors when file paths contain spaces

**Solutions:**

- Use quotes around paths with spaces:
  ```r
  protein_expr_file <- "./data/my data/origin_data.txt"
  ```
- Or avoid spaces in folder/file names
- Use underscores instead: `my_data` instead of `my data`

#### 19. Excel File Compatibility

**Issue:** Cannot read Excel files

**Solutions:**

- Save as tab-separated text (.txt) instead
- Or install readxl package:
  ```r
  install.packages("readxl")
  ```
- Ensure Excel file is not password-protected

#### 20. Analysis Takes Too Long

**Issue:** Analysis runs for hours without completing

**Solutions:**

- Check which step is running (look at console messages)
- For large datasets (> 5000 proteins), imputation can take 10-30 minutes
- Enrichment analysis can take 5-15 minutes per comparison
- Consider using faster imputation method:
  ```r
  imputation_method <- "perseus"  # Faster than knn
  ```
- Run overnight for very large datasets

### Getting Additional Help

If you encounter issues not covered here:

1. **Check Error Messages**: Read the full error message carefully - it often indicates the problem
2. **Check R Version**: Ensure R >= 4.4.0 with `R.version.string`
3. **Check Package Versions**: Run `renv::status()` to check package status
4. **Search Online**: Search the error message on Google or Stack Overflow
5. **Submit an Issue**:
   - Gitee: https://gitee.com/yukun-r/easy-proteomics-analysis/issues
   - GitHub: https://github.com/YukunR/easy-proteomics-analysis/issues

**When submitting an issue, include:**

- Full error message
- Your R version (`R.version.string`)
- Your operating system
- Steps to reproduce the error
- Relevant configuration from `main.R`

---

## FAQ

### General Questions

**Q1: Do I need to know R programming to use this tool?**

A: No. This tool is designed for users without programming experience. You only need to modify configuration values in `main.R` - no coding required.

**Q2: What file formats are supported?**

A: The tool requires tab-separated text files (.txt). You can create these from Excel by using "Save As" and selecting "Text (Tab delimited)".

**Q3: Can I use this tool on Windows/Mac/Linux?**

A: Yes. The tool works on all platforms that support R (>= 4.4.0).

**Q4: How long does the analysis take?**

A: Depends on dataset size:

- Small datasets (< 1000 proteins, < 10 samples): 2-5 minutes
- Medium datasets (1000-3000 proteins, 10-20 samples): 5-15 minutes
- Large datasets (> 3000 proteins, > 20 samples): 15-30 minutes

### Data Preparation

**Q5: What should my data look like?**

A: You need two files:

1. **Protein expression data**: Tab-separated, first column = protein IDs, other columns = sample abundance values
2. **Sample information**: Tab-separated, first column = sample names (matching expression data), second column = group names

See example files in `data/` directory.

**Q6: Should my data be log-transformed?**

A: Either is fine. The tool works with both raw abundance values and log-transformed data. The analysis will handle the transformation internally.

**Q7: How many samples do I need?**

A: Minimum 3 samples per group (recommended). With only 2 samples per group, statistical power is limited but analysis is still possible.

**Q8: Can I have missing values in my data?**

A: Yes. The tool automatically handles missing values through imputation. Proteins with too many missing values (> 60% by default) will be filtered out.

**Q9: What protein ID formats are supported?**

A: Any format is supported (UniProt IDs, gene names, etc.). Just be consistent throughout your files.

### Analysis Configuration

**Q10: How do I choose the right FC threshold?**

A: Use `fc_threshold_mode <- "auto"` (recommended). The tool will automatically calculate an optimal threshold based on your data using coverage analysis.

**Q11: What's the difference between "global" and "within_group" normalization?**

A:

- **Global**: Normalizes all samples together (recommended for most cases)
- **Within_group**: Normalizes each group separately (use only if groups have very different protein compositions)

**Q12: Which imputation method should I use?**

A: Use `imputation_method <- "auto"` (recommended). The tool will automatically select the best method based on your data characteristics.

**Q13: Can I compare more than two groups?**

A: Yes. You can define multiple pairwise comparisons in the `comparisons` list. You can also compare one group against multiple groups combined.

**Q14: Can I analyze multiple experiments in one run?**

A: Not directly. Run separate analyses for each experiment. However, you can include multiple comparisons within one experiment.

### Results Interpretation

**Q15: What does "significant" mean?**

A: A protein is significant if:

1. Fold change exceeds the threshold (e.g., FC > 2 or FC < 0.5)
2. Adjusted P-value < 0.05 (default)

Both criteria must be met.

**Q16: What's the difference between P-value and adjusted P-value?**

A:

- **P-value**: Raw statistical significance
- **Adjusted P-value**: P-value after correction for multiple testing (FDR)
- Always use adjusted P-value for determining significance

**Q17: Why are there so few/many significant proteins?**

A:

- **Too few**: Small biological difference, stringent thresholds, or high variability
- **Too many**: Large biological difference, lenient thresholds, or technical artifacts

Check PCA plot to assess overall group separation.

**Q18: How do I interpret the PCA plot?**

A:

- Each point = one sample
- Samples from the same group should cluster together
- Distance between points reflects similarity
- If groups don't separate, biological difference may be small

**Q19: What if my PCA plot shows no group separation?**

A: This suggests:

- Small biological difference between groups
- High within-group variability
- Possible batch effects or technical issues
- May still find some significant proteins, but fewer than expected

### Technical Issues

**Q20: The analysis was interrupted. Do I need to start over?**

A: No. Simply run `source("main.R")` again. The checkpoint system will automatically resume from where it stopped.

**Q21: Can I change parameters and rerun?**

A: Yes. Modify parameters in `main.R` and run again. The tool will automatically detect changes and rerun only the affected steps.

**Q22: How do I completely restart the analysis?**

A: Run `clean_project()` to delete all checkpoints and results, then run `source("main.R")`.

**Q23: Where are the results saved?**

A: All results are in the `res/` directory. Each comparison has its own subdirectory with volcano plots, differential results, and enrichment analyses.

**Q24: Can I delete the checkpoint files?**

A: Yes, after analysis is complete. Checkpoints are in `res/.checkpoints/` and can be safely deleted to save disk space.

### Advanced Usage

**Q25: Can I customize the plots?**

A: Yes. You can specify custom colors in `main.R`:

```r
custom_colors <- c("#FF5733", "#33C3FF", "#75FF33")
```

For more extensive customization, you can modify the plotting code in `analysis_steps.R` or regenerate plots using the exported data.

**Q26: How do I perform enrichment analysis?**

A: Provide GO and KEGG background files in `main.R`:

```r
go_background_file <- "./data/all_uniprot_go_background.csv"
kegg_background_file <- "./data/pathfromKegg_mmu.txt"
```

The tool will automatically perform enrichment analysis for significant proteins.

**Q27: Can I use this tool for other omics data?**

A: Yes, with caution. The tool is designed for proteomics but can work with other quantitative omics data (metabolomics, lipidomics) that have similar structure. However, some assumptions (e.g., imputation methods) are proteomics-specific.

**Q28: How do I cite this tool?**

A: Please cite the GitHub/Gitee repository:

```
Easy Proteomics Analysis. https://github.com/YukunR/easy-proteomics-analysis
```

---

## License / 许可证

**English:** This project is licensed under the MIT License. See LICENSE file for details.

**中文：** 本项目采用 MIT 许可证。详见 LICENSE 文件。

---

## Contributing / 贡献

**English:** Contributions welcome! Please submit Issues or Pull Requests via Gitee or GitHub.

**中文：** 欢迎贡献！请通过 Gitee 或 GitHub 提交 Issue 或 Pull Request。

---

**Last Updated / 最后更新**: 2026-01-29
