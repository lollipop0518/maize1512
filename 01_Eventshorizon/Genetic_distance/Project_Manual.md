# 基因型遗传距离计算与分析流程说明书

本说明书旨在帮助用户理解并使用项目中的 R 脚本，完成从 HapMap 基因型数据到遗传距离计算、结果筛选及报表生成的全流程。文档涵盖了环境准备、参数设置、脚本运行步骤及结果文件解读。

## 1. 项目背景

在作物育种和种质资源研究中，评估材料间的遗传距离（Genetic Distance, GD）是核心环节。本项目通过解析 HapMap 格式的基因型数据，计算所有样本间的成对遗传距离（基于 Manhattan 距离），并筛选出每个样本遗传距离“最近”（最相似）和“最远”（差异最大）的特定数量的种质，最终输出直观的 Excel 报表和聚类热图，辅助育种家进行亲本选配。

## 2. 环境准备与依赖安装

本项目基于 R 语言运行。请确保您的计算机已安装 R 环境（建议版本 >= 4.0.0）。

### 安装依赖包
项目依赖 `data.table`、`snpReady`、`ComplexHeatma  p` 等包。请复制以下代码块至 R 控制台或 RStudio 中运行，以完成环境配置：

```r
# 定义所需包列表
pkgs <- c("data.table", "snpReady", "openxlsx", "ComplexHeatmap", "circlize", "readxl", "writexl", "dplyr", "tidyr")

# 检查并安装 CRAN 上的包
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

# ComplexHeatmap 属于 Bioconductor 包，若上述安装失败，请使用以下命令单独安装
if (!require("ComplexHeatmap", quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
}
```

## 3. 详细操作流程

整个分析分为两个主要步骤：核心计算 (`distance.R`) 和 报表格式化 (`final_excel.R`)。

### 第一步：核心计算 (运行 `distance.R`)

该脚本负责读取原始数据，执行质量控制 (QC)，计算遗传距离矩阵，并导出初步结果和热图。

**1. 准备数据**
将您的 HapMap 格式文件（如 `365geno.hmp.txt`）放置在脚本同级目录下。

**2. 关键参数设置**
打开 `distance.R` 文件，根据您的实际需求调整以下参数：

```r
### 输入文件设置 ###
# 将此处修改为您实际的文件名
hmp_file <- "365geno.hmp.txt"

### 质量控制 (QC) 参数 (位于 snpReady::raw.data 函数内) ###
# sweep.sample: 样本过滤阈值。缺失率 > 0.8 (80%) 的样本将被剔除
# call.rate:    SNP 过滤阈值。检出率 < 0.95 (95%) 的位点将被剔除
# maf:          最小等位基因频率。MAF < 0.05 的位点将被剔除
rd <- snpReady::raw.data(
  ...
  sweep.sample = 0.80, 
  call.rate    = 0.95, 
  maf          = 0.05,
  ...
)

### 筛选邻居数量 ###
# 设定为每个样本筛选多少个最近/最远的邻居
k <- 20  
```

**3. 运行脚本与输出**
运行 `distance.R` 后，当前目录下将生成：
*   `GeneticDistance_full_365lines.xlsx`: 完整的遗传距离矩阵。
*   `GeneticDistance_nearest20.xlsx`: 最近邻居列表。
*   `GeneticDistance_farthest20.xlsx`: 最远邻居列表。
*   `GeneticDistance_heatmap.pdf`: 遗传距离聚类热图。

### 第二步：报表格式化 (运行 `final_excel.R`)

该脚本将第一步生成的标准表格转换为“一行样本名、一行距离值”的交错格式，便于人工查阅。

**⚠️ 重要操作提醒：文件路径管理**
`final_excel.R` 默认从 `GD_output` 子文件夹读取数据。为了确保脚本正常运行，请在运行此脚本前执行以下任一操作：

*   **方案 A（推荐）**：在项目文件夹下新建一个名为 `GD_output` 的文件夹，并将第一步生成的三个 `.xlsx` 文件（full, nearest, farthest）移动进去。
*   **方案 B**：修改 `final_excel.R` 脚本中的读取路径，删除 `GD_output/` 前缀。

**关键参数设置**

```r
### 选择输出模式 ###
# 默认输出“最远”邻居的详细数据。
# 如需输出“最近”邻居，请修改变量 far_or_near
far_or_near <- farthest_df  # 输出最远邻居
# far_or_near <- nearest_df   # 若需输出最近邻居，请取消此行注释并注释上一行
```

**运行结果**
运行成功后，生成最终文件：`Formatted_Genotype_Distance.xlsx`。

## 4. 结果文件解读

### 1. 遗传距离热图 (`GeneticDistance_heatmap.pdf`)
*   **颜色含义**：蓝色代表遗传距离近（相似），红色代表遗传距离远（差异大），白色为中间值。
*   **聚类树**：图表上方和左侧的树状图展示了样本间的亲缘关系结构，分支聚集在一起的样本属于同一亚群。

### 2. 格式化报表 (`Formatted_Genotype_Distance.xlsx`)
这是最终交付的详细数据表，结构如下：

| 序号 | ID | 材料.1 | 材料.2 | ... |
| :--- | :--- | :--- | :--- | :--- |
| 1 | **样本A** | 邻居X | 邻居Y | ... |
| | | 0.4521 | 0.4489 | ... |
| 2 | **样本B** | 邻居M | 邻居N | ... |
| | | 0.5102 | 0.5098 | ... |

*   **奇数行**：显示核心样本 ID 及其筛选出的 Top 20 邻居样本 ID。
*   **偶数行**：对应显示核心样本与各邻居之间的具体遗传距离数值（保留 4 位小数）。
*   **数值含义**：数值范围 0~1。0 表示完全相同，1 表示完全不同。数值越大，说明两个材料的遗传背景差异越大。

## 5. 常见问题 (FAQ)

**Q1: 运行脚本时提示 `there is no package called 'xxx'`？**
A: 请参考第 2 节“环境准备”，运行代码块中的安装命令。

**Q2: 运行 `final_excel.R` 时报错 `path does not exist`？**
A: 请检查是否已创建 `GD_output` 文件夹并将第一步生成的 Excel 文件移动进去。参考第 3 节中的“重要操作提醒”。

**Q3: 如何处理自己的数据？**
A: 将您的 HapMap 文件重命名为 `365geno.hmp.txt` 替换原文件，或者在 `distance.R` 中修改 `hmp_file` 变量为您自己的文件名。
