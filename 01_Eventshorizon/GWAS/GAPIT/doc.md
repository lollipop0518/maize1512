# GAPIT HapMap 全流程分析脚本使用手册

## 写在前面

这是一份为植物遗传学研究人员准备的GWAS（全基因组关联分析）流程使用指南。如果你手上有表型数据和基因型数据，想找出哪些SNP位点与性状相关，但不太熟悉R语言或GWAS分析，这份文档会手把手教你完成整个流程。脚本基于GAPIT包开发，专门处理HapMap格式的基因型数据，已经过多次测试和bug修复，可以稳定运行。

## 什么是GWAS和这个脚本能做什么

GWAS是一种通过比较大量个体的基因型和表型差异来寻找性状相关基因位点的方法。比如你种了500株玉米，测量了它们的株高、产量等性状，同时对这些玉米进行了全基因组SNP分型。GWAS就是帮你找出"哪些SNP位点的变异与株高或产量显著相关"。

这个脚本把GAPIT包的复杂操作封装成了一个简单的命令行工具，你只需要准备好数据文件，运行一条命令就能完成分析。脚本会自动处理数据对齐、样本匹配、多性状分析、多模型比较等复杂操作，最后输出标准的GWAS结果文件。

## 开始之前需要准备什么

你需要在Linux、Mac或Windows系统上安装R语言（版本4.1或更高）。如果你的电脑还没有R，可以去R官网（https://www.r-project.org/）下载安装。Windows用户建议同时安装RStudio，这样操作起来更直观。

这个脚本依赖几个R包，其中最重要的是GAPIT包（用于GWAS计算）、data.table包（用于快速读写大文件）和optparse包（用于解析命令行参数）。脚本启动时会自动检查并安装data.table和optparse，但GAPIT包需要手动安装。

### 安装GAPIT包

打开R或RStudio，在命令行输入以下代码安装GAPIT。GAPIT目前主要通过GitHub分发，推荐使用remotes包安装最新版本：

```r
# 先安装remotes包（如果还没有）
install.packages("remotes")

# 从GitHub安装GAPIT
remotes::install_github("jiabowang/GAPIT")

# 安装完成后测试是否成功
library(GAPIT)
```

如果GitHub访问速度慢或者安装失败，可以尝试安装GAPIT3包作为替代：

```r
install.packages("GAPIT3")
library(GAPIT3)
```

安装过程中可能会提示安装其他依赖包，选择"Yes"或"All"安装即可。如果遇到编译错误，Windows用户需要先安装Rtools，Mac用户需要安装Xcode Command Line Tools。

## 准备你的数据文件

脚本需要两类文件：表型文件和基因型文件。这两个文件必须有共同的样本ID才能正确匹配。

### 表型文件格式

表型文件可以是CSV或TSV格式，第一行是列名，第一列必须是样本ID（列名可以是Taxa、taxa、Sample、sample、ID或id中的任意一个），后面的列是各种性状的测量值。举个例子，如果你测量了500个玉米品系的株高和产量：

```
Taxa,PlantHeight,YieldPerPlot,FloweringTime
Line001,185.5,4.2,65
Line002,192.3,4.8,68
Line003,178.9,3.9,63
Line004,201.2,5.1,70
...
```

这个文件有几个要求：样本ID列不能有重复值，性状列必须是数值类型（缺失值用NA表示或留空），样本ID不能包含特殊字符。如果你的数据在Excel里，保存成CSV格式时选择"UTF-8 CSV"可以避免中文乱码问题。

### 基因型HapMap文件格式

HapMap是植物遗传学中常用的基因型文件格式，由Tassel软件定义。标准的HapMap文件前11列是SNP的注释信息，第12列开始是各个样本的基因型。第一行是标题行，内容类似这样：

```
rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	Line001	Line002	Line003...
snp001	A/G	1	12345	+	NA	NA	NA	NA	NA	NA	AA	AG	GG...
snp002	C/T	1	23456	+	NA	NA	NA	NA	NA	NA	CC	CT	TT...
```

前11列的含义是：SNP编号、等位基因、染色体、物理位置、正负链、基因组版本、测序中心等信息（有些列可以是NA）。从第12列开始每列代表一个样本，列名是样本ID，内容是该样本在这个SNP位置的基因型（如AA、AG、GG等）。

重要提醒：HapMap文件的样本ID（第一行第12列开始的那些名字）必须与表型文件的Taxa列中的样本ID完全一致才能匹配。如果表型文件用的是"Line001"，HapMap里也必须是"Line001"，多个空格或大小写不同都会导致匹配失败。

如果你有多个染色体的HapMap文件（比如Chr1.hmp.txt、Chr2.hmp.txt等），脚本可以自动合并它们，但要确保所有文件的样本列完全一致且顺序相同。

## 运行脚本的基本方法

假设你已经把脚本文件保存为`gapit_hmp_pipeline_v2_4_fixed.r`，表型文件是`phenotype.csv`，基因型文件是`genotype.hmp.txt`，最简单的运行命令是：

```bash
Rscript gapit_hmp_pipeline_v24.r \
  --pheno phenotype.csv \
  --geno_hmp genotype.hmp.txt \
  --outdir results
```

这条命令会分析表型文件中所有的数值型性状，使用BLINK模型（GAPIT默认推荐的模型），结果保存在results文件夹里。运行过程中终端会打印实时日志，包括读取了多少样本、匹配成功多少个样本、每个性状的分析进度等。完整的日志也会保存在输出目录下的一个txt文件里。

### 多个基因型文件的合并

如果基因型数据分散在多个文件中，比如每条染色体一个文件，可以用逗号分隔多个文件路径：

```bash
Rscript gapit_hmp_pipeline_v24.r \
  --pheno phenotype.csv \
  --geno_hmp chr1.hmp.txt,chr2.hmp.txt,chr3.hmp.txt,chr4.hmp.txt,chr5.hmp.txt \
  --outdir results_all_chrs
```

脚本会自动把这些文件按行拼接（rbind）成一个大的基因型矩阵。注意所有文件的列数和样本顺序必须完全一致，否则会报错。

## 参数详细说明

脚本提供了丰富的参数选项来控制分析行为，下面逐个解释每个参数的含义和用法。

**--pheno**（必需参数）

指定表型文件的路径，支持CSV、TSV或TXT格式。文件第一列必须包含样本ID，列名可以是Taxa、taxa、Sample、sample、ID或id。如果文件路径包含空格，需要用引号括起来。

```bash
--pheno "/path/to/my phenotype data.csv"
```

**--geno_hmp**（必需参数）

指定HapMap格式的基因型文件路径。可以是单个文件或多个文件（用逗号分隔）。文件必须包含标准的HapMap标题行和前11列注释信息。

```bash
--geno_hmp genotype.hmp.txt
# 或多文件
--geno_hmp "chr1.hmp.txt,chr2.hmp.txt,chr3.hmp.txt"
```

**--outdir**

指定输出目录路径，默认是当前目录下的GAPIT_out文件夹。如果目录不存在会自动创建。建议为每次分析设置不同的输出目录名称以便区分。

```bash
--outdir my_gwas_results_2024
```

**--traits**

指定要分析哪些性状，多个性状用逗号分隔。默认值是"all"，表示自动识别表型文件中所有数值型列并全部分析。如果你只想分析其中几个性状：

```bash
--traits PlantHeight,YieldPerPlot
```

性状名称必须与表型文件中的列名完全一致（区分大小写）。如果指定的性状名在文件中不存在或不是数值类型，会被自动跳过。

**--models**

指定要运行的GAPIT模型，多个模型用逗号分隔。GAPIT支持的常用模型包括：

- GLM：广义线性模型，最简单快速但可能假阳性较高
- MLM：混合线性模型，经典方法但计算较慢
- CMLM：压缩混合线性模型，MLM的加速版本
- FarmCPU：固定和随机模型循环概率统一，速度和检验力平衡较好
- BLINK：贝叶斯信息和连锁不平衡迭代嵌套keyway，目前推荐的默认方法
- MLMM：多位点混合线性模型

默认使用BLINK模型。如果想比较多个模型的结果：

```bash
--models GLM,MLM,FarmCPU,BLINK
```

每个模型都会单独运行并输出结果文件。不同模型适用于不同的遗传结构和关联模式，多模型比较有助于提高发现的可靠性。

**--n_pcs**

指定PCA主成分数量，用于控制群体结构。默认是3。主成分分析可以捕捉样本间的亲缘关系和群体分层，将其作为协变量纳入模型可以减少假阳性。一般来说3-5个主成分足够，但如果你的群体结构特别复杂可以增加到10。

```bash
--n_pcs 5
```

**--maf**

最小等位基因频率（Minor Allele Frequency）阈值，默认0.05。低于这个频率的SNP会被过滤掉。稀有变异虽然可能有大效应，但统计检验力不足容易产生假阳性，所以通常会过滤掉。如果样本量特别大（比如超过1000个）或者你特别关注稀有变异，可以降低这个阈值：

```bash
--maf 0.01
```

**--seed**

随机数种子，默认是1。设置固定的种子可以保证每次运行结果完全一致（可重复性）。如果你想测试结果的稳定性，可以尝试不同的种子值：

```bash
--seed 12345
```

**--hmp_header_cols**

HapMap文件前导列数，默认是11（标准HapMap格式）。如果你的HapMap文件经过修改，前面的注释列不是11列，需要手动指定正确的数字。比如有些简化的HapMap文件可能只有前5列：

```bash
--hmp_header_cols 5
```

## 一个完整的运行示例

假设你在做玉米GWAS分析，有500个自交系的数据：

- 表型文件：maize_traits.csv，包含株高、产量、开花时间三个性状
- 基因型文件：10条染色体的HapMap文件（maize_chr1.hmp.txt到maize_chr10.hmp.txt）
- 想用FarmCPU和BLINK两个模型进行分析
- MAF阈值设为0.02（因为样本量较大）
- 使用5个主成分控制群体结构

完整的命令如下：

```bash
Rscript gapit_hmp_pipeline_v24.r \
  --pheno maize_traits.csv \
  --geno_hmp maize_chr1.hmp.txt,maize_chr2.hmp.txt,maize_chr3.hmp.txt,maize_chr4.hmp.txt,maize_chr5.hmp.txt,maize_chr6.hmp.txt,maize_chr7.hmp.txt,maize_chr8.hmp.txt,maize_chr9.hmp.txt,maize_chr10.hmp.txt \
  --outdir maize_gwas_results \
  --traits PlantHeight,YieldPerPlot,FloweringTime \
  --models FarmCPU,BLINK \
  --n_pcs 5 \
  --maf 0.02 \
  --seed 2024
```

运行开始后，终端会显示类似这样的日志：

```
[2024-12-17 10:30:15] [INFO] ===== GAPIT HapMap 流程启动 (精简稳定版 v2.4) =====
[2024-12-17 10:30:15] [INFO] 输出目录: /home/user/maize_gwas_results
[2024-12-17 10:30:15] [INFO] 使用包: GAPIT
[2024-12-17 10:30:16] [INFO] 表型: 样本数=500, 性状数=3
[2024-12-17 10:30:16] [INFO] HapMap文件数: 10
[2024-12-17 10:30:16] [INFO]   读入: maize_chr1.hmp.txt
[2024-12-17 10:30:18] [INFO]   读入: maize_chr2.hmp.txt
...
[2024-12-17 10:31:05] [INFO] HapMap前导列数: 11
[2024-12-17 10:31:05] [INFO] 检测到HapMap样本数: 500
[2024-12-17 10:31:05] [INFO] 样本交集数: 498 (表型=500, HapMap=500)
[2024-12-17 10:31:05] [INFO] 样本对齐完成: 最终样本数=498
[2024-12-17 10:31:05] [INFO] 运行模型: FarmCPU, BLINK
[2024-12-17 10:31:05] [INFO] >>> 开始分析性状: PlantHeight
...
```

整个流程可能需要几分钟到几小时不等，取决于数据规模和计算机性能。SNP数量越多、样本数越多、性状和模型越多，耗时越长。

## 理解输出结果

分析完成后，输出目录会包含以下内容：

```
maize_gwas_results/
├── run_log_20241217_103015.txt          # 完整运行日志
├── GWAS_AllTraits_AllModels_Combined.csv # 所有性状所有模型的汇总结果
├── trait_PlantHeight/                    # 株高性状的结果目录
│   ├── GWAS_Combined_PlantHeight.csv     # 该性状所有模型的合并结果
│   ├── GAPIT.Association.GWAS_Results.FarmCPU.Trait.csv  # FarmCPU模型结果
│   ├── GAPIT.Association.GWAS_Results.BLINK.Trait.csv    # BLINK模型结果
│   ├── GAPIT.Manhattan.FarmCPU.Trait.pdf                 # Manhattan图
│   ├── GAPIT.QQ-Plot.FarmCPU.Trait.pdf                   # QQ图
│   └── ... (其他GAPIT输出文件)
├── trait_YieldPerPlot/
└── trait_FloweringTime/
```

### 核心结果文件解读

最重要的结果是各个CSV文件，包含每个SNP的关联统计信息。打开`GWAS_Combined_PlantHeight.csv`，你会看到类似这样的列：

- **SNP**：SNP标记ID
- **CHR**：染色体编号
- **BP**：物理位置（碱基对坐标）
- **P**：P值，表示该SNP与性状关联的显著性，越小表示越显著
- **FDR_BH**：经过Benjamini-Hochberg多重检验校正后的FDR（假发现率），控制假阳性用
- **TRAIT**：性状名称
- **SOURCE_FILE**：来源于哪个模型的结果文件

GWAS分析的核心目标是找出P值足够小的SNP。传统上使用P < 5e-8作为全基因组显著性阈值（对应Bonferroni校正），但这个标准在植物中往往过于严格。实际应用中可以根据QQ图和Manhattan图来设定阈值，或者使用FDR < 0.05作为标准。

汇总文件`GWAS_AllTraits_AllModels_Combined.csv`包含了所有性状所有模型的结果，方便你一次性筛选跨性状或跨模型重复检测到的显著SNP，这些往往是最可靠的候选位点。

### 可视化结果

每个性状目录下的PDF文件提供了标准的GWAS可视化图表：

**Manhattan图**（曼哈顿图）：横轴是SNP在染色体上的位置，纵轴是-log10(P值)，每个点代表一个SNP。显著关联的SNP会形成高耸的"摩天大楼"，正常的SNP则像曼哈顿天际线一样低矮分散。这张图能直观展示全基因组范围内哪些区域与性状相关。

**QQ图**（分位数-分位数图）：横轴是期望的P值分布，纵轴是实际观测到的P值分布。如果模型控制得好，大部分点应该落在对角线上，只有最显著的少数点偏离对角线。如果大量点在对角线上方，说明可能存在未控制的群体结构（假阳性膨胀）；如果点在对角线下方，可能存在过度校正问题。

这些图表可以用来评估分析质量和选择合适的显著性阈值。

## 常见问题排查

**问题：提示"HapMap与表型样本交集过小，无法继续分析"**

这是最常见的错误，说明样本ID匹配失败。日志文件会打印出表型和HapMap的前5个样本名，仔细对比它们的差异。常见原因包括：前后有空格（"Line001"和" Line001 "被认为是不同的）、大小写不同（"Line001"和"line001"）、命名格式不一致（表型用"L001"但HapMap用"Line001"）、有不可见字符（从Excel复制粘贴时容易出现）。解决方法是统一两个文件的样本ID格式，可以在R里用`trimws()`和`tolower()`函数批量处理。

**问题：某个性状没有输出结果**

日志里会说明原因。如果提示"insufficient_samples"，说明该性状有效样本太少（剔除NA值后不足10个），需要检查表型数据质量。如果提示"gapit_error"，可能是该性状方差太小（所有值都差不多）或者存在极端异常值导致模型崩溃，可以尝试对数据进行对数转换或标准化。

**问题：分析运行极其缓慢**

如果SNP数量超过100万且使用MLM或CMLM模型，计算亲缘关系矩阵会非常耗时。建议先用BLINK或FarmCPU模型，它们针对大数据做了优化。另外可以考虑提高MAF阈值来减少SNP数量，比如从0.05提高到0.1。如果服务器内存不足，可能需要分染色体单独分析。

**问题：输出的P值全部都很大，没有显著关联**

这有几种可能：一是性状确实是受微效多基因控制，单个SNP效应很小；二是样本量不足（一般至少需要200个以上才有合理的统计检验力）；三是基因型质量差或者测序深度低导致噪音太多。可以尝试降低MAF阈值、增加样本量或者检查基因型数据质量。

**问题：运行到一半中断或报内存错误**

大型数据集（比如500个样本×100万SNP）会占用大量内存。可以尝试：只分析部分染色体、提高MAF阈值减少SNP数量、关闭其他占内存的程序、或者在配置更高的服务器上运行。Windows系统在处理大文件时更容易出现内存问题，建议使用Linux系统。

## 结果的下游分析建议

GWAS只是第一步，找到显著SNP后还需要后续分析来确认它们的生物学意义。这里提供一些思路：

首先是基因定位。根据显著SNP的染色体位置和物理坐标，在基因组注释数据库中查找附近的基因。通常认为SNP上下游100kb范围内的基因都是候选基因。可以用Ensembl Plants、Phytozome等数据库查询。

然后是功能注释。对候选基因进行GO富集分析、KEGG通路分析，看看它们的功能是否与目标性状相关。如果找到的基因在已知的代谢通路或调控网络中，可信度会大大提高。

连锁不平衡分析（LD analysis）可以帮助你确定每个关联峰值代表的独立关联位点数量。如果一个区域内多个SNP都显著，可能是因为它们相互连锁，实际上只对应一个因果变异。

表达数据整合很有价值。如果有相应的转录组数据，可以做eQTL分析看看这些SNP是否影响附近基因的表达水平，这能帮助从相关关联推进到因果关系。

最后是跨群体验证。在不同的育种群体或自然群体中重复GWAS分析，如果同一个位点被多次检测到，说明它是稳定遗传的效应位点而不是假阳性。

## 引用和致谢

如果你在研究中使用了这个脚本，请引用GAPIT的相关文献：

- GAPIT原始文献：Lipka et al. (2012) Bioinformatics, 28(18):2397-2399
- GAPIT3更新：Wang and Zhang (2021) Genomics, Proteomics & Bioinformatics, 19(4):587-591
- BLINK模型：Huang et al. (2019) GigaScience, 8(2):giy154
- FarmCPU模型：Liu et al. (2016) PLoS Genetics, 12(2):e1005767

这个脚本基于GAPIT包开发，针对HapMap格式数据的实际使用场景进行了大量优化和错误处理，力求让GWAS分析对新手更友好。如果遇到使用问题或发现bug，建议先仔细阅读错误日志找出问题原因，大多数问题都与数据格式或样本匹配有关。祝你的GWAS分析顺利！