### ============================================================
### 0. 安装并加载所需 R 包
### ============================================================
pkgs <- c("data.table", "snpReady", "openxlsx", "ComplexHeatmap", "circlize")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(pkgs, require, character.only = TRUE))

### ============================================================
### 1. 读取 HapMap hmp.txt，整理为 “样本 × SNP” 矩阵
### ============================================================
hmp_file <- "365geno.hmp.txt"

# HapMap：行 = SNP，前 11 列是注释，12 列之后是各材料基因型
hmp <- data.table::fread(hmp_file, data.table = FALSE)

# 样本名称（第 12 列以后）
sample_ids <- colnames(hmp)[-(1:11)]

# SNP map 信息：rs#, chrom, pos
map_raw <- as.matrix(hmp[, c(1, 3, 4)])
colnames(map_raw) <- c("rs", "chrom", "pos")

# Genotype 部分：SNP × 样本
geno_long <- as.matrix(hmp[, -(1:11)])

# 把常见缺失编码设为 NA
geno_long[geno_long %in% c("N", "NN", "NA", "./.")] <- NA

# 暂未处理IUPAC 编码格式（R,Y,M,S 等）


# 转成 snpReady 需要的 wide 格式：样本 × SNP
geno_wide <- t(geno_long)      # 现在：行 = 样本，列 = SNP
rownames(geno_wide) <- sample_ids
colnames(geno_wide) <- map_raw[, "rs"]  # 保证与 map 第一列一致
# 转换为矩阵并清洗非法碱基
geno_mat <- as.matrix(geno_wide) %>% toupper() %>% trimws()
allowed_bases <- c("AA", "AC", "AG", "AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
invalid_mask <- !(geno_mat %in% allowed_bases) & !is.na(geno_mat)
cat("替换", sum(invalid_mask), "个非法基因型为NA\n")
geno_mat[invalid_mask] <- NA

### ============================================================
### 2. 用 snpReady 做 SNP QC + 0/1/2 编码
### ============================================================
# 典型阈值：call.rate = 0.95, maf = 0.05，可根据项目需求调整
rd <- snpReady::raw.data(
  data         = geno_mat,
  frame        = "wide",
  hapmap       = map_raw,
  base         = TRUE,      # 基因型为碱基（ACGT / AB）
  sweep.sample = 0.80,      # 个体最大缺失率阈值（>0.80 的个体被删除）
  call.rate    = 0.95,      # SNP 最低 call rate
  maf          = 0.05,      # SNP 最低 MAF
  imput        = TRUE,      # 对剩余缺失进行填补
  imput.type   = "mean",    # 用 SNP 均值填补（稳健且快速）
  outfile      = "012",     # 输出 0/1/2 编码矩阵
  plot         = TRUE       # 额外输出一个 QC 报告 PDF
)

# QC + 编码之后的 0/1/2 矩阵：行 = 样本，列 = SNP
geno012   <- rd$M.clean
clean_ids <- rownames(geno012)   # QC 后保留的样本名
n_ind     <- nrow(geno012)
n_snp     <- ncol(geno012)

cat("QC 后样本数:", n_ind, "；SNP 数:", n_snp, "\n")

### ============================================================
### 3. 计算成对遗传距离矩阵（等位基因差异比例）
### ============================================================
# Manhattan 距离 = sum |g_i - g_j|
# 归一化到 [0,1]：d_ij = sum|g_i - g_j| / (2 * L)
dist_raw <- as.matrix(dist(geno012, method = "manhattan"))
dist_mat <- dist_raw / (2 * n_snp)

rownames(dist_mat) <- clean_ids
colnames(dist_mat) <- clean_ids

### ============================================================
### 4. 导出完整的遗传距离矩阵（Excel）
### ============================================================
# 这里直接把对称矩阵写出，行/列都是样本名
openxlsx::write.xlsx(
  as.data.frame(dist_mat),
  file     = "GeneticDistance_full_365lines.xlsx",
  rowNames = TRUE
)

### ============================================================
### 5. 为每个样本找最近 20 和最远 20 个样本，并导出 Excel
### ============================================================
k <- 20  

nearest_names  <- matrix(NA_character_, nrow = n_ind, ncol = k)
farthest_names <- matrix(NA_character_, nrow = n_ind, ncol = k)

for (i in seq_len(n_ind)) {
  this_dist <- dist_mat[i, ]
  
  # 排序（从小到大），去掉自身（距离 0）
  ord_inc <- order(this_dist, decreasing = FALSE)
  ord_inc <- ord_inc[ord_inc != i]
  
  # 从大到小（用于最远）
  ord_dec <- order(this_dist, decreasing = TRUE)
  
  near_idx <- head(ord_inc, k)
  far_idx  <- head(ord_dec, k)
  
  nearest_names[i, ]  <- clean_ids[near_idx]
  farthest_names[i, ] <- clean_ids[far_idx]
}

colnames(nearest_names)  <- paste0("Nearest_",  seq_len(k))
colnames(farthest_names) <- paste0("Farthest_", seq_len(k))

nearest_df  <- data.frame(Sample = clean_ids,
                          nearest_names,
                          check.names = FALSE)
farthest_df <- data.frame(Sample = clean_ids,
                          farthest_names,
                          check.names = FALSE)

openxlsx::write.xlsx(
  nearest_df,
  file     = "GeneticDistance_nearest20.xlsx",
  rowNames = FALSE
)

openxlsx::write.xlsx(
  farthest_df,
  file     = "GeneticDistance_farthest20.xlsx",
  rowNames = FALSE
)

### ============================================================
### 6. 绘制高质量遗传距离热图
### ============================================================
# 6.1 用距离矩阵做行/列的层次聚类
d_for_clust <- as.dist(dist_mat)
hc <- hclust(d_for_clust, method = "average")

# 6.2 为了让颜色梯度更均匀，把对角线设为 NA，只看样本间距离
dist_heat <- dist_mat
diag(dist_heat) <- NA

min_val <- min(dist_heat, na.rm = TRUE)
max_val <- max(dist_heat, na.rm = TRUE)
med_val <- median(dist_heat[!is.na(dist_heat)])

# 6.3 定义蓝-白-红渐变色标
col_fun <- circlize::colorRamp2(
  breaks = c(min_val, med_val, max_val),
  colors = c("#2166AC", "#F7F7F7", "#B2182B")
)

# 6.4 构建热图对象
ht <- ComplexHeatmap::Heatmap(
  dist_heat,
  name              = "Genetic distance",
  col               = col_fun,
  cluster_rows      = as.dendrogram(hc),
  cluster_columns   = as.dendrogram(hc),
  show_row_dend     = TRUE,
  show_column_dend  = TRUE,
  show_row_names    = FALSE,  # 样本过多时名字太密，图中一般不显示
  show_column_names = FALSE,
  row_title         = "Maize inbred lines (rows)",
  column_title      = "",
  border            = NA,
  heatmap_legend_param = list(
    title           = " ",
    title_position  = "topcenter",
    legend_direction = "vertical"
  )
)

# 6.5 导出高分辨率 PDF
pdf("GeneticDistance_heatmap.pdf", width = 8, height = 8)  # 可改成 10×10 inch
ComplexHeatmap::draw(ht)
dev.off()

# 直接预览
# ComplexHeatmap::draw(ht)
