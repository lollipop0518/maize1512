########################################
# 日期：2025/12/14
# 内容：GBLUP (SNPReady + rrBLUP)
# 说明：
#   1. 输入：HapMap基因型 (train/pred) + TXT表型
#   2. 处理：集成用户提供的 tidyverse 合并流程与 IUPAC 转码逻辑
#   3. 模型：GBLUP (80%训练/20%验证) + 独立预测
########################################
rm(list=ls())
# 设置工作路径
setwd("D:/maize1512/01_Eventshorizon//GS/GBLUP")

# 依赖包检查与安装
install_required_packages <- function() {
  packages <- c("rrBLUP", "data.table","impute", "tidyverse", "readr","snpReady")
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      if (pkg == "impute") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install("impute")
      } else {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
}

install_required_packages()

if(!dir.exists("./result")) dir.create("./result", showWarnings=FALSE)
t1 <- proc.time()

dataset_dir <- "./dataset"

# -------------------------------------------------------------------------
# 1. 辅助函数定义
# -------------------------------------------------------------------------

# 等位基因规范化：只保留ACGT/并排序
canonical_alleles <- function(x) {
  x_clean <- gsub("[^ACGT/]", "", x)
  parts <- strsplit(x_clean, "/")
  sapply(parts, function(v) {
    v <- unique(v)
    paste(sort(v), collapse = "/")
  })
}

# IUPAC 转双碱基编码 (如 R -> AG, A -> AA)
recode_iupac_to_diploid <- function(x) {
  x <- toupper(trimws(as.character(x)))
  
  # 明确标记缺失
  x[x %in% c("", "-", "N", "NA", ".", "NN")] <- NA
  
  # 已经是两碱基且只含 A/C/G/T 的，原样保留
  is_diploid <- grepl("^[ACGT]{2}$", x)
  
  res <- x
  res[is_diploid] <- x[is_diploid]
  
  # 单字符 IUPAC 映射表
  map <- c(
    "A" = "AA", "C" = "CC", "G" = "GG", "T" = "TT",
    "R" = "AG", "Y" = "CT", "S" = "CG", "W" = "AT",
    "K" = "GT", "M" = "AC"
  )
  
  # 需要处理的：不是两碱基，又不是 NA
  single <- !is_diploid & !is.na(x)
  
  # 映射
  if (any(single)) {
    res[single] <- unname(map[x[single]])
  }
  
  # 对于不在 map 里的字母，保守处理为 NA
  unmapped <- single & is.na(res)
  if (any(unmapped)) res[unmapped] <- NA
  
  return(res)
}

# -------------------------------------------------------------------------
# 2. 读取 HapMap 并执行合并逻辑 (tidyverse)
# -------------------------------------------------------------------------

check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("缺失文件:", file_path))
  }
}

train_hmp_file <- file.path(dataset_dir, "train_data.hmp.txt")
pred_hmp_file <- file.path(dataset_dir, "pred_data.hmp.txt")

check_file_exists(train_hmp_file)
check_file_exists(pred_hmp_file)

cat("读取 HapMap 文件...\n")
# 使用 read_tsv 保持列为字符型
old_hmp <- read_tsv(train_hmp_file, col_types = cols(.default = "c"))
new_hmp <- read_tsv(pred_hmp_file, col_types = cols(.default = "c"))

snp_col_name    <- "rs#"
allele_col_name <- "alleles"
chr_col_name    <- "chrom"
pos_col_name    <- "pos"

# 检查列名是否存在
req_cols <- c(snp_col_name, allele_col_name, chr_col_name, pos_col_name)
if (!all(req_cols %in% colnames(old_hmp))) stop("建模 HapMap 列名不匹配，需包含 rs#, alleles, chrom, pos")
if (!all(req_cols %in% colnames(new_hmp))) stop("预测 HapMap 列名不匹配")

cat("提取并比对 SNP 信息...\n")

# 从旧 hmp (建模集) 抽信息
old_info <- old_hmp %>%
  transmute(
    rs_old      = .data[[snp_col_name]],
    chrom       = .data[[chr_col_name]],
    pos         = as.integer(.data[[pos_col_name]]),
    alleles_old = .data[[allele_col_name]]
  ) %>%
  mutate(al_set_old = canonical_alleles(alleles_old))

# 从新 hmp (预测集) 抽信息
new_info <- new_hmp %>%
  transmute(
    rs_new      = .data[[snp_col_name]],
    chrom       = .data[[chr_col_name]],
    pos         = as.integer(.data[[pos_col_name]]),
    alleles_new = .data[[allele_col_name]]
  ) %>%
  mutate(al_set_new = canonical_alleles(alleles_new))

# 按 chr + pos 做内连接
matched <- inner_join(old_info, new_info, by = c("chrom", "pos")) %>%
  mutate(same_alleles = al_set_old == al_set_new)

cat("SNP 等位基因一致性统计:\n")
print(table(matched$same_alleles))

# 筛选可信交集
good_snps <- matched %>% filter(same_alleles)

if (nrow(good_snps) == 0) stop("未找到位置与等位基因均一致的 SNP")

# 提取保留的行
old_keep <- old_hmp[match(good_snps$rs_old, old_hmp[[snp_col_name]]), ]
new_keep <- new_hmp[match(good_snps$rs_new, new_hmp[[snp_col_name]]), ]

# 统一新数据的 ID 与 alleles
new_keep[[snp_col_name]] <- old_keep[[snp_col_name]]
new_keep[[allele_col_name]] <- old_keep[[allele_col_name]]

# 合并数据
meta_cols <- 1:11
sample_cols_old <- (max(meta_cols) + 1):ncol(old_keep)
samples_old <- colnames(old_keep)[sample_cols_old]

sample_cols_new <- (max(meta_cols) + 1):ncol(new_keep)
samples_new <- colnames(new_keep)[sample_cols_new]

# 检查重复样本
dups <- intersect(samples_old, samples_new)
if (length(dups) > 0) cat("警告: 发现重复样本名:", paste(dups, collapse=", "), "将在合并中保留两者\n")

cat("合并 HapMap 矩阵...\n")
combined_hmp <- cbind(old_keep, new_keep[, sample_cols_new, drop = FALSE])

# -------------------------------------------------------------------------
# 3. 转码与 SNPReady 清洗
# -------------------------------------------------------------------------

cat("转码 IUPAC 至双碱基格式...\n")
geno_char <- as.matrix(combined_hmp[, -(meta_cols)]) # 行=SNP, 列=样本
# apply 按列操作 (margin=2)，这里需要对每个SNP(行)操作吗？
geno_recoded <- apply(geno_char, 2, recode_iupac_to_diploid)

# SNPReady 需要 行=样本, 列=Marker
geno_ready_input <- t(geno_recoded)
colnames(geno_ready_input) <- combined_hmp[[snp_col_name]]

# 构造 map 信息
hapmap_info <- as.matrix(combined_hmp[, c(snp_col_name, chr_col_name, pos_col_name)])
colnames(hapmap_info) <- c("rs", "chrom", "pos")

cat("运行 SNPReady QC 与填充...\n")
raw_out <- raw.data(
  data       = geno_ready_input,
  frame      = "wide",
  hapmap     = hapmap_info,
  base       = TRUE,          
  sweep.sample = 1,           # 不按缺失率删样本
  call.rate  = 0.95,          
  maf        = 0.05,          
  imput      = TRUE,
  imput.type = "mean",
  outfile    = "012",          # 不输出文件
  plot       = FALSE
)

M_clean <- raw_out$M.clean # 0/1/2 Matrix (Sample x Marker)
cat("QC 后剩余 SNP:", ncol(M_clean), "\n")

# -------------------------------------------------------------------------
# 4. 读取表型与对齐
# -------------------------------------------------------------------------

phen_file <- file.path(dataset_dir, "trait_imputed.txt")
cat("读取表型:", phen_file, "\n")
ph <- fread(phen_file, header=TRUE, data.table=FALSE, sep="\t")

trait_col_name <- "fenzhi_BLUP"
if (!trait_col_name %in% colnames(ph)) colnames(ph) <- trimws(colnames(ph))
if (!trait_col_name %in% colnames(ph)) stop(paste("未找到表型列:", trait_col_name))

# 区分样本
all_clean_ids <- rownames(M_clean)
# 建模集样本
tr_samples_clean <- intersect(samples_old, all_clean_ids)
# 预测集样本
pr_samples_clean <- intersect(samples_new, all_clean_ids)

# 建模样本需与表型取交集
ph_id_col <- colnames(ph)[1]
ph_ids <- as.character(ph[[ph_id_col]])

common_ids <- intersect(tr_samples_clean, ph_ids)
if (length(common_ids) == 0) stop("建模基因型与表型样本无交集")

cat("建模样本数 (Final):", length(common_ids), "\n")
cat("预测样本数 (Final):", length(pr_samples_clean), "\n")

# 提取矩阵
M_tr_final <- M_clean[common_ids, , drop=FALSE]
M_pr_final <- M_clean[pr_samples_clean, , drop=FALSE]

# 提取表型
y_final <- ph[match(common_ids, ph_ids), trait_col_name]

# -------------------------------------------------------------------------
# 5. GBLUP 模型构建
# -------------------------------------------------------------------------

cat("计算 A.mat...\n")
M_total <- rbind(M_tr_final, M_pr_final)
# SNPReady: 0,1,2 -> A.mat: -1,0,1
K <- A.mat(M_total - 1)

# 划分 80/20
set.seed(2023)
n_model <- length(common_ids)
n_pred <- nrow(M_pr_final)

train_idx <- sample(seq_len(n_model), size = floor(0.8 * n_model))
val_idx <- setdiff(seq_len(n_model), train_idx)

# 构建输入数据
y_vec_model <- y_final
y_vec_model[val_idx] <- NA
y_vec_all <- c(y_vec_model, rep(NA, n_pred))
gid_vec_all <- c(common_ids, rownames(M_pr_final))

data_blup <- data.frame(y = y_vec_all, gid = gid_vec_all, stringsAsFactors=FALSE)

cat("运行 kin.blup...\n")
ans <- kin.blup(data=data_blup, geno="gid", pheno="y", K=K)

# -------------------------------------------------------------------------
# 6. 结果输出
# -------------------------------------------------------------------------

# (1) 验证集
pred_val <- ans$pred[val_idx]
true_val <- y_final[val_idx]
ids_val <- common_ids[val_idx]

PCC <- suppressWarnings(cor(true_val, pred_val, use="complete.obs"))
RMSE <- sqrt(mean((true_val - pred_val)^2, na.rm=TRUE))
MAE <- mean(abs(true_val - pred_val), na.rm=TRUE)
MSE <- mean((true_val - pred_val)^2, na.rm=TRUE)

cat("验证集 PCC:", PCC, "\n")

val_res <- data.frame(SampleID=ids_val, TrueValue=true_val, Prediction=pred_val)

output_filename_val <- paste0("./result/gblup_validation_pred_", Sys.Date(), ".csv")
write.csv(val_res, output_filename_val, row.names=FALSE)

metrics_filename_val <- paste0("./result/gblup_validation_metrics_", Sys.Date(), ".csv")
write.csv(data.frame(PCC=PCC, RMSE=RMSE, MAE=MAE, MSE=MSE), metrics_filename_val, row.names=FALSE)

# (2) 预测集
pred_idx_range <- (n_model + 1):(n_model + n_pred)
pred_preds <- ans$pred[pred_idx_range]
ids_preds <- rownames(M_pr_final)

pred_res <- data.frame(SampleID=ids_preds, Prediction=pred_preds)

output_filename_pred <- paste0("./result/gblup_prediction_pred_", Sys.Date(), ".csv")
write.csv(pred_res, output_filename_pred, row.names=FALSE)

# save.image(paste0("./result/GBLUP_SNPReady_", Sys.Date(), ".RData"))
