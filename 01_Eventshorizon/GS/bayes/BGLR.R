# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------蒙特卡洛交叉验证（Monte Carlo Cross-Validation）---------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

library(BGLR)
library(dplyr)
library(snpReady)
# ========== 1. 数据读取与预处理 ==========
# 表型数据读取
pheno_file <- 'D:/Bioinfo/LH_gapit/lh_gapit/01_cleaned_data/Height_BLUP_SpATS.txt' 
pheno_df <- read.csv(pheno_file, stringsAsFactors = FALSE, sep = '\t')
colnames(pheno_df)[1:2] <- c("Taxa", "TraitValue")  # 统一列名
pheno_df <- pheno_df %>% filter(!is.na(TraitValue))  # 过滤表型缺失

# 基因型数据读取
geno_raw <- read.csv(
  "D:/Bioinfo/LH_gapit/lh_gapit/00_original_data/genome_data/taxa-prder.hmp.txt",
  stringsAsFactors = FALSE, check.names = FALSE, sep = '\t'
)

sample_ids <- colnames(geno_raw)[12:ncol(geno_raw)]  # 提取样本ID
snp_ids <- geno_raw[, 1]                             # 提取SNP ID
geno_calls <- geno_raw[, 12:ncol(geno_raw)]          # 提取纯基因型矩阵
rownames(geno_calls) <- snp_ids
geno_t <- t(geno_calls)                              # 转置：行=样本，列=SNP

# 样本匹配（表型+基因型交集）
common_taxa <- intersect(pheno_df$Taxa, rownames(geno_t))
cat("共匹配到", length(common_taxa), "个样本（表型+基因型）\n")
if (length(common_taxa) == 0) stop("无匹配样本！检查Taxa列命名")

# 对齐表型和基因型数据
pheno_final <- pheno_df %>% filter(Taxa %in% common_taxa) %>% arrange(Taxa)
geno_t_final <- geno_t[common_taxa, ]
if (!all(pheno_final$Taxa == rownames(geno_t_final))) stop("数据对齐失败")

# ========== 2. 基因型清洗（非法字符+质控） ==========
# 转换为矩阵并清洗非法碱基
geno_mat <- as.matrix(geno_t_final) %>% toupper() %>% trimws()
allowed_bases <- c("AA", "AC", "AG", "AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
invalid_mask <- !(geno_mat %in% allowed_bases) & !is.na(geno_mat)
cat("替换", sum(invalid_mask), "个非法基因型为NA\n")
geno_mat[invalid_mask] <- NA

# SNP质控（过滤+填补）
geno_processed <- raw.data(
  data = geno_mat,
  frame = "wide",        # 数据格式：宽矩阵
  base = TRUE,           # 输入为碱基格式
  sweep.sample = 0.5,    # 过滤缺失率>50%的样本
  call.rate = 0.90,      # 过滤检出率<90%的SNP
  maf = 0.05,            # 过滤MAF<0.05的SNP
  imput = TRUE,          # 填补缺失值
  imput.type = "mean",   # 均值填补
  outfile = FALSE        # 不输出中间文件
)

# 提取质控后基因型矩阵（-1/0/1编码）
M <- geno_processed$M.clean
cat("基因型质控完成：", dim(M)[1], "样本 ×", dim(M)[2], "SNP\n")
cat("过滤掉", ncol(geno_t_final) - ncol(M), "个低质量SNP\n\n")

# ========== 3. 分析数据准备 ==========
# 提取表型和基因型矩阵
y <- pheno_final$TraitValue
X <- as.matrix(M)

# 过滤表型缺失样本（保持表型-基因型同步）
complete_idx <- !is.na(y)
y <- y[complete_idx]
X <- X[complete_idx, ]

# 基因型矩阵标准化（BGLR推荐）
X <- scale(X, center = TRUE, scale = TRUE)


# === 确保样本名存在 用于后续保存每次抽样结果 ===
sample_ids <- pheno_final$Taxa

# === 预处理：标准化 (可做) ===
y_raw <- Y[,1]
y <- as.vector(scale(y_raw)) 
n <- length(y)

# 2. 设置保存路径
output_dir <- "BGLR_Results"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 3. 设置重复次数
n_repeats <- 50
correlations <- numeric(n_repeats)

cat("开始运行循环，结果将保存在:", output_dir, "\n")

# 4. 开始循环
for(i in 1:n_repeats) {
  
  # --- A. 随机划分 ---
  tst <- sample(1:n, size = round(0.2 * n)) # 20% 测试集
  
  y_NA <- y
  y_NA[tst] <- NA # 遮挡测试集真实值
  
  # --- B. 运行 BGLR---
  skip <- FALSE
  tryCatch({
    fm <- BGLR(y = y_NA, 
               ETA = list(list(X = X, model = "BRR")), 
               nIter = 2000, burnIn = 500, verbose = FALSE)
  }, error = function(e) {
    message(paste("Round", i, "Error:", e$message))
    skip <<- TRUE
  })
  
  if(skip) { correlations[i] <- NA; next }
  
  # --- C. 计算准确度 ---
  pred_val <- fm$yHat[tst]
  true_val <- y[tst] 
  correlations[i] <- cor(pred_val, true_val)
  
  # --- D. 整理并保存详细结果到文件---
  
  # 构建数据框：包含 样本名、真实值、预测值、分组信息
  result_df <- data.frame(
    SampleID = sample_ids,
    Observed = y,        # 真实值 (标准化后的)
    Predicted = fm$yHat, # 预测值/育种值
    Group = "Training"   # 默认为训练集
  )
  
  # 将测试集的 Group 标记修改为 "Testing"
  result_df$Group[tst] <- "Testing"
  
  # 生成文件名
  file_name <- paste0(output_dir, "/Round_", i, ".csv")
  
  # 写入 CSV 文件
  write.csv(result_df, file_name, row.names = FALSE, quote = FALSE)
  
  # 打印进度
  cat("Round", i, "Done. Accuracy:", round(correlations[i], 3), "-> Saved to", file_name, "\n")
}

# --- 1. 保存准确度汇总表格 ---
# 构建汇总数据框
summary_df <- data.frame(
  Round = 1:n_repeats,
  Accuracy = correlations,
  Status = ifelse(is.na(correlations), "Failed", "Success") # 标记是否有失败的轮次
)

# 保存为 CSV 文件
summary_file <- paste0(output_dir, "/All_Rounds_Accuracy.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("准确度汇总表至:", summary_file, "\n")

# --- 2. 可视化并保存图片 ---
# 设置图片保存路径 (PNG格式)
plot_file <- paste0(output_dir, "/Accuracy_Distribution.png")
png(filename = plot_file, width = 800, height = 600, res = 120)

# 绘图布局参数
par(mar = c(5, 5, 4, 2))

# A. 绘制箱线图 (Boxplot) 展示整体分布
boxplot(correlations, 
        main = paste("Prediction Accuracy over", n_repeats, "Repeats"),
        ylab = "Pearson Correlation (r)",
        col = "lightblue",
        border = "darkblue",
        outline = FALSE, # 不显示离群点，下面用散点代替
        ylim = c(min(correlations, na.rm=TRUE) * 0.9, max(correlations, na.rm=TRUE) * 1.1))

# B. 添加抖动散点 (Jitter Points) - 让每个点都显示出来
stripchart(correlations, 
           method = "jitter", 
           vertical = TRUE, 
           pch = 19, 
           col = rgb(0, 0, 0, 0.5), # 半透明黑
           add = TRUE)

# C. 添加平均值红线
mean_acc <- mean(correlations, na.rm = TRUE)
abline(h = mean_acc, col = "red", lwd = 2, lty = 2)
text(x = 1.3, y = mean_acc, labels = paste("Mean:", round(mean_acc, 3)), col = "red", pos = 3)

# 关闭绘图设备（保存文件）
dev.off()
cat("已保存可视化统计图至:", plot_file, "\n")
cat("--------------------------------\n")
