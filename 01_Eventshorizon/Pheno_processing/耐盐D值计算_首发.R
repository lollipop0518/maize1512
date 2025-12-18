############################################################
## 0. 环境准备
############################################################
# 如未安装请先运行 install.packages("tidyverse") 等
library(tidyverse)
library(factoextra)  # 用于PCA可视化和聚类辅助
library(cluster)

############################################################
## 1. 读入数据
############################################################
# 假设你的数据为 csv，列名为：
# samples, SL_idx, RL_idx, SFW_idx, RFW_idx, SDW_idx, RDW_idx
# 每一行是一个材料/样本
dat_raw <- read_xlsx('C:/Users/夙玉/Desktop/三次预实验表型相对值.xlsx',sheet = 1)

# 指标列名称
idx_cols <- c("SL_idx","RL_idx","SFW_idx","RFW_idx","SDW_idx","RDW_idx")

# 只取数值矩阵做后续分析
dat_idx <- dat_raw[, idx_cols]
rownames(dat_idx) <- dat_raw$samples

############################################################
## 2. 主成分分析（PCA）
############################################################
# 标准化后做 PCA
pca_res <- prcomp(dat_idx, scale. = TRUE)

# 计算各主成分方差贡献率
eig_vals <- pca_res$sdev^2
var_exp  <- eig_vals / sum(eig_vals)         # 单个主成分贡献率
cum_var  <- cumsum(var_exp)                  # 累积贡献率

# 选取累积贡献率 ≥ 85% 的前 k 个主成分
k_pc <- which(cum_var >= 0.85)[1]
if (is.na(k_pc)) k_pc <- length(var_exp)     # 防止极端情况

cat("使用前", k_pc, "个主成分，累积解释率 =",
    round(cum_var[k_pc] * 100, 2), "%\n")

# 取前 k_pc 个主成分得分矩阵 (Xij)
pc_scores <- as.data.frame(pca_res$x[, 1:k_pc, drop = FALSE])

############################################################
## 3. 计算隶属函数 μ(Xij)
##   μ(Xij) = (Xij - Xjmin) / (Xjmax - Xjmin)
############################################################
mins   <- apply(pc_scores, 2, min)
maxs   <- apply(pc_scores, 2, max)
ranges <- maxs - mins
# 避免除以 0
ranges[ranges == 0] <- 1

mu_mat <- sweep(pc_scores, 2, mins, "-")
mu_mat <- sweep(mu_mat, 2, ranges, "/")
mu_mat <- as.data.frame(mu_mat)
colnames(mu_mat) <- paste0("mu_PC", 1:k_pc)

############################################################
## 4. 计算权重 Wj（各主成分的方差贡献率）
##   Wj = 第 j 个主成分的方差贡献率 / 选中主成分贡献率之和
############################################################
Wj <- var_exp[1:k_pc]
Wj <- Wj / sum(Wj)  # 归一化，权重之和=1
names(Wj) <- paste0("PC", 1:k_pc)
Wj

############################################################
## 5. 计算综合耐盐评价值 D
##   对每个材料 i：D_i = Σ_j μ(Xij) * Wj
############################################################
D_vec <- as.matrix(mu_mat) %*% Wj
D_vec <- as.numeric(D_vec)

# 合并回原始数据
dat_result <- dat_raw %>%
  mutate(D_value = D_vec)

############################################################
## 6. 基于 D 值的 K-means 聚类
############################################################
# 这里默认分成 3 类：高耐盐、中耐盐、低耐盐
# 你可以改成 2、4 等
set.seed(123)
k_cluster <- 5

km_res <- kmeans(D_vec, centers = k_cluster, nstart = 25)

# 按各簇 D 均值从小到大重新标号（1=低, 2=中, 3=高）
cluster_mean <- tapply(D_vec, km_res$cluster, mean)
order_id     <- order(cluster_mean)                     # 从小到大
new_cluster  <- match(km_res$cluster, order_id)

# 给出中文标签
label_vec <- c("低耐盐","中等耐盐","高耐盐")
label_vec <- c("高敏(HS)","敏感(S)","中耐(MT)","耐(T)","高耐(HT)")
if (length(label_vec) < k_cluster) {
  # 若簇数不是3，可简单用 C1, C2, ... 命名
  label_vec <- paste0("C", 1:k_cluster)
}
cluster_label <- factor(new_cluster,
                        levels = 1:k_cluster,
                        labels = label_vec[1:k_cluster])

dat_result$Cluster <- cluster_label

############################################################
## 7. 输出结果表格
############################################################
# 1）主成分信息表
pca_info <- data.frame(
  PC          = paste0("PC", 1:length(eig_vals)),
  Eigenvalue  = eig_vals,
  VarExpl     = round(var_exp, 4),
  CumVarExpl  = round(cum_var, 4)
)
write.csv(pca_info, "PCA_eigen_info.csv", row.names = FALSE)

# 2）样本 D 值与聚类结果
result_table <- dat_result %>%
  select(samples, all_of(idx_cols), D_value, Cluster) %>%
  arrange(desc(D_value))
write.csv(result_table, "SaltTolerance_D_and_Cluster.csv", row.names = FALSE)
write_xlsx(result_table, "D值计算结果.xlsx")
print(head(result_table))

############################################################
## 8. 可视化结果
############################################################

## 8.1 PCA 碎石图
fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, max(var_exp)*100 + 5))

## 8.2 D 值柱状图（按 D 值从高到低排序）
ggplot(result_table,
       aes(x = reorder(samples, D_value),
           y = D_value,
           fill = Cluster)) +
  geom_col() +
  coord_flip() +
  labs(x = "样本", y = "综合耐盐评价值 D",
       title = "各样本综合耐盐评价值 D 及 K-means 聚类") +
  theme_bw()

ggsave("D_value_barplot.png", width = 8, height = 6, dpi = 300)

## 8.3 D 值一维散点聚类图
ggplot(result_table,
       aes(x = D_value, y = 0, color = Cluster, label = samples)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "综合耐盐评价值 D",
       title = "基于 D 值的 K-means 聚类结果")

ggsave("D_value_cluster_scatter.png", width = 7, height = 4, dpi = 300)

############################################################
## 结束：SaltTolerance_D_and_Cluster.csv 为结果汇总表，
##       D_value_barplot.png 和 D_value_cluster_scatter.png 为主要图形
############################################################
