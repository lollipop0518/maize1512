# Load required packages
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)

# ==============================================================================
# 1. 数据读取 (请修改为你的实际文件名)
# ==============================================================================

# 读取距离矩阵
# 第一列是样本名，需要把它转换成行名，方便后续查找
raw_dist <- read_xlsx('GD_output/GeneticDistance_full_365lines.xlsx',sheet = 1)
dist_mat <- as.data.frame(raw_dist)

# 规定第一列是样本ID，将其设为行名，并删除第一列
rownames(dist_mat) <- dist_mat[[1]]
dist_mat <- dist_mat[, -1] 

# 读取最近最远邻列表 (365行, Sample + 20个est)
nearest_df <- read_excel("GD_output/GeneticDistance_nearest20.xlsx")
farthest_df <- read_excel("GD_output/GeneticDistance_farthest20.xlsx")

# ==============================================================================
# 2. 核心处理逻辑
# ==============================================================================

# 初始化一个空的列表来存储结果，比循环rbind速度快
output_list <- list()

far_or_near <- farthest_df # make u choice

# 遍历每一个样本 (每一行)
for(i in 1:nrow(far_or_near)) {
  
  # --- A. 获取数据 ---
  # 当前核心样本 (Focal Sample)
  focal_id <- as.character(far_or_near[i, 1]) 
  
  # 当前样本的20个邻居 ID (第2列到第21列)
  neighbor_ids <- as.character(far_or_near[i, 2:ncol(far_or_near)])
  
  # --- B. 构建第一行：名字行 (Name Row) ---
  # 格式: [核心ID, 邻居1_ID, 邻居2_ID, ...]
  row_names <- c(focal_id, neighbor_ids)
  
  # --- C. 构建第二行：数值行 (Value Row) ---
  # 从距离矩阵中查找 focal_id 和每一个 neighbor_id 之间的距离
  
  # 使用 sapply 批量查找
  distances <- sapply(neighbor_ids, function(nb_id) {
    # 错误处理：防止矩阵中找不到该ID报错
    if(focal_id %in% rownames(dist_mat) && nb_id %in% colnames(dist_mat)) {
      # 查找矩阵对应行列的值
      val <- dist_mat[focal_id, nb_id]
      # 格式化，保留3位或4位小数
      return(sprintf("%.4f", val)) 
    } else {
      return(NA) # 找不到则填NA
    }
  })
  
  # 格式: [空字符串, 距离1, 距离2, ...]
  row_values <- c("", distances)
  
  # --- D. 存入列表 ---
  output_list[[length(output_list) + 1]] <- row_names
  output_list[[length(output_list) + 1]] <- row_values
}

# ==============================================================================
# 3. 结果整理与格式化
# ==============================================================================

# 将列表转换为数据框
final_df <- do.call(rbind, output_list)
final_df <- as.data.frame(final_df, stringsAsFactors = FALSE)

# 设置列名 (参考图片: ID, 材料.1, 材料.2 ...)
col_headers <- c("ID", paste0("材料.", 1:20))
colnames(final_df) <- col_headers

# 添加序号列
seq_column <- 1:nrow(final_df)
final_df <- cbind(序号 = seq_column, final_df)

# ==============================================================================
# 4. 导出结果
# ==============================================================================

# 写入 Excel
write_xlsx(final_df, "Formatted_Genotype_Distance.xlsx")

print("done！file saved with Formatted_Genotype_Distance.xlsx")

# preview
head(final_df)

