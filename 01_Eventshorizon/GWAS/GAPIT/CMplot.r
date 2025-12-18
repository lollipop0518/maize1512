#---------------------- HEAD INFO --------------------------
project_info <- list(
  Project   = "gwas",
  Script    = "gwas_analysis.r",
  Author    = " Hanabi <anonymous401@163.com>",
  Date      = "2025-08-17",
  Purpose   = "gwas analysis for GWAS output file"
)
options(myProject.info = project_info)

# ctrl + enter Run code block
library(qqman)
library(tidyverse)
library(data.table)
library(readxl)  


gwas <- fread("result_hmp.txt") %>%          # 或 MLMstats.txt
        select(Trait, Marker, Chr, Pos, p) %>%
        mutate(Chr = as.numeric(Chr))

manhattan(gwas,
          chr = "Chr",
          bp  = "Pos",
          p   = "P.value",
          col   = myCol, 
          snp = "SNP",
          main = "Maize GWAS (MLM)",
          genomewideline = -log10(0.05/nrow(gwas)),   # Bonferroni
          suggestiveline = FALSE
)

# 曼哈顿图
manhattan(gwas,
          chr = "Chr",
          bp  = "Pos",
          p   = "p",
          snp = "Marker",
          main = "Maize GWAS (GLM)",
          genomewideline = -log10(0.05/nrow(gwas))   # Bonferroni
          #suggestiveline = -log10(1e-4))
        )
# QQ 图

qq(gwas$p, main = "Q-Q plot (GLM)")


# 提取最显著 SNP（Top 5）
top5 <- gwas %>% arrange(p) %>% slice(1:5)
print(top5)

#曼哈顿美化

library(qqman)

# 1. 先把染色体列转成 numeric 或 factor，确保 qqman 认识
gwas$Chr_num <- factor(gwas$Chr, levels = unique(gwas$Chr))
gwas$Chr_num <- as.numeric(gwas$Chr)

# 2. 自定义 20 个颜色循环（可再扩充）
myCol <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
           "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
           "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5",
           "#c49c94","#f7b6d2","#c7c7c7","#dbdb8d","#9edae5")

# 3. 画图
manhattan(gwas,
          chr   = "Chr_num",   # 用刚才的 factor
          bp    = "Pos",
          p     = "p",
          snp   = "Marker",
          col   = myCol,       # 关键：按染色体循环配色
          main  = "Maize GWAS (GLM)",
          genomewideline = -log10(0.05/nrow(gwas)),
          suggestiveline = -log10(1e-4),
          cex   = 0.6,         # 点大小
          las   = 1)           # y 轴刻度水平



# CMplot

# 若未安装：install.packages("CMplot")
library(CMplot)

# 1. 保证列名符合 CMplot 约定
#    chr、BP、P、SNP 必须对应

colnames(gwas)[match(c("Chr","Pos","p","Marker"), colnames(gwas))] <- 
  c("chr","BP","P","SNP")

# 2. 按染色体准备颜色向量（长度=染色体数即可）
chr_col <- c("#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33",
             "#a65628","#f781bf","#999999","#66c2a5","#fc8d62")

# 3. 画图


# 读取数据并转换
gwas_wide <- gwas %>%
  pivot_wider(
    names_from = Trait,   # 性状列名
    values_from = p,      # p值列名
    names_prefix = "p_"   # 新列名前缀（可选）
  ) %>%
  rename(SNP = Marker, Chromosome = Chr, Position = Pos) %>%
  select(SNP, Chromosome, Position, starts_with("p_"))

# 检查数据结构
head(gwas_wide)


CMplot(
  gwas_wide,
  plot.type = "m",         # 曼哈顿图类型
  col = chr_col,
  LOG10 = T,            # 对p值取-log10转换
  threshold = 5e-8,        # 显著性阈值（常用GWAS阈值）
  threshold.col = "red",   # 阈值线颜色
  amplify = TRUE,          # 放大显著点
  signal.col = "red",      # 显著点颜色
  ylab = "-log10(p)",      # Y轴标签
  width = 26,              # 图宽度
  height = 6               # 图高度
)

# SNP密度图

CMplot(gwas_wide,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),
        file="jpg",dpi=300,file.output=TRUE, verbose=TRUE)


# 弦曼哈顿图

CMplot(gwas_wide,plot.type="c",r=0.4,
        outward=FALSE,cir.chr.h=1.3,chr.den.col="black",file="jpg",
        dpi=300,file.output=TRUE,verbose=TRUE)


# SNP密度+曼哈顿弦图

CMplot(gwas_wide,plot.type="c",r=0.4,col=c("grey30","grey60"),
      threshold=c(1e-5,1e-4),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red",
      "blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
      bin.size=1e6,outward=FALSE,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE)


# 连锁不平衡 plot


