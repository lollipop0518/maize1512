#!/usr/bin/env Rscript

# =========================================================
# GAPIT 全流程脚本（多性状、多模型、并行、质控、汇总、作图）
# 作者：<你可以在这里写自己的信息>
# 版本：v1.0
# 依赖：R >= 4.1；GAPIT3；data.table；optparse；stringr；qqman；foreach；doParallel
# =========================================================


library(data.table)
library(stringr)
library(optparse)
library(here)

# ---------- 工具函数：缺包自动安装 ----------
install_if_missing <- function(pkgs){
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if(length(to_install)){
    message("[依赖] 将安装缺失包: ", paste(to_install, collapse=", "))
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

# 基础依赖先确保
install_if_missing(c("data.table","optparse","stringr"))

# 尝试加载 GAPIT3 或 GAPIT
gapit_pkg <- NULL
if (requireNamespace("GAPIT3", quietly = TRUE)) {
  library(GAPIT3)
  gapit_pkg <- "GAPIT3"
} else if (requireNamespace("GAPIT", quietly = TRUE)) {
  library(GAPIT)  # 某些环境包名为 GAPIT
  gapit_pkg <- "GAPIT"
} else {
  message("[错误] 未检测到 GAPIT3/GAPIT 包。请先安装：")
  message("  CRAN（若可用）：install.packages('GAPIT3')")
  message("  或者 GitHub：remotes::install_github('jiabowang/GAPIT3')")
  quit(status = 1)
}

# 作图与并行等可选依赖
optional_pkgs <- c("qqman","doParallel","foreach")
install_if_missing(optional_pkgs)

suppressPackageStartupMessages({
  library(qqman)
  library(doParallel)
  library(foreach)
})

# ---------- 命令行参数 ----------
option_list <- list(
  make_option(c("-p","--pheno"), type="character", help="表型文件路径(CSV/TSV)，需含样本列(推荐列名 Taxa)。", metavar="FILE"),
  make_option(c("--geno_hmp"), type="character", default=NULL,
              help="HapMap 基因型文件路径；多个文件用逗号分隔（支持按染色体分文件）。", metavar="FILE(S)"),
  make_option(c("--GD"), type="character", default=NULL,
              help="GD 数值型基因型矩阵文件（CSV/TSV/RDS）；行=样本，列=SNP；有行名或首列为样本名。"),
  make_option(c("--GM"), type="character", default=NULL,
              help="GM 位点注释文件（CSV/TSV/RDS）；需含 SNP, Chr/Chromosome, Pos/Position 列。"),
  make_option(c("-c","--covar"), type="character", default=NULL,
              help="协变量文件（CSV/TSV），第一列为 Taxa。"),
  make_option(c("--kinship"), type="character", default=NULL,
              help="亲缘矩阵（CSV/TSV），行列名为样本名（可为表型样本的超集）。"),
  make_option(c("-o","--outdir"), type="character", default="GAPIT_out",
              help="输出目录（默认: %default）。"),
  make_option(c("--traits"), type="character", default="all",
              help="要分析的性状列名，逗号分隔；默认 all 为全部数值列。"),
  make_option(c("--models"), type="character", default="GLM,MLM,CMLM,FarmCPU,BLINK,SUPER",
              help="GAPIT 模型列表，逗号分隔（可选：GLM,MLM,CMLM,MLMM,SUPER,FarmCPU,BLINK）。"),
  make_option(c("--n_pcs"), type="integer", default=3, help="PCA 主成分数（默认 %default）。"),
  make_option(c("--maf"), type="double", default=0.05, help="最小等位基因频率阈值 MAF（默认 %default）。"),
  make_option(c("--snp_miss"), type="double", default=0.2, help="位点允许最大缺失率（仅对GD预过滤；默认 %default）。"),
  make_option(c("--ind_miss"), type="double", default=0.2, help="样本允许最大缺失率（仅对GD预过滤；默认 %default）。"),
  make_option(c("--impute"), type="character", default="mean",
              help="GD 简单插补方法：mean(每位点均值)、major(主等位基因)；默认 %default。"),
  make_option(c("--parallel"), action="store_true", default=FALSE, help="是否按性状并行运行。"),
  make_option(c("--cores"), type="integer", default=2, help="并行核数（默认 %default）。"),
  make_option(c("--alpha_fdr"), type="double", default=0.05, help="结果汇总时的 FDR 阈值（默认 %default）。"),
  make_option(c("--plot_extra"), action="store_true", default=FALSE, help="使用 qqman 额外绘制每(性状,模型) Manhattan/QQ 图。"),
  make_option(c("--seed"), type="integer", default=1, help="随机种子（默认 %default）。")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---------- 基础检查与输出目录 ----------
if (is.null(opt$pheno)) {
  stop("必须提供 --pheno 表型文件。")
}
if (is.null(opt$geno_hmp) && (is.null(opt$GD) || is.null(opt$GM))) {
  stop("必须提供 HapMap (--geno_hmp) 或 GD/GM (--GD 与 --GM) 二者之一。")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# 日志
log_file <- file.path(opt$outdir, paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
zz <- file(log_file, open="wt"); sink(zz, split=TRUE)

message("===== GAPIT 流程启动 =====")
message("使用包：", gapit_pkg)
message("输出目录：", normalizePath(opt$outdir))

set.seed(opt$seed)

# ---------- I/O 工具 ----------
fread_any <- function(path){
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("rds","rda","rdata")) {
    obj <- readRDS(path)
    return(as.data.table(obj))
  } else {
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", tolower(path)), "\t", ",")
    return(fread(path, sep=sep, data.table = TRUE, showProgress = FALSE))
  }
}

# ---------- 读入表型 ----------
Y_raw <- fread_any(opt$pheno)
# 统一样本列名为 Taxa
taxa_col <- c("Taxa","taxa","Sample","sample","ID","id")
hit <- taxa_col[taxa_col %in% names(Y_raw)]
if (length(hit) == 0) stop("表型文件中未找到样本列（如 Taxa/Sample/ID）")
setnames(Y_raw, hit[1], "Taxa")

# 选择性状列
num_cols <- names(Y_raw)[sapply(Y_raw, is.numeric)]
trait_list <- if (tolower(opt$traits) == "all") num_cols else str_trim(unlist(strsplit(opt$traits, ",")))
trait_list <- setdiff(trait_list, "Taxa")
if (length(trait_list) == 0) stop("未识别到需要分析的数值性状列。")
message("[表型] 样本数=", nrow(Y_raw), "；候选性状数=", length(trait_list))

# ---------- 读入协变量（可选） ----------
CV <- NULL
if (!is.null(opt$covar)) {
  CV <- fread_any(opt$covar)
  hit <- taxa_col[taxa_col %in% names(CV)]
  if (length(hit) == 0) stop("协变量文件缺少样本列（如 Taxa）。")
  setnames(CV, hit[1], "Taxa")
  # 仅保留与表型有交集的样本
  CV <- CV[Taxa %in% Y_raw$Taxa]
  message("[协变量] 样本数=", nrow(CV), "；协变量列数=", ncol(CV)-1)
}

# ---------- 读入亲缘矩阵（可选） ----------
KI <- NULL
if (!is.null(opt$kinship)) {
  KI_dt <- fread_any(opt$kinship)
  # 亲缘矩阵需行为样本，列为样本；优先使用第一列为行名
  if (!is.null(KI_dt[[1]]) && !is.numeric(KI_dt[[1]])) {
    rownames_mat <- KI_dt[[1]]
    KI_mat <- as.matrix(KI_dt[,-1,with=FALSE])
    rownames(KI_mat) <- rownames_mat
    colnames(KI_mat) <- colnames(KI_dt)[-1]
  } else {
    KI_mat <- as.matrix(KI_dt)
  }
  KI <- KI_mat
  message("[亲缘矩阵] 维度：", paste(dim(KI), collapse=" x "))
}

# ---------- 读入基因型 ----------
use_hapmap <- !is.null(opt$geno_hmp)
G <- NULL; GD <- NULL; GM <- NULL

if (use_hapmap) {
  hmp_files <- str_trim(unlist(strsplit(opt$geno_hmp, ",")))
  message("[基因型] HapMap 文件数：", length(hmp_files))
  # 逐个读入并按行合并（所有文件列结构需一致）
  hap_list <- lapply(hmp_files, function(f){
    message("  - 读入 ", f)
    fread(f, data.table = TRUE, showProgress = FALSE,header = FALSE)
  })
  # 基本检查：前11列为 HapMap 固定信息列，其余为样本
  base_cols <- 11
  G <- rbindlist(hap_list, use.names=TRUE, fill=TRUE)
  message("[HapMap] SNP 数=", nrow(G), "；样本数（列）=", ncol(G)-base_cols)
} else {
  # GD/GM
  GD_dt <- fread_any(opt$GD)
  GM_dt <- fread_any(opt$GM)
  
  # 统一 GM 列名
  gm_names <- tolower(names(GM_dt))
  snp_col <- which(gm_names %in% c("snp","marker","rs","id"))
  chr_col <- which(gm_names %in% c("chr","chromosome"))
  pos_col <- which(gm_names %in% c("pos","position","bp"))
  if (length(snp_col)==0 || length(chr_col)==0 || length(pos_col)==0) {
    stop("GM 文件需包含 SNP、Chr/Chromosome、Pos/Position 列")
  }
  setnames(GM_dt, c(snp_col[1],chr_col[1],pos_col[1]), c("SNP","Chr","Pos"))
  GM <- as.data.frame(GM_dt)
  # 统一 GD 行名（样本名）
  first_col_name <- tolower(names(GD_dt)[1])
  has_taxa_col <- first_col_name %in% c("taxa","sample","id")
  if (has_taxa_col) {
    taxa_vec <- GD_dt[[1]]
    GD_mat  <- as.matrix(GD_dt[,-1,with=FALSE])
    rownames(GD_mat) <- taxa_vec
  } else if (!is.null(rownames(GD_dt))) {
    GD_mat <- as.matrix(GD_dt)
  } else {
    stop("GD 需要行名或首列为样本名（Taxa/Sample/ID）")
  }
  # 列名应为 SNP，与 GM$SNP 对应
  if (!all(colnames(GD_mat) %in% GM$SNP)) {
    warn_miss <- setdiff(colnames(GD_mat), GM$SNP)
    if (length(warn_miss) > 0) {
      message("[警告] 有 ", length(warn_miss), " 个 GD 列未在 GM$SNP 中找到，将被丢弃。")
      keep <- intersect(colnames(GD_mat), GM$SNP)
      GD_mat <- GD_mat[, keep, drop=FALSE]
    }
  }
  # 基础质控（仅 GD 路线：HapMap 可直接交给 GAPIT 在内部根据 MAF 过滤）
  message("[GD/GM] 初始维度：样本=", nrow(GD_mat), "；SNP=", ncol(GD_mat))
  # 计算样本缺失率并过滤
  ind_miss <- rowMeans(is.na(GD_mat))
  keep_ind <- ind_miss <= opt$ind_miss
  if (any(!keep_ind)) message("  - 移除样本：", sum(!keep_ind), "（缺失率超过 ", opt$ind_miss, "）")
  GD_mat <- GD_mat[keep_ind, , drop=FALSE]
  # 计算位点缺失率并过滤
  snp_miss <- colMeans(is.na(GD_mat))
  keep_snp <- snp_miss <= opt$snp_miss
  if (any(!keep_snp)) message("  - 移除位点：", sum(!keep_snp), "（缺失率超过 ", opt$snp_miss, "）")
  GD_mat <- GD_mat[, keep_snp, drop=FALSE]
  GM <- GM[GM$SNP %in% colnames(GD_mat), , drop=FALSE]
  
  # 计算 MAF 并过滤（仅 GD）
  calc_maf <- function(x){
    # x 为某 SNP 向量：0/1/2/NA
    x <- x[!is.na(x)]
    if (length(x)==0) return(NA_real_)
    p <- mean(x)/2  # 次等位基因频率（假设编码为次等位基因计数，不确定/混乱时取 min(p,1-p)）
    p <- min(p, 1-p)
    return(p)
  }
  maf_vec <- apply(GD_mat, 2, calc_maf)
  keep_maf <- (!is.na(maf_vec)) & (maf_vec >= opt$maf)
  if (any(!keep_maf)) message("  - 移除位点：", sum(!keep_maf), "（MAF <", opt$maf, "）")
  GD_mat <- GD_mat[, keep_maf, drop=FALSE]
  GM <- GM[GM$SNP %in% colnames(GD_mat), , drop=FALSE]
  
  # 简单插补
  if (opt$impute %in% c("mean","major")) {
    message("  - 对 GD 缺失进行简单插补：", opt$impute)
    if (opt$impute == "mean") {
      col_means <- colMeans(GD_mat, na.rm=TRUE)
      idx <- which(is.na(GD_mat), arr.ind = TRUE)
      if (nrow(idx) > 0) GD_mat[idx] <- col_means[idx[,2]]
    } else {
      # 主等位基因插补：将 NA 填为各列众数（取最接近 0/1/2 的频数最高值）
      fill_major <- function(v){
        tab <- table(v, useNA="no")
        if (length(tab)==0) return(v)
        maj <- as.numeric(names(tab)[which.max(tab)])
        v[is.na(v)] <- maj
        v
      }
      GD_mat <- apply(GD_mat, 2, fill_major)
    }
  }
  
  GD <- as.matrix(GD_mat)
  # 确保 GM 与 GD 列一致且顺序一致
  GM <- GM[match(colnames(GD), GM$SNP), ]
  message("[GD/GM] 过滤后维度：样本=", nrow(GD), "；SNP=", ncol(GD))
}

# ---------- 样本对齐 ----------
all_taxa <- Y_raw$Taxa
if (use_hapmap) {
  # HapMap 的样本在列（前11列为注释）
  taxa_hmp <- as.character(G[1, ])[-(1:base_cols)]  
  inter_taxa <- intersect(all_taxa, taxa_hmp)
  if (length(inter_taxa) < 10) warning("HapMap 与表型样本交集过小：", length(inter_taxa))
  Y_raw <- Y_raw[Taxa %in% inter_taxa]
  # 按表型样本顺序重排列
  ord <- match(Y_raw$Taxa, taxa_hmp)
  # 构建一个新的 G，仅保留所需样本列
  G <- cbind(G[, 1:base_cols, with=FALSE], G[, (base_cols + ord), with=FALSE])
  # 保留标题行，剔除含“/*/”的多等位行
  header  <- G[1, ]
  cleanG  <- G[-1][!grepl("/.*/", V2)]
  G_my     <- rbind(header, cleanG)
  rm(header, cleanG)
  #setnames(G, c(names(G)[1:base_cols], Y_raw$Taxa))
} else {
  # GD 行为样本
  inter_taxa <- intersect(all_taxa, rownames(GD))
  if (length(inter_taxa) < 10) warning("GD 与表型样本交集过小：", length(inter_taxa))
  Y_raw <- Y_raw[Taxa %in% inter_taxa]
  GD <- GD[inter_taxa, , drop=FALSE]
  rownames(GD) <- inter_taxa
  # 如果给了 KI，需要裁剪/重排
  if (!is.null(KI)) {
    keep <- intersect(rownames(KI), inter_taxa)
    KI <- KI[keep, keep, drop=FALSE]
    # GAPIT 会自己处理 KI 与 Y 的对齐；此处确保自洽
  }
  # 协变量对齐
  if (!is.null(CV)) {
    CV <- CV[Taxa %in% inter_taxa]
  }
}

message("[对齐] 最终样本数：", nrow(Y_raw))

# ---------- 解析模型 ----------
model_vec <- str_trim(unlist(strsplit(opt$models, ",")))
# 合法模型集合
valid_models <- c("GLM","MLM","CMLM","MLMM","SUPER","FarmCPU","BLINK")
if (!all(model_vec %in% valid_models)) {
  bad <- setdiff(model_vec, valid_models)
  stop("发现不支持的模型：", paste(bad, collapse=", "))
}
message("[模型] 将运行：", paste(model_vec, collapse=", "))

# ---------- 并行设置（按性状分发） ----------
if (isTRUE(opt$parallel)) {
  n_cores <- max(1, as.integer(opt$cores))
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  message("[并行] 开启，核数=", n_cores)
}

# ---------- 运行 GAPIT 的函数 ----------
run_gapit_one_trait <- function(trait_name){
  message(">>> 开始性状：", trait_name)
  Y2 <- data.frame(Taxa = Y_raw$Taxa, Trait = Y_raw[[trait_name]])
  # 去除缺失性状
  keep <- !is.na(Y2$Trait)
  Y2 <- Y2[keep, , drop=FALSE]
  
  # 协变量、亲缘矩阵对齐（若提供）
  CV2 <- if(!is.null(CV)) { CV[CV$Taxa %in% Y2$Taxa, , drop=FALSE] } else { NULL }
  KI2 <- if(!is.null(KI)) {
    taxa_use <- Y2$Taxa
    if (all(taxa_use %in% rownames(KI))) {
      KI[taxa_use, taxa_use, drop=FALSE]
    } else {
      NULL
    }
  } else { NULL }
  
  # 为每个性状单独建子目录以防输出覆盖
  trait_dir <- file.path()
  trait_dir <- here(opt$outdir, paste0("trait_", trait_name))
  dir.create(trait_dir, showWarnings = FALSE, recursive = TRUE)
  owd <- setwd(trait_dir); on.exit(setwd(owd), add=TRUE)
  
  # 调用 GAPIT
  # 说明：
  # - 使用 PCA.total 控制主成分个数
  # - 使用 SNP.MAF 控制最小等位基因频率（HapMap 走此过滤；GD 已预过滤）
  # - file.output=TRUE 输出全套图表与表格
  # - model 接收向量，可同时跑多个模型
  g_res <- try({
    if (use_hapmap) {
      GAPIT(
        Y = Y2,
        G = as.data.frame(G_my),
        CV = CV2,
        KI = KI2,
        PCA.total = opt$n_pcs,
        model = model_vec,
        file.output = TRUE,
        SNP.MAF = opt$maf
      )
    } else {
      GAPIT(
        Y = Y2,
        GD = GD,
        GM = GM,
        CV = CV2,
        KI = KI2,
        PCA.total = opt$n_pcs,
        model = model_vec,
        file.output = TRUE,
        SNP.MAF = opt$maf
      )
    }
  }, silent = TRUE)
  
  if (inherits(g_res, "try-error")) {
    message("[错误] GAPIT 调用失败：", attr(g_res, "condition")$message)
    return(list(trait=trait_name, ok=FALSE, result_file=NULL, combined=NULL))
  }
  
  # 收集该性状下的 GWAS 结果（各模型）
  # 一般为 *GWAS.Results.csv；不同版本命名略有差异，这里用模糊匹配
  files <- list.files(trait_dir, pattern = "GWAS.*Results|Results.*GWAS", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    # GAPIT 有时写在主目录（很少见）；再宽松搜一次
    files <- list.files(trait_dir, pattern = "^GAPIT\\.Association.*Results.*\\.csv$", recursive = TRUE, full.names = TRUE)
  }
  
  if (length(files) == 0) {
    message("[警告] 未找到 GWAS 结果表（CSV）。")
    message("trait_dir: ", trait_dir)
    message("目录下所有CSV文件：")
    print(list.files(trait_dir, pattern = "\\.csv$", recursive = TRUE))
    return(list(trait=trait_name, ok=TRUE, result_file=NULL, combined=NULL))
  }
  
  # 合并并添加模型名与性状名
  comb <- rbindlist(lapply(files, function(fp){
    dt <- fread(fp, data.table = TRUE)
    # 猜测模型名：路径或文件名中包含 GLM/MLM/...
    model_guess <- paste(valid_models[valid_models %in% toupper(strsplit(fp, "[/_.]")[[1]])], collapse="+")
    if (model_guess == "") model_guess <- NA_character_
    dt[, MODEL := model_guess]
    dt
  }), use.names = TRUE, fill = TRUE)
  
  # 标准化列名，确保至少有：SNP/Marker, Chromosome/Chr, Position/Pos, P/P.value
  nm <- tolower(names(comb))
  col_snp <- which(nm %in% c("snp","marker","rs","id","rs.#"))
  col_chr <- which(nm %in% c("chr","chromosome"))
  col_pos <- which(nm %in% c("pos","position","bp"))
  col_p   <- which(nm %in% c("p","p.value","pvalue","p.value.","p_wald","p-lrt"))
  if (length(col_snp)) setnames(comb, col_snp[1], "SNP")
  if (length(col_chr)) setnames(comb, col_chr[1], "CHR")
  if (length(col_pos)) setnames(comb, col_pos[1], "BP")
  if (length(col_p))   setnames(comb, col_p[1],   "P")
  
  comb[, TRAIT := trait_name]
  # 计算 FDR（BH）
  if (!is.null(comb$P)) {
    comb[, FDR_BH := p.adjust(P, method = "BH")]
  } else {
    comb[, FDR_BH := NA_real_]
  }
  
  # 保存该性状的合并表
  out_trait_csv <- file.path(trait_dir, paste0("GWAS_Combined_", trait_name, ".csv"))
  fwrite(comb, out_trait_csv)
  
  # 可选：生成额外 Manhattan/QQ（按模型分别绘制）
  if (isTRUE(opt$plot_extra) && all(c("CHR","BP","P") %in% names(comb))) {
    models_here <- unique(comb$MODEL)
    for(m in models_here){
      dtm <- comb[MODEL == m]
      # 为 qqman 准备列名：CHR, BP, P, SNP
      manh_png <- file.path(trait_dir, paste0("Manhattan_", trait_name, "_", m, ".png"))
      qq_png   <- file.path(trait_dir, paste0("QQ_", trait_name, "_", m, ".png"))
      try({
        png(manh_png, width=1600, height=900, res=150)
        manhattan(dtm[, .(CHR = as.numeric(CHR), BP = as.numeric(BP), P = as.numeric(P), SNP = as.character(SNP))],
                  main = paste0("Manhattan - ", trait_name, " (", m, ")"))
        dev.off()
        png(qq_png, width=1200, height=900, res=150)
        qq(as.numeric(dtm$P), main = paste0("QQ - ", trait_name, " (", m, ")"))
        dev.off()
      }, silent = TRUE)
    }
  }
  
  list(trait=trait_name, ok=TRUE, result_file=out_trait_csv, combined=comb)
}

# ---------- 主循环（可并行） ----------
runner <- if (isTRUE(opt$parallel)) {
  function(f, items) foreach(x = items, .packages = c(gapit_pkg,"data.table","qqman")) %dopar% f(x)
} else {
  function(f, items) lapply(items, f)
}

res_list <- runner(run_gapit_one_trait, trait_list)

# ---------- 全部性状的汇总 ----------
comb_all <- rbindlist(lapply(res_list, `[[`, "combined"), use.names = TRUE, fill = TRUE)
if (nrow(comb_all) > 0) {
  # 统一排序：按性状、模型、p 值
  if (!is.null(comb_all$P)) {
    setorder(comb_all, TRAIT, MODEL, P)
  }
  out_all <- file.path(opt$outdir, "GWAS_AllTraits_AllModels_Combined.csv")
  fwrite(comb_all, out_all)
  message("[汇总] 全部结果：", out_all)
  
  # 导出显著位点（FDR 阈值）
  if ("FDR_BH" %in% names(comb_all)) {
    res_dir <- "result"
    if (!dir.exists(res_dir)) {
      dir.create(res_dir, recursive = TRUE)                   # 关键！！！
      message("创建汇总输出目录：", res_dir)
    }
    sig <- comb_all[!is.na(FDR_BH) & FDR_BH <= opt$alpha_fdr]
    out_sig <- file.path(res_dir, paste0("GWAS_Significant_FDR", opt$alpha_fdr, ".csv"))
                                      # 你想放汇总结果的目录
    fwrite(sig, out_sig,row.names = FALSE, quote = FALSE)
    message("[汇总] 显著位点（FDR<=", opt$alpha_fdr, "）：", out_sig, "（行数=", nrow(sig), "）")
  }
} else {
  message("[提示] 未汇总到任何 GWAS 结果。")
}

# ---------- 清理并输出会话信息 ----------
if (isTRUE(opt$parallel)) {
  parallel::stopCluster(cl)
}

session_info_txt <- file.path(opt$outdir, "sessionInfo.txt")
writeLines(c(capture.output(sessionInfo())), con = session_info_txt)
message("[完成] 运行完毕。日志文件：", log_file)
sink(NULL); close(zz)
