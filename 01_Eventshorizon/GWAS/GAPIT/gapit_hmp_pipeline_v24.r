#!/usr/bin/env Rscript

# =========================================================
# GAPIT HapMap 全流程脚本（精简稳定版 v2.4）
# 仅支持 HapMap (G) 输入；删除 GD/GM 支持；删除并行相关代码
#
# 关键修复：
# 1) HapMap 必须用 header=FALSE 方式读入：让“标题行”保留在第1行（这是 GAPIT 推荐做法）
# 2) HapMap 样本名从第1行、第12列开始提取并用于对齐（不是从列名提取）
#
# 依赖：R >= 4.1；GAPIT 包；data.table；optparse
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# ---------- 工具函数 ----------
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
}

safe_stop <- function(msg) {
  log_message(msg, "ERROR")
  quit(status = 1)
}

# 自动安装缺失的基础依赖（可按需注释掉）
install_if_missing <- function(pkgs){
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if(length(miss)){
    message("[依赖] 将安装缺失包: ", paste(miss, collapse = ", "))
    install.packages(miss, repos = "https://cloud.r-project.org")
  }
}
install_if_missing(c("data.table","optparse"))

# ---------- 命令行参数 ----------
option_list <- list(
  make_option(c("-p","--pheno"), type="character", help="表型文件路径(CSV/TSV/TXT)，必须包含样本列 Taxa/taxa/Sample/ID", metavar="FILE"),
  make_option(c("-g","--geno_hmp"), type="character", help="HapMap 基因型文件路径；多文件用逗号分隔（文件本身有标题行，但我们会用 header=FALSE 读入）", metavar="FILE(S)"),
  make_option(c("-o","--outdir"), type="character", default="GAPIT_out", help="输出目录（默认: %default）"),
  make_option(c("--traits"), type="character", default="all", help="分析性状列名，逗号分隔；all=自动识别数值列"),
  make_option(c("--models"), type="character", default="BLINK", help="GAPIT模型列表，逗号分隔（如 GLM,MLM,FarmCPU,BLINK）"),
  make_option(c("--n_pcs"), type="integer", default=3, help="PCA主成分数（默认 %default）"),
  make_option(c("--maf"), type="double", default=0.05, help="最小等位基因频率 MAF（默认 %default）"),
  make_option(c("--seed"), type="integer", default=1, help="随机种子（默认 %default）"),
  make_option(c("--hmp_header_cols"), type="integer", default=11, help="HapMap 前导列数（标准为11；默认 %default）")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---------- 参数验证 ----------
if (is.null(opt$pheno) || !nzchar(opt$pheno)) safe_stop("必须提供 --pheno 表型文件")
if (!file.exists(opt$pheno)) safe_stop(paste0("表型文件不存在: ", opt$pheno))

if (is.null(opt$geno_hmp) || !nzchar(opt$geno_hmp)) safe_stop("必须提供 --geno_hmp HapMap 文件")
hmp_files <- trimws(unlist(strsplit(opt$geno_hmp, ",")))
hmp_files <- hmp_files[nzchar(hmp_files)]
for (f in hmp_files) if (!file.exists(f)) safe_stop(paste0("HapMap文件不存在: ", f))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- 日志设置 ----------
log_file <- file.path(opt$outdir, paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
log_conn <- file(log_file, open="wt")
sink(log_conn, split=TRUE)
on.exit({ try(sink(), silent=TRUE); try(close(log_conn), silent=TRUE) }, add=TRUE)

log_message("===== GAPIT HapMap 流程启动 (精简稳定版 v2.4) =====")
log_message(paste("输出目录:", normalizePath(opt$outdir, winslash="/", mustWork = FALSE)))
set.seed(opt$seed)

# ---------- 加载 GAPIT ----------
if (requireNamespace("GAPIT", quietly = TRUE)) {
  suppressPackageStartupMessages(library(GAPIT))
  log_message("使用包: GAPIT")
} else if (requireNamespace("GAPIT3", quietly = TRUE)) {
  suppressPackageStartupMessages(library(GAPIT3))
  log_message("使用包: GAPIT3")
} else {
  safe_stop("未检测到 GAPIT/GAPIT3 包。请先安装：remotes::install_github('jiabowang/GAPIT') 或 install.packages('GAPIT')（若有）")
}

# ---------- 读入表型 ----------
fread_any <- function(path){
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("txt","tsv")) "\t" else ","
  # check.names=FALSE 防止列名被改写（表型一般无所谓，但稳妥）
  as.data.frame(fread(path, sep=sep, header=TRUE, data.table=FALSE, check.names=FALSE, showProgress=FALSE))
}

Y_raw <- fread_any(opt$pheno)

# 统一样本列名
taxa_col_options <- c("Taxa","taxa","Sample","sample","ID","id")
hit <- intersect(taxa_col_options, names(Y_raw))
if (length(hit) == 0) safe_stop("表型文件缺少样本列（需要: Taxa/taxa/Sample/ID）")
names(Y_raw)[names(Y_raw) == hit[1]] <- "Taxa"

# 样本名清理
Y_raw$Taxa <- trimws(as.character(Y_raw$Taxa))

# 去重检查
if (any(duplicated(Y_raw$Taxa))) {
  dupn <- sum(duplicated(Y_raw$Taxa))
  log_message(sprintf("警告: 表型样本名存在重复（%d个重复行），将保留第一次出现", dupn), "WARN")
  Y_raw <- Y_raw[!duplicated(Y_raw$Taxa), , drop=FALSE]
}

# 选择性状列
num_cols <- names(Y_raw)[sapply(Y_raw, is.numeric)]
trait_list <- if (tolower(opt$traits) == "all") {
  num_cols
} else {
  trimws(unlist(strsplit(opt$traits, ",")))
}
trait_list <- setdiff(trait_list, "Taxa")
trait_list <- trait_list[trait_list %in% names(Y_raw)]

if (length(trait_list) == 0) safe_stop("未识别到需要分析的性状列（--traits 或自动识别均为空）")
log_message(sprintf("表型: 样本数=%d, 性状数=%d", nrow(Y_raw), length(trait_list)))

# ---------- 读入 HapMap（关键：header=FALSE） ----------
read_hmp_one <- function(path){
  # 说明：
  # GAPIT 官方教程与维护者回复均建议：HapMap 用 head=FALSE 读入，使“标题行”作为第1行数据存在。
  # 这样样本ID在第1行、第12列以后，不会受到 R 对列名 make.names() 的影响。
  fread(path, header=FALSE, data.table=TRUE, check.names=FALSE, showProgress=FALSE)
}

log_message(sprintf("HapMap文件数: %d", length(hmp_files)))
hmp_list <- vector("list", length(hmp_files))
for (i in seq_along(hmp_files)) {
  f <- hmp_files[i]
  log_message(paste0("  读入: ", basename(f)))
  dt <- read_hmp_one(f)
  if (nrow(dt) < 2) safe_stop(paste0("HapMap文件行数异常(不足2行): ", f))
  if (i > 1) {
    # 多文件合并时去掉重复标题行
    dt <- dt[-1]
  }
  hmp_list[[i]] <- dt
}

# 列数一致性检查
ncols <- sapply(hmp_list, ncol)
if (length(unique(ncols)) != 1) {
  safe_stop(sprintf("多个HapMap文件列数不一致: %s；请确保样本集合与列数完全一致后再合并",
                    paste(ncols, collapse=", ")))
}

G <- rbindlist(hmp_list, use.names = FALSE, fill = FALSE)
rm(hmp_list); gc(verbose = FALSE)

base_cols <- as.integer(opt$hmp_header_cols)
if (ncol(G) <= base_cols) {
  safe_stop(sprintf("HapMap列数(%d) <= 前导列数(%d)，未检测到样本列；请检查 --hmp_header_cols 是否正确（标准=11）",
                    ncol(G), base_cols))
}

# ---------- HapMap 样本名从“第一行”提取（关键修复） ----------
taxa_hmp <- as.character(unlist(G[1, (base_cols+1):ncol(G), with=FALSE], use.names=FALSE))
taxa_hmp <- trimws(taxa_hmp)

log_message(sprintf("HapMap前导列数: %d", base_cols))
log_message(sprintf("检测到HapMap样本数: %d", length(taxa_hmp)))
log_message(sprintf("前5个HapMap样本名: %s", paste(head(taxa_hmp, 5), collapse=", ")))
log_message(sprintf("前5个表型样本名: %s", paste(head(Y_raw$Taxa, 5), collapse=", ")))

# 样本交集
inter_taxa <- intersect(Y_raw$Taxa, taxa_hmp)
log_message(sprintf("样本交集数: %d (表型=%d, HapMap=%d)",
                    length(inter_taxa), nrow(Y_raw), length(taxa_hmp)))

if (length(inter_taxa) < 10) {
  log_message("样本名匹配失败诊断:", "ERROR")
  log_message(sprintf("  - 表型样本示例: %s", paste(head(Y_raw$Taxa, 5), collapse=", ")), "ERROR")
  log_message(sprintf("  - HapMap样本示例: %s", paste(head(taxa_hmp, 5), collapse=", ")), "ERROR")
  log_message("  - 关键检查：HapMap 是否按 head=FALSE 方式读入/传入 GAPIT；样本名是否有空格/前缀/后缀差异", "ERROR")
  safe_stop("HapMap与表型样本交集过小，无法继续分析")
}

# 过滤表型到交集
Y_raw <- Y_raw[Y_raw$Taxa %in% inter_taxa, , drop=FALSE]

# 预先裁剪 HapMap 到交集样本（降低内存）
ord_all <- match(Y_raw$Taxa, taxa_hmp)
keep_cols <- c(1:base_cols, base_cols + ord_all)
G <- G[, keep_cols, with=FALSE]

log_message(sprintf("样本对齐完成: 最终样本数=%d", nrow(Y_raw)))

# ---------- 解析模型 ----------
model_vec <- trimws(unlist(strsplit(opt$models, ",")))
model_vec <- model_vec[nzchar(model_vec)]
if (length(model_vec) == 0) safe_stop("未指定模型 --models")

log_message(sprintf("运行模型: %s", paste(model_vec, collapse=", ")))

# ---------- 单性状运行函数 ----------
run_gapit_one_trait <- function(trait_name, Y_all, G_all, base_cols, model_list, n_pc, maf_th, outdir_base){

  log_message(sprintf(">>> 开始分析性状: %s", trait_name))

  # 组装Y（GAPIT 要求：第一列是样本ID）
  Y2 <- data.frame(
    Taxa  = Y_all$Taxa,
    Trait = Y_all[[trait_name]],
    stringsAsFactors = FALSE
  )
  Y2 <- Y2[!is.na(Y2$Trait), , drop=FALSE]

  if (nrow(Y2) < 10) {
    log_message(sprintf("跳过 %s: 有效样本不足10个", trait_name), "WARN")
    return(list(trait=trait_name, ok=FALSE, reason="insufficient_samples"))
  }

  # 从（当前）HapMap 第一行提取样本名并按Y2过滤/重排列
  taxa_hmp2 <- as.character(unlist(G_all[1, (base_cols+1):ncol(G_all), with=FALSE], use.names=FALSE))
  taxa_hmp2 <- trimws(taxa_hmp2)

  idx <- match(Y2$Taxa, taxa_hmp2)
  if (any(is.na(idx))) {
    miss <- Y2$Taxa[is.na(idx)]
    log_message(sprintf("部分样本不在基因型中 [%s]: %s", trait_name, paste(head(miss, 5), collapse=", ")), "ERROR")
    return(list(trait=trait_name, ok=FALSE, reason="genotype_mismatch"))
  }

  keep_cols2 <- c(1:base_cols, base_cols + idx)
  G2 <- G_all[, keep_cols2, with=FALSE]

  log_message(sprintf("[%s] Y2样本数: %d, G2样本数: %d",
                     trait_name, nrow(Y2), ncol(G2) - base_cols))

  # 性状输出目录
  trait_dir <- file.path(opt$outdir, sprintf("trait_%s", trait_name))
  # trait_dir <- file.path(outdir_base, sprintf("trait_%s", trait_name))

  dir.create(trait_dir, showWarnings = FALSE, recursive = TRUE)
  setwd(trait_dir)
  # 调用 GAPIT
  g_res <- tryCatch({
    GAPIT(
      Y = Y2,
      G = as.data.frame(G2),   
      PCA.total = n_pc,
      model = model_list,
      SNP.MAF = maf_th,
      file.output = TRUE
    )
  }, error = function(e) {
    log_message(sprintf("GAPIT运行错误 [%s]: %s", trait_name, e$message), "ERROR")
    return(NULL)
  })

  if (is.null(g_res)) {
    return(list(trait=trait_name, ok=FALSE, reason="gapit_error"))
  }

  # 收集结果文件
  setwd(trait_dir)
  result_files <- list.files(trait_dir,
                             pattern = "GAPIT\\.Association\\.GWAS_Results\\.[A-Za-z0-9_]+\\.Trait\\.csv$",
                             full.names = TRUE,
							 recursive = TRUE,
                             ignore.case = TRUE)

  # 调试输出，查看文件列表
  cat("查看文件列表: \n")
  print(result_files)
  if (length(result_files) == 0) {	
	log_message(sprintf("GWAS结果文件保存于 [%s] ", trait_dir), "WARN")
    return(list(trait=trait_name, ok=FALSE, reason="no_results_file"))
  }

  combined <- rbindlist(lapply(result_files, function(fp){
    dt <- fread(fp, data.table = TRUE, showProgress = FALSE)
    dt[, SOURCE_FILE := basename(fp)]
    dt
  }), use.names = TRUE, fill = TRUE)

  # 标准化常用列名（尽量兼容不同GAPIT版本输出）
  col_map <- list(
    SNP = c("snp","marker","rs","id","rs.#"),
    CHR = c("chr","chromosome"),
    BP  = c("pos","position","bp"),
    P   = c("p","p.value","pvalue","p_wald")
  )
  for (target in names(col_map)) {
    candidates <- col_map[[target]]
    hit_idx <- which(tolower(names(combined)) %in% candidates)
    if (length(hit_idx) > 0) setnames(combined, hit_idx[1], target)
  }

  combined[, TRAIT := trait_name]
  if ("P" %in% names(combined)) combined[, FDR_BH := p.adjust(P, method="BH")]

  out_csv <- file.path(trait_dir, sprintf("GWAS_Combined_%s.csv", trait_name))
  fwrite(combined, out_csv)

  log_message(sprintf("<<< 完成性状: %s (结果数: %d)", trait_name, nrow(combined)))
  return(list(trait=trait_name, ok=TRUE, result_file=out_csv, combined=combined))
}

# ---------- 主循环 ----------
res_list <- lapply(trait_list, run_gapit_one_trait,
                   Y_all = Y_raw,
                   G_all = G,
                   base_cols = base_cols,
                   model_list = model_vec,
                   n_pc = opt$n_pcs,
                   maf_th = opt$maf,
                   outdir_base = opt$outdir)

# ---------- 汇总 ----------
ok_list <- Filter(function(x) isTRUE(x$ok) && !is.null(x$combined), res_list)

if (length(ok_list) == 0) {
  log_message("无任何性状成功完成分析", "ERROR")
} else {
  log_message(sprintf("成功完成 %d/%d 个性状", length(ok_list), length(trait_list)))

  all_dt <- rbindlist(lapply(ok_list, `[[`, "combined"), use.names = TRUE, fill = TRUE)
  if (nrow(all_dt) > 0) {
    out_all <- file.path(opt$outdir, "GWAS_AllTraits_AllModels_Combined.csv")
    fwrite(all_dt, out_all)
    log_message(sprintf("汇总结果已保存: %s", out_all))
  }
}

log_message("===== 流程完成 =====")
log_message(sprintf("日志文件: %s", log_file))
