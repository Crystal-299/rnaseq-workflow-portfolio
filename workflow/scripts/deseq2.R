#!/usr/bin/env Rscript
# DESeq2 差异表达分析 + PCA 质控图
# 输入：featureCounts 计数矩阵（counts.txt）、样本表（samples.tsv）
# 输出：results/deseq2/deseq2_results.csv、results/deseq2/pca.pdf
#
# 设计矩阵从样本表推导：condition 列 → design(~ condition)
# 对照方向：contrast = c(分子, 分母)，来自 config 的 deseq2.contrast
# 注意：featureCounts 的列名形如 ".../CASE_rep1_Aligned.sortedByCoord.out.bam"，
#       这里用 sub 剥掉路径与 BAM 后缀，只留样本名与 samples.tsv 对齐。

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# ---------- 极简命令行参数解析 ----------
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name) {
  i <- match(name, args)
  if (is.na(i)) stop(sprintf("缺少参数 %s", name))
  return(args[i + 1])
}

counts_file <- parse_arg("--counts")
samples_file <- parse_arg("--samples")
design_str  <- parse_arg("--design")
contrast_str <- parse_arg("--contrast")
outdir      <- parse_arg("--outdir")

contrast <- strsplit(contrast_str, ",")[[1]]
if (length(contrast) != 2) stop("--contrast 应为 'A,B' 两个元素")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------- 读计数矩阵（featureCounts 首行注释以 # 开头） ----------
cnt <- read.table(counts_file, header = TRUE, row.names = 1,
                  check.names = FALSE, comment.char = "#")
# 去掉非计数列（featureCounts 末 5 列为 Chr/Start/End/Strand/Length）
count_cols <- grepl("\\.bam$", colnames(cnt))
counts <- cnt[, count_cols, drop = FALSE]
colnames(counts) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(colnames(counts)))
counts <- round(as.matrix(counts))

# ---------- 读样本表并对齐 ----------
samples <- read.delim(samples_file, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE)
rownames(samples) <- samples$sample
samples <- samples[colnames(counts), , drop = FALSE]
samples$condition <- factor(samples$condition)
if (any(is.na(samples$condition))) {
  stop("样本表 condition 列有缺失，无法构建设计矩阵")
}
if (length(unique(samples$condition)) < 2) {
  stop("condition 种类少于 2，DESeq2 需要对照")
}

# ---------- 构建 DESeq2 对象 ----------
# design 字符串（如 "~ condition"）由 config 传入
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = as.formula(design_str))
# 过滤：至少一半样本有 ≥1 个计数，减少低表达噪声
dds <- dds[rowSums(counts(dds) >= 1) >= ceiling(ncol(counts) / 2), ]
dds <- DESeq(dds)

# ---------- 结果与输出 ----------
res <- results(dds, contrast = c("condition", contrast[1], contrast[2]))
res <- res[order(res$padj, na.last = TRUE), ]
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
                     "stat", "pvalue", "padj")]
write.csv(res_df, file.path(outdir, "deseq2_results.csv"),
          row.names = FALSE, quote = FALSE)

# ---------- PCA（vst 变换，观察重复是否聚在一起） ----------
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct <- attr(pca_data, "percentVar")
p <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  xlab(paste0("PC1: ", round(pct[1] * 100, 1), "% variance")) +
  ylab(paste0("PC2: ", round(pct[2] * 100, 1), "% variance")) +
  theme_bw()
ggsave(file.path(outdir, "pca.pdf"), p, width = 6, height = 5)

cat(sprintf("[OK] DESeq2 完成：%d 个基因，结果写入 %s\n",
            nrow(res_df), file.path(outdir, "deseq2_results.csv")))
