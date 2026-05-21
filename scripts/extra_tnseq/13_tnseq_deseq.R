# Tn-seq DESeq2: conditional essentiality screen — HI Serum vs BHI.
# Counts from HTSeq represent transposon insertion abundance per gene.
# Genes with fewer reads in Serum (lower padj) are conditionally essential.

suppressPackageStartupMessages({
  library(DESeq2)
})

repo <- getwd()
count_dir <- file.path(repo, "data/tmp/12_tnseq_htseq")
out_dir  <- file.path(repo, "data/tmp/13_tnseq_deseq")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Tn-seq samples: HI Serum vs BHI (3 replicates each).
samples <- data.frame(
  sample = c(
    "HI_Serum_trim_ERR1801009_pass", "HI_Serum_trim_ERR1801010_pass", "HI_Serum_trim_ERR1801011_pass",
    "BHI_trim_ERR1801012_pass",      "BHI_trim_ERR1801013_pass",       "BHI_trim_ERR1801014_pass"
  ),
  condition = factor(
    c(rep("HI_Serum", 3L), rep("BHI", 3L)),
    levels = c("BHI", "HI_Serum")
  ),
  stringsAsFactors = FALSE
)

samples$file <- file.path(count_dir, paste0(samples$sample, ".counts.tsv"))
missing <- samples$file[!file.exists(samples$file)]
if (length(missing) > 0L) {
  stop("Missing count files:\n", paste(missing, collapse = "\n"), call. = FALSE)
}

read_htseq <- function(path) {
  x <- read.table(path, sep = "\t", header = FALSE, row.names = 1L, check.names = FALSE)
  feat <- !grepl("^__", row.names(x))
  x[feat, , drop = FALSE]
}

lst <- lapply(samples$file, read_htseq)
genes <- Reduce(intersect, lapply(lst, row.names))
mat <- sapply(seq_along(lst), function(i) as.integer(lst[[i]][genes, 1L]))
colnames(mat) <- samples$sample
rownames(mat) <- genes
mode(mat) <- "integer"

coldata <- DataFrame(condition = samples$condition)
rownames(coldata) <- samples$sample

dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "HI_Serum", "BHI"))
res_ord <- res[order(res$padj), ]

write.csv(as.data.frame(res_ord),
  file = file.path(out_dir, "Tnseq_DESeq2_results_HI_Serum_vs_BHI.csv"),
  row.names = TRUE)

write.csv(as.data.frame(counts(dds, normalized = TRUE)),
  file = file.path(out_dir, "Tnseq_DESeq2_normalized_counts.csv"),
  row.names = TRUE)

sink(file.path(out_dir, "Tnseq_DESeq2_summary.txt"))
cat("Contrast: HI_Serum vs BHI (reference = BHI)\n")
cat("sizeFactors:\n")
print(sizeFactors(dds))
cat("\nTn-seq genes with padj < 0.05:",
    sum(!is.na(res$padj) & res$padj < 0.05, na.rm = TRUE), "\n")
cat("padj < 0.05 & log2FC < -1 (depleted = essential in Serum):",
    sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE), "\n")
cat("padj < 0.05 & log2FC > 1 (enriched = higher fitness in Serum):",
    sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("\n")
sessionInfo()
sink()

cat("Tn-seq DESeq2 outputs under ", out_dir, "\n", sep = "")
