# Paper I (Zhang et al.): DESeq2 on HTSeq-count tables — Serum vs BHI (rich medium).
# Student Manual: count matrix then differential expression with DESeq2.

#  1. set up
#for readability
suppressPackageStartupMessages({
  library(DESeq2)
})

# directories stuff
repo <- getwd() # we're in root coming from 9 bash script
count_dir <- file.path(repo, "data/tmp/8_htseq_counts")
out_dir <- file.path(repo, "data/tmp/9_deseq")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. experiment table
samples <- data.frame(
  sample = c(
    "BH_ERR1797972", "BH_ERR1797973", "BH_ERR1797974",
    "Serum_ERR1797969", "Serum_ERR1797970", "Serum_ERR1797971"
  ),
  condition = factor(
    c(rep("BH", 3L), rep("Serum", 3L)),
    levels = c("BH", "Serum")
  ),
  stringsAsFactors = FALSE
)

# 3. read count files
samples$file <- file.path(count_dir, paste0(samples$sample, ".counts.tsv"))

missing <- samples$file[!file.exists(samples$file)]
if (length(missing) > 0L) {
  stop("Missing count files:\n", paste(missing, collapse = "\n"), call. = FALSE)
}

read_htseq_feature_counts <- function(path) {
  x <- read.table(path, sep = "\t", header = FALSE, row.names = 1L, check.names = FALSE)
  feat <- !grepl("^__", row.names(x))
  x[feat, , drop = FALSE]
}
# 4. build count matrix MAT
# apply read for all samples in samples$file (6 files)
lst <- lapply(samples$file, read_htseq_feature_counts)
# find shared gene IDs across all count files
genes <- Reduce(intersect, lapply(lst, row.names))
if (length(genes) == 0L) {
  stop("No shared gene IDs across count files.", call. = FALSE)
}
mat <- sapply(seq_along(lst), function(i) as.integer(lst[[i]][genes, 1L]))
colnames(mat) <- samples$sample
rownames(mat) <- genes
mode(mat) <- "integer"

coldata <- DataFrame(condition = samples$condition)
rownames(coldata) <- samples$sample

# 5. statistical model
dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData = coldata,
  design = ~ condition
)
dds <- DESeq(dds)

# 6. save results
res <- results(dds, contrast = c("condition", "Serum", "BH"))
res_ord <- res[order(res$padj), ]

write.csv(
  as.data.frame(res_ord),
  file = file.path(out_dir, "DESeq2_results_Serum_vs_BH.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(counts(dds, normalized = TRUE)),
  file = file.path(out_dir, "DESeq2_normalized_counts.csv"),
  row.names = TRUE
)

sink(file.path(out_dir, "DESeq2_summary.txt"))
cat("Contrast: Serum vs BH (reference = BH)\n")
cat("sizeFactors:\n")
print(sizeFactors(dds))
cat("\nGenes with padj < 0.05:", sum(!is.na(res$padj) & res$padj < 0.05, na.rm = TRUE), "\n")
cat("\n")
sessionInfo()
sink()

cat("Wrote outputs under ", out_dir, "\n", sep = "")
