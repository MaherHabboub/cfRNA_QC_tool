# 01_normalize_log.R
# DESeq2 normalization + log2(normalized + 1)
# Input: HTSeq count matrix with gene_id as first column and samples as columns
# Output:
#   - counts_deseq2_normalized.tsv
#   - counts_deseq2_log2norm_plus1.tsv
#   - QC plots before/after normalization

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(data.table)
  library(ggplot2)
})

# ---- command-line arguments ----
option_list <- list(
  make_option(
    c("--counts"),
    type = "character",
    default = "./data/silverseq_htseq_counts_combined.tsv",
    help = "Input HTSeq count matrix TSV. First column should be gene IDs, remaining columns samples."
  ),
  make_option(
    c("--out"),
    type = "character",
    default = "./results/norm_log",
    help = "Output directory for normalized/log counts and QC plots."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

counts_path <- opt$counts
out_dir <- opt$out

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(counts_path)) {
  stop("Counts file not found: ", counts_path)
}

# ---- load counts ----
dt <- fread(counts_path, sep = "\t", header = TRUE)

if (ncol(dt) < 2) {
  stop("Counts matrix must have at least 2 columns: gene_id + at least one sample.")
}

# Rename first column to gene_id for consistency
setnames(dt, 1, "gene_id")

gene_col <- "gene_id"
sample_cols <- setdiff(names(dt), gene_col)

if (length(sample_cols) < 1) {
  stop("No sample columns found after gene_id column.")
}

# Drop HTSeq summary rows if present
dt <- dt[!grepl("^__", dt[[gene_col]]), ]

# Remove empty/missing gene IDs
dt <- dt[!is.na(dt[[gene_col]]) & dt[[gene_col]] != "", ]

# Make sure sample columns are numeric
for (sc in sample_cols) {
  dt[[sc]] <- suppressWarnings(as.numeric(dt[[sc]]))
}

if (anyNA(dt[, ..sample_cols])) {
  stop("NA values detected after converting count columns to numeric. Check input file formatting.")
}

# ---- build count matrix ----
count_mat <- as.matrix(dt[, ..sample_cols])
rownames(count_mat) <- dt[[gene_col]]
colnames(count_mat) <- sample_cols

# DESeq2 expects integer counts
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

if (any(count_mat < 0, na.rm = TRUE)) {
  stop("Negative counts detected. DESeq2 requires non-negative integer counts.")
}

# ---- DESeq2 size-factor normalization ----
coldata <- data.frame(dummy = rep(1, ncol(count_mat)))
rownames(coldata) <- colnames(count_mat)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ 1
)

dds <- dds[rowSums(counts(dds)) > 0, ]

dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)
log_norm <- log2(norm_counts + 1)

# ---- save outputs ----
fwrite(
  as.data.table(norm_counts, keep.rownames = "gene_id"),
  file.path(out_dir, "counts_deseq2_normalized.tsv"),
  sep = "\t"
)

fwrite(
  as.data.table(log_norm, keep.rownames = "gene_id"),
  file.path(out_dir, "counts_deseq2_log2norm_plus1.tsv"),
  sep = "\t"
)

# ---- plot 1: before/after histogram for first sample ----
samp <- colnames(count_mat)[1]

df_before <- data.frame(
  x = log2(count_mat[, samp] + 1),
  which = "Before: log2(raw + 1)"
)

df_after <- data.frame(
  x = log_norm[, samp],
  which = "After: log2(DESeq2-normalized + 1)"
)

df <- rbind(df_before, df_after)

df$which <- factor(
  df$which,
  levels = c(
    "Before: log2(raw + 1)",
    "After: log2(DESeq2-normalized + 1)"
  )
)

p_hist <- ggplot(df, aes(x = x)) +
  geom_histogram(bins = 60) +
  facet_wrap(~which, ncol = 1, scales = "free_y") +
  theme_minimal() +
  labs(
    title = paste("Distribution check:", samp),
    x = "Expression value",
    y = "Number of genes"
  )

ggsave(
  file.path(out_dir, "qc_before_after_hist.png"),
  p_hist,
  width = 7,
  height = 6
)

# ---- plot 2: library size raw vs normalized ----
lib_raw <- colSums(count_mat)
lib_norm <- colSums(norm_counts)

lib_df <- data.frame(
  sample = names(lib_raw),
  raw = as.numeric(lib_raw),
  normalized = as.numeric(lib_norm)
)

p_lib <- ggplot(lib_df, aes(x = raw, y = normalized, label = sample)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Library size check: raw vs normalized",
    x = "Raw library size",
    y = "Normalized library size (sum of normalized counts)"
  )

ggsave(
  file.path(out_dir, "qc_librarysize_raw_vs_norm.png"),
  p_lib,
  width = 6,
  height = 5
)

# ---- plot 3A: density before normalization, nonzero genes only ----
raw_log <- log2(count_mat + 1)

raw_dt <- as.data.table(raw_log, keep.rownames = "gene_id")
raw_long <- melt(
  raw_dt,
  id.vars = "gene_id",
  variable.name = "sample",
  value.name = "value"
)

# Remove zero-count genes for readability
raw_long_nz <- raw_long[value > 0]

p_density_before <- ggplot(raw_long_nz, aes(x = value, color = sample)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Density curves before normalization",
    subtitle = "Nonzero genes only; values are log2(raw counts + 1)",
    x = "log2(raw counts + 1)",
    y = "Density"
  )

ggsave(
  file.path(out_dir, "qc_density_before_nonzero.png"),
  p_density_before,
  width = 8,
  height = 5
)

# ---- plot 3B: density after normalization, nonzero genes only ----
norm_dt <- as.data.table(log_norm, keep.rownames = "gene_id")
norm_long <- melt(
  norm_dt,
  id.vars = "gene_id",
  variable.name = "sample",
  value.name = "value"
)

# Remove zero-expression genes for readability
norm_long_nz <- norm_long[value > 0]

p_density_after <- ggplot(norm_long_nz, aes(x = value, color = sample)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Density curves after DESeq2 normalization",
    subtitle = "Nonzero genes only; values are log2(normalized counts + 1)",
    x = "log2(DESeq2-normalized counts + 1)",
    y = "Density"
  )

ggsave(
  file.path(out_dir, "qc_density_after_nonzero.png"),
  p_density_after,
  width = 8,
  height = 5
)

message("Normalization complete.")
message("Input counts: ", counts_path)
message("Output directory: ", out_dir)
message("Samples: ", ncol(count_mat))
message("Genes after filtering zero-total genes: ", nrow(norm_counts))