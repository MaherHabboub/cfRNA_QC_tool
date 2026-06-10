# 02_biotype_distribution.R
# Biotype distribution QC using raw HTSeq counts and cached Ensembl annotations.
# Input:
#   - HTSeq count matrix: gene_id as first column, samples as remaining columns
#   - cached Ensembl annotation table from 00_cache_ensembl_biotypes.R
# Output:
#   - biotype_distribution_counts_and_fraction.tsv
#   - biotype_distribution_percent.png, if sample count is manageable
#   - biotype_distribution_chunks/biotype_distribution_percent_part_XXX.png

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ---- command-line arguments ----
option_list <- list(
  make_option(
    c("--counts"),
    type = "character",
    default = "./data/platelet_htseq_counts_matrix.tsv",
    help = "Input HTSeq count matrix TSV. First column should be gene IDs, remaining columns samples."
  ),
  make_option(
    c("--cache"),
    type = "character",
    default = "./resources/ensembl_cache/hsapiens_gene_biotypes.tsv",
    help = "Cached Ensembl annotation TSV with ensembl_gene_id, hgnc_symbol, gene_biotype, chromosome_name."
  ),
  make_option(
    c("--out"),
    type = "character",
    default = "./results/biotype_qc",
    help = "Output directory for biotype QC results."
  ),
  make_option(
    c("--samples_per_plot"),
    type = "integer",
    default = 20,
    help = "Maximum number of samples per chunked biotype plot."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

counts_path <- opt$counts
cache_path  <- opt$cache
out_dir     <- opt$out
samples_per_plot <- opt$samples_per_plot

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(counts_path)) {
  stop("Counts file not found: ", counts_path)
}

if (!file.exists(cache_path)) {
  stop(
    "Cache file not found:\n  ", cache_path,
    "\nRun: r_scripts/00_cache_ensembl_biotypes.R"
  )
}

if (samples_per_plot < 1) {
  stop("--samples_per_plot must be at least 1.")
}

# ---- load cached annotation ----
all_genes <- fread(cache_path)

required_cols <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "chromosome_name")
missing_cols <- setdiff(required_cols, names(all_genes))
if (length(missing_cols) > 0) {
  stop("Annotation cache is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- collapse biotypes into broader bins (lab logic) ----
all_genes$gene_biotype[grep("pseudogene", all_genes$gene_biotype)] <- "pseudogene"

all_genes$gene_biotype[grep("^TR",  all_genes$gene_biotype)] <- "other"
all_genes$gene_biotype[grep("^IG",  all_genes$gene_biotype)] <- "other"
all_genes$gene_biotype[grep("^TEC", all_genes$gene_biotype)] <- "other"
all_genes$gene_biotype[grep("^sc",  all_genes$gene_biotype)] <- "other"
all_genes$gene_biotype[grep("ribozyme", all_genes$gene_biotype)] <- "other"
all_genes$gene_biotype[grep("miRNA",    all_genes$gene_biotype)] <- "other"
all_genes$gene_biotype[grep("vault",    all_genes$gene_biotype)] <- "other"

all_genes$gene_biotype[grep("antisense_RNA",          all_genes$gene_biotype)] <- "lncRNA"
all_genes$gene_biotype[grep("lincRNA",                all_genes$gene_biotype)] <- "lncRNA"
all_genes$gene_biotype[grep("lncRNA",                 all_genes$gene_biotype)] <- "lncRNA"
all_genes$gene_biotype[grep("overlapping",            all_genes$gene_biotype)] <- "lncRNA"
all_genes$gene_biotype[grep("^sense",                 all_genes$gene_biotype)] <- "lncRNA"
all_genes$gene_biotype[grep("non_coding",             all_genes$gene_biotype)] <- "lncRNA"
all_genes$gene_biotype[grep("processed_transcript",   all_genes$gene_biotype)] <- "lncRNA"

all_genes$gene_biotype[grep("^s.*RNA$", all_genes$gene_biotype)] <- "s(no)RNA"

all_genes$gene_biotype[grep("^MT$", all_genes$chromosome_name)] <- "MT_gene"
all_genes$gene_biotype[grep("^MT-RNR", all_genes$hgnc_symbol)]  <- "MT_RNRgene"

# Make sure one Ensembl gene ID maps to one row for joining
all_genes <- all_genes %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ---- load counts ----
dt <- fread(counts_path, sep = "\t", header = TRUE)

if (ncol(dt) < 2) {
  stop("Counts matrix must have at least 2 columns: gene_id + at least one sample.")
}

# Rename first column to gene_id for consistency
setnames(dt, 1, "gene_id")

# Drop HTSeq summary rows if present
dt <- dt[!grepl("^__", dt$gene_id), ]

# Remove empty/missing gene IDs
dt <- dt[!is.na(gene_id) & gene_id != "", ]

# Sample columns are all remaining columns
sample_cols <- setdiff(names(dt), "gene_id")

if (length(sample_cols) < 1) {
  stop("No sample columns found after gene_id column.")
}

# Ensure numeric count columns
for (sc in sample_cols) {
  dt[[sc]] <- suppressWarnings(as.numeric(dt[[sc]]))
}

if (anyNA(dt[, ..sample_cols])) {
  stop("NA values detected after converting count columns to numeric. Check input file formatting.")
}

# Strip Ensembl version suffix if present
dt$ensembl_gene_id <- sub("\\..*$", "", dt$gene_id)

# ---- join counts with biotypes ----
counts_bt <- dt %>%
  left_join(all_genes, by = "ensembl_gene_id")

# Extra label used in lab code/reference
counts_bt$gene_biotype[grep("45S", counts_bt$gene_id)] <- "rRNA_45S"

# Any unmatched genes become unknown
counts_bt$gene_biotype[is.na(counts_bt$gene_biotype)] <- "unknown"

# ---- summarise counts per biotype per sample ----
bt_sum <- counts_bt %>%
  group_by(gene_biotype) %>%
  summarise(
    across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample_id",
    values_to = "count"
  ) %>%
  group_by(sample_id) %>%
  mutate(frac = count / sum(count)) %>%
  ungroup()

# Keep sample order from the input matrix
bt_sum$sample_id <- factor(bt_sum$sample_id, levels = sample_cols)

# ---- save table ----
fwrite(
  as.data.table(bt_sum),
  file.path(out_dir, "biotype_distribution_counts_and_fraction.tsv"),
  sep = "\t"
)

# ---- shared plotting function ----
make_biotype_plot <- function(plot_data, title_suffix = NULL) {
  p <- ggplot(plot_data, aes(x = sample_id, y = frac, fill = gene_biotype)) +
    geom_col(width = 0.85, color = "black", linewidth = 0.2) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = element_blank()
    ) +
    labs(
      title = ifelse(
        is.null(title_suffix),
        "% counts per biotype (raw HTSeq counts)",
        paste0("% counts per biotype (raw HTSeq counts): ", title_suffix)
      )
    )
  
  p
}

# ---- plot 1: full plot only if sample count is manageable ----
# For many samples, the full plot becomes unreadable, so we still save chunked plots below.
if (length(sample_cols) <= samples_per_plot) {
  p_full <- make_biotype_plot(bt_sum)
  
  ggsave(
    file.path(out_dir, "biotype_distribution_percent.png"),
    p_full,
    width = 8,
    height = 5
  )
} else {
  message(
    "Skipping single full biotype plot because there are ",
    length(sample_cols),
    " samples. Chunked plots will be generated instead."
  )
}

# ---- plot 2: chunked plots, max N samples each ----
chunk_dir <- file.path(out_dir, "biotype_distribution_chunks")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

# Remove old chunked plots to avoid stale files from previous runs
old_chunks <- list.files(
  chunk_dir,
  pattern = "^biotype_distribution_percent_part_.*\\.png$",
  full.names = TRUE
)
if (length(old_chunks) > 0) {
  file.remove(old_chunks)
}

sample_chunks <- split(
  sample_cols,
  ceiling(seq_along(sample_cols) / samples_per_plot)
)

chunk_manifest <- data.frame(
  part = integer(),
  file = character(),
  first_sample = character(),
  last_sample = character(),
  n_samples = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(sample_chunks)) {
  chunk_samples <- sample_chunks[[i]]
  
  chunk_data <- bt_sum %>%
    filter(as.character(sample_id) %in% chunk_samples)
  
  chunk_data$sample_id <- factor(
    as.character(chunk_data$sample_id),
    levels = chunk_samples
  )
  
  title_suffix <- paste0(
    "part ", i, " of ", length(sample_chunks),
    " (", chunk_samples[1], " to ", chunk_samples[length(chunk_samples)], ")"
  )
  
  p_chunk <- make_biotype_plot(chunk_data, title_suffix = title_suffix)
  
  out_file <- file.path(
    chunk_dir,
    paste0("biotype_distribution_percent_part_", sprintf("%03d", i), ".png")
  )
  
  # Width scales mildly with number of samples, but stays manageable
  plot_width <- max(8, min(14, 0.45 * length(chunk_samples)))
  
  ggsave(
    out_file,
    p_chunk,
    width = plot_width,
    height = 5
  )
  
  chunk_manifest <- rbind(
    chunk_manifest,
    data.frame(
      part = i,
      file = out_file,
      first_sample = chunk_samples[1],
      last_sample = chunk_samples[length(chunk_samples)],
      n_samples = length(chunk_samples),
      stringsAsFactors = FALSE
    )
  )
}

fwrite(
  as.data.table(chunk_manifest),
  file.path(out_dir, "biotype_distribution_chunk_manifest.tsv"),
  sep = "\t"
)

message("Biotype QC complete.")
message("Input counts: ", counts_path)
message("Annotation cache: ", cache_path)
message("Output directory: ", out_dir)
message("Samples: ", length(sample_cols))
message("Samples per chunked plot: ", samples_per_plot)
message("Chunked plots generated: ", length(sample_chunks))
message("Biotype categories: ", length(unique(bt_sum$gene_biotype)))