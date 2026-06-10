# 07_hpc_metric_pc_correlation.R
# Correlates HPC QC metrics with expression PCs.
# Expected HPC summary format:
#   sample    condition    layout    numeric_qc_metric_1    numeric_qc_metric_2 ...
#
# Outputs:
#   - hpc_qc_metrics_zscore.tsv
#   - hpc_qc_pc_correlations.tsv
#   - hpc_qc_pc_correlation_heatmap.png
#   - top_metric_pc_scatterplots/

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
    c("--hpc_summary"),
    type = "character",
    default = "./data/hpc_qc_transfer_bundle/hpc_qc_summary.tsv",
    help = "HPC QC summary TSV. Must contain sample or sample_id and numeric QC metric columns."
  ),
  make_option(
    c("--pca"),
    type = "character",
    default = "./results/pca/pca_scores.tsv",
    help = "PCA scores TSV produced by 03_pca_top_variable_genes.R."
  ),
  make_option(
    c("--out"),
    type = "character",
    default = "./results/hpc_metric_pc_qc",
    help = "Output directory."
  ),
  make_option(
    c("--max_pcs"),
    type = "integer",
    default = 10,
    help = "Maximum number of PCs to test. Script will use fewer if fewer PCs exist."
  ),
  make_option(
    c("--top_scatter"),
    type = "integer",
    default = 10,
    help = "Number of strongest metric-PC associations to plot as scatterplots."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

hpc_path <- opt$hpc_summary
pca_path <- opt$pca
out_dir <- opt$out
max_pcs <- opt$max_pcs
top_scatter <- opt$top_scatter

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(hpc_path)) stop("HPC summary file not found: ", hpc_path)
if (!file.exists(pca_path)) stop("PCA scores file not found: ", pca_path)
if (max_pcs < 1) stop("--max_pcs must be at least 1.")
if (top_scatter < 1) stop("--top_scatter must be at least 1.")

# ---- load files ----
hpc <- fread(hpc_path)
pca <- fread(pca_path)

# ---- standardize sample ID columns ----
if ("sample" %in% names(hpc)) {
  setnames(hpc, "sample", "sample_id")
} else if (!"sample_id" %in% names(hpc)) {
  stop(
    "HPC summary must contain either 'sample' or 'sample_id'. Found columns: ",
    paste(names(hpc), collapse = ", ")
  )
}

if (!"sample_id" %in% names(pca)) {
  stop("PCA scores must contain a sample_id column.")
}

hpc$sample_id <- as.character(hpc$sample_id)
pca$sample_id <- as.character(pca$sample_id)

# ---- detect PC columns ----
pc_cols_all <- names(pca)[grepl("^PC[0-9]+$", names(pca))]

if (length(pc_cols_all) < 1) {
  stop("No PC columns found in PCA file. Expected columns like PC1, PC2, ...")
}

pc_nums <- as.integer(sub("^PC", "", pc_cols_all))
pc_cols_all <- pc_cols_all[order(pc_nums)]
pc_cols <- pc_cols_all[seq_len(min(max_pcs, length(pc_cols_all)))]

# ---- join by sample_id ----
common_samples <- intersect(hpc$sample_id, pca$sample_id)

if (length(common_samples) < 3) {
  stop("Need at least 3 overlapping samples between HPC summary and PCA scores. Found: ", length(common_samples))
}

hpc <- hpc %>% filter(sample_id %in% common_samples)
pca <- pca %>% filter(sample_id %in% common_samples)

# ---- detect numeric HPC metric columns ----
metadata_cols <- c("sample_id", "condition", "layout")
candidate_metric_cols <- setdiff(names(hpc), metadata_cols)

if (length(candidate_metric_cols) < 1) {
  stop("No candidate QC metric columns found after excluding metadata columns.")
}

for (col in candidate_metric_cols) {
  hpc[[col]] <- suppressWarnings(as.numeric(hpc[[col]]))
}

numeric_metric_cols <- candidate_metric_cols[
  sapply(hpc[, ..candidate_metric_cols], is.numeric)
]

metric_keep <- c()

for (m in numeric_metric_cols) {
  vals <- hpc[[m]]
  vals <- vals[is.finite(vals)]
  
  if (length(vals) >= 3 && sd(vals, na.rm = TRUE) > 0) {
    metric_keep <- c(metric_keep, m)
  }
}

metric_cols <- metric_keep

if (length(metric_cols) < 1) {
  stop("No usable numeric QC metrics found in HPC summary after filtering.")
}

# ---- build joined table ----
joined <- hpc %>%
  select(sample_id, all_of(metric_cols)) %>%
  inner_join(
    pca %>% select(sample_id, all_of(pc_cols)),
    by = "sample_id"
  )

# ---- z-score HPC metrics ----
z_df <- joined %>%
  select(sample_id, all_of(metric_cols))

for (m in metric_cols) {
  mu <- mean(z_df[[m]], na.rm = TRUE)
  sig <- sd(z_df[[m]], na.rm = TRUE)
  
  if (is.na(sig) || sig == 0) {
    z_df[[m]] <- NA_real_
  } else {
    z_df[[m]] <- (z_df[[m]] - mu) / sig
  }
}

fwrite(
  as.data.table(z_df),
  file.path(out_dir, "hpc_qc_metrics_zscore.tsv"),
  sep = "\t"
)

# ---- correlation helper ----
safe_cor_test <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 3) {
    return(data.frame(
      correlation = NA_real_,
      p_value = NA_real_,
      n_samples = sum(ok)
    ))
  }
  
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(data.frame(
      correlation = NA_real_,
      p_value = NA_real_,
      n_samples = sum(ok)
    ))
  }
  
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  
  data.frame(
    correlation = unname(ct$estimate),
    p_value = ct$p.value,
    n_samples = sum(ok)
  )
}

# ---- calculate correlations: each metric vs each PC ----
cor_rows <- list()

for (m in metric_cols) {
  for (pc in pc_cols) {
    pear <- safe_cor_test(joined[[m]], joined[[pc]], method = "pearson")
    spear <- safe_cor_test(joined[[m]], joined[[pc]], method = "spearman")
    
    cor_rows[[length(cor_rows) + 1]] <- data.frame(
      metric = m,
      PC = pc,
      pearson_r = pear$correlation,
      pearson_p = pear$p_value,
      spearman_rho = spear$correlation,
      spearman_p = spear$p_value,
      n_samples = pear$n_samples,
      stringsAsFactors = FALSE
    )
  }
}

cor_df <- rbindlist(cor_rows, fill = TRUE)

cor_df[, pearson_padj := p.adjust(pearson_p, method = "BH")]
cor_df[, spearman_padj := p.adjust(spearman_p, method = "BH")]

cor_df[, sig := fifelse(
  is.na(pearson_padj), "",
  fifelse(
    pearson_padj < 0.001, "***",
    fifelse(
      pearson_padj < 0.01, "**",
      fifelse(pearson_padj < 0.05, "*", "")
    )
  )
)]

cor_df[, label := ifelse(
  is.na(pearson_r),
  "",
  paste0(round(pearson_r, 2), sig)
)]

fwrite(
  cor_df,
  file.path(out_dir, "hpc_qc_pc_correlations.tsv"),
  sep = "\t"
)

# ---- heatmap ----
metric_order <- cor_df %>%
  group_by(metric) %>%
  summarise(max_abs_r = max(abs(pearson_r), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_abs_r) %>%
  pull(metric)

cor_df$metric <- factor(cor_df$metric, levels = metric_order)
cor_df$PC <- factor(cor_df$PC, levels = pc_cols)

p_heat <- ggplot(cor_df, aes(x = PC, y = metric, fill = pearson_r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 3) +
  coord_fixed() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9)
  ) +
  labs(
    caption = "* adjusted p < 0.05    ** adjusted p < 0.01    *** adjusted p < 0.001"
  )

# Dynamic dimensions so heatmap tiles remain square and readable.
# Width adapts to the number of PCs.
# Height adapts to the number of QC metrics.
heat_width <- max(7, 0.75 * length(pc_cols) + 3)
heat_height <- max(5, 0.35 * length(metric_cols) + 2.5)

ggsave(
  file.path(out_dir, "hpc_qc_pc_correlation_heatmap.png"),
  p_heat,
  width = heat_width,
  height = heat_height
)

# ---- top scatterplots ----
top_df_all <- cor_df %>%
  filter(!is.na(pearson_r)) %>%
  arrange(desc(abs(pearson_r)))

n_top <- min(top_scatter, nrow(top_df_all))

top_df <- top_df_all %>%
  slice_head(n = n_top)

if (nrow(top_df) == 0) {
  message("No valid metric-PC correlations available for scatterplots.")
} else {
  scatter_dir <- file.path(out_dir, "top_metric_pc_scatterplots")
  dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (i in seq_len(nrow(top_df))) {
    metric_i <- as.character(top_df$metric[i])
    pc_i <- as.character(top_df$PC[i])
    
    plot_df <- joined %>%
      select(sample_id, all_of(metric_i), all_of(pc_i))
    
    names(plot_df) <- c("sample_id", "metric_value", "pc_value")
    
    r_i <- top_df$pearson_r[i]
    p_i <- top_df$pearson_p[i]
    padj_i <- top_df$pearson_padj[i]
    
    lab <- paste0(
      "R = ", round(r_i, 2),
      "\nP = ", format.pval(p_i, digits = 2),
      "\nBH padj = ", format.pval(padj_i, digits = 2)
    )
    
    p_scatter <- ggplot(plot_df, aes(x = metric_value, y = pc_value, label = sample_id)) +
      geom_point(size = 2) +
      geom_smooth(method = "lm", se = FALSE) +
      annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.2, size = 4) +
      theme_minimal() +
      labs(
        title = paste0(metric_i, " vs ", pc_i),
        x = metric_i,
        y = pc_i
      )
    
    safe_metric <- gsub("[^A-Za-z0-9_]+", "_", metric_i)
    out_png <- file.path(
      scatter_dir,
      paste0(sprintf("%02d", i), "_", safe_metric, "_vs_", pc_i, ".png")
    )
    
    ggsave(out_png, p_scatter, width = 7, height = 4.5)
  }
}

message("HPC metric-PC correlation QC complete.")
message("HPC summary: ", hpc_path)
message("PCA scores: ", pca_path)
message("Output directory: ", out_dir)
message("Samples used: ", length(common_samples))
message("Metrics tested: ", length(metric_cols))
message("PCs tested: ", paste(pc_cols, collapse = ", "))
message("Heatmap size: ", round(heat_width, 2), " x ", round(heat_height, 2), " inches")