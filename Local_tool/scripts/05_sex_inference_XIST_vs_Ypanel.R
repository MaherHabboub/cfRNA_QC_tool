# 05_sex_inference_XIST_vs_Ypanel.R
# Sex-chromosome QC plot using XIST expression vs mean expression of 5 Y-linked genes.
# This is a QC visualization, not a formal gender/sex classifier.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ---- command-line arguments ----
option_list <- list(
  make_option(
    c("--expr"),
    type = "character",
    default = "./results/norm_log/counts_deseq2_log2norm_plus1.tsv",
    help = "Input log-normalized expression matrix TSV with gene_id + sample columns."
  ),
  make_option(
    c("--panel"),
    type = "character",
    default = "./resources/gene_sets/sex_qc_panel_XIST_plus_5Y.tsv",
    help = "Sex QC panel TSV with hgnc_symbol and ensembl_gene_id columns."
  ),
  make_option(
    c("--out"),
    type = "character",
    default = "./results/sex_inference",
    help = "Output directory for sex inference QC results."
  ),
  make_option(
    c("--min_y_genes"),
    type = "integer",
    default = 2,
    help = "Minimum number of Y-panel genes required to compute the Y score."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

expr_path  <- opt$expr
panel_path <- opt$panel
out_dir    <- opt$out
min_y_genes <- opt$min_y_genes

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- input checks ----
if (!file.exists(expr_path)) {
  stop("Expression file not found: ", expr_path)
}

if (!file.exists(panel_path)) {
  stop("Sex QC panel file not found: ", panel_path)
}

if (min_y_genes < 1) {
  stop("--min_y_genes must be at least 1.")
}

# ---- load expression ----
expr_dt <- fread(expr_path)

if (!"gene_id" %in% names(expr_dt)) {
  stop("Expected a 'gene_id' column in expression file: ", expr_path)
}

sample_cols <- setdiff(names(expr_dt), "gene_id")

if (length(sample_cols) < 2) {
  stop("Need >=2 samples to make the sex inference plot. Found: ", length(sample_cols))
}

expr_mat <- as.matrix(expr_dt[, ..sample_cols])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- sub("\\..*$", "", expr_dt$gene_id)
colnames(expr_mat) <- sample_cols

if (anyNA(expr_mat)) {
  stop("NA values detected in expression matrix after numeric conversion.")
}

# ---- load panel ----
panel <- fread(panel_path)

required_cols <- c("hgnc_symbol", "ensembl_gene_id")
missing_cols <- setdiff(required_cols, names(panel))

if (length(missing_cols) > 0) {
  stop("Sex QC panel missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- identify XIST and Y genes ----
y_symbols <- c("EIF1AY", "KDM5D", "UTY", "DDX3Y", "RPS4Y1")

xist_ensg <- panel$ensembl_gene_id[panel$hgnc_symbol == "XIST"]
y_ensg <- panel$ensembl_gene_id[panel$hgnc_symbol %in% y_symbols]

xist_ensg <- unique(xist_ensg[!is.na(xist_ensg) & xist_ensg != ""])
y_ensg <- unique(y_ensg[!is.na(y_ensg) & y_ensg != ""])

if (length(xist_ensg) != 1) {
  stop("Could not resolve XIST to a single Ensembl gene ID from the panel file.")
}

if (length(y_ensg) < min_y_genes) {
  stop("Too few Y genes resolved from the panel file: ", length(y_ensg))
}

# ---- intersect with expression matrix ----
xist_in <- intersect(xist_ensg, rownames(expr_mat))
y_in <- intersect(y_ensg, rownames(expr_mat))

if (length(xist_in) != 1) {
  stop("XIST not found in expression matrix after version stripping.")
}

if (length(y_in) < min_y_genes) {
  stop(
    "Too few Y genes found in expression matrix: ",
    length(y_in),
    " found; minimum required is ",
    min_y_genes
  )
}

# ---- compute per-sample scores ----
xist_vec <- expr_mat[xist_in, ]
y_mean <- colMeans(expr_mat[y_in, , drop = FALSE], na.rm = TRUE)

scores <- data.frame(
  sample_id = colnames(expr_mat),
  XIST_log2norm = as.numeric(xist_vec),
  Y_panel_mean_log2norm = as.numeric(y_mean),
  Y_genes_used = length(y_in),
  stringsAsFactors = FALSE
)

fwrite(
  scores,
  file.path(out_dir, "sex_inference_scores.tsv"),
  sep = "\t"
)

# ---- marker match stats ----
matched_y_symbols <- panel %>%
  filter(ensembl_gene_id %in% y_in) %>%
  pull(hgnc_symbol) %>%
  unique()

match_stats <- data.frame(
  marker_type = c("XIST", "Y_panel"),
  markers_expected = c(1, length(y_symbols)),
  markers_resolved_to_ensg = c(length(xist_ensg), length(y_ensg)),
  markers_found_in_matrix = c(length(xist_in), length(y_in)),
  genes_found = c(
    "XIST",
    paste(matched_y_symbols, collapse = ",")
  ),
  stringsAsFactors = FALSE
)

fwrite(
  match_stats,
  file.path(out_dir, "sex_inference_marker_match_stats.tsv"),
  sep = "\t"
)

# ---- plot ----
x_max <- max(scores$Y_panel_mean_log2norm, na.rm = TRUE)
y_max <- max(scores$XIST_log2norm, na.rm = TRUE)

# Add a little padding above the maximum so labels do not get cut off
x_upper <- ifelse(is.finite(x_max) && x_max > 0, x_max * 1.10, 1)
y_upper <- ifelse(is.finite(y_max) && y_max > 0, y_max * 1.15, 1)

p <- ggplot(
  scores,
  aes(x = Y_panel_mean_log2norm, y = XIST_log2norm, label = sample_id)
) +
  geom_point(size = 2) +
  geom_text(vjust = -0.6, size = 3) +
  coord_cartesian(
    xlim = c(0, x_upper),
    ylim = c(0, y_upper),
    clip = "off"
  ) +
  theme_minimal() +
  labs(
    title = "Sex chromosome QC: XIST vs mean(Y panel)",
    x = "Mean log2(norm+1) of Y genes (EIF1AY, KDM5D, UTY, DDX3Y, RPS4Y1)",
    y = "log2(norm+1) of XIST"
  )

ggsave(
  file.path(out_dir, "sex_inference_scatter.png"),
  p,
  width = 7.5,
  height = 5
)
message("Sex chromosome QC complete.")
message("Input expression: ", expr_path)
message("Panel file: ", panel_path)
message("Output directory: ", out_dir)
message("Samples: ", length(sample_cols))
message("Y genes found/used: ", length(y_in), " / ", length(y_symbols))