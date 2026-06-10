# 03_pca_top_variable_genes.R
# PCA on top variable genes from log-normalized expression matrix.
# Input:
#   - gene_id + sample columns expression matrix
#   - optional metadata TSV with sample_id + annotation columns
# Output:
#   - pca_scores.tsv
#   - pca_scores_with_metadata.tsv, if metadata is provided
#   - pca_variance_explained.tsv
#   - pca_object.rds
#   - pca_scree.png
#   - pca_scatter_PC1_PC2.png
#   - pca_scatter_PC1_PC2_by_<metadata_column>.png, if metadata is provided

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- command-line arguments ----
option_list <- list(
  make_option(
    c("--expr"),
    type = "character",
    default = "./results/norm_log/counts_deseq2_log2norm_plus1.tsv",
    help = "Input log-normalized expression matrix TSV. Must contain gene_id column and sample columns."
  ),
  make_option(
    c("--out"),
    type = "character",
    default = "./results/pca",
    help = "Output directory for PCA results."
  ),
  make_option(
    c("--top_n"),
    type = "integer",
    default = 1000,
    help = "Number of top variable genes to use for PCA."
  ),
  make_option(
    c("--metadata"),
    type = "character",
    default = "",
    help = "Optional sample metadata TSV. Must contain sample_id column matching expression sample names."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

expr_path <- opt$expr
out_dir <- opt$out
top_n <- opt$top_n
metadata_path <- opt$metadata

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_path)) {
  stop("Expression file not found: ", expr_path)
}

if (top_n < 2) {
  stop("--top_n must be at least 2.")
}

# ---- load expression matrix ----
dt <- fread(expr_path)

if (!"gene_id" %in% names(dt)) {
  stop("Expected a 'gene_id' column in: ", expr_path)
}

gene_ids <- dt$gene_id
sample_cols <- setdiff(names(dt), "gene_id")

if (length(sample_cols) < 2) {
  stop(
    "PCA requires >= 2 samples (columns). Found: ", length(sample_cols),
    "\nYour file currently looks like a single-sample table."
  )
}

mat <- as.matrix(dt[, ..sample_cols])
rownames(mat) <- gene_ids
colnames(mat) <- sample_cols

# Ensure numeric
mode(mat) <- "numeric"

if (anyNA(mat)) {
  stop("NA values detected in expression matrix after numeric conversion.")
}

# ---- filter genes with zero variance ----
gene_var <- apply(mat, 1, var, na.rm = TRUE)
keep <- is.finite(gene_var) & gene_var > 0

mat2 <- mat[keep, , drop = FALSE]
gene_var2 <- gene_var[keep]

if (nrow(mat2) < 2) {
  stop("Not enough variable genes after filtering to run PCA.")
}

# ---- select top variable genes ----
n_use <- min(top_n, nrow(mat2))
top_genes <- names(sort(gene_var2, decreasing = TRUE))[seq_len(n_use)]
mat_top <- mat2[top_genes, , drop = FALSE]

# ---- PCA on samples ----
# center = TRUE means each gene is mean-centered across samples.
# scale. = FALSE means genes are not scaled to unit variance.
pca <- prcomp(t(mat_top), center = TRUE, scale. = FALSE)

# ---- variance explained ----
ve <- (pca$sdev^2) / sum(pca$sdev^2)

ve_df <- data.frame(
  PC = paste0("PC", seq_along(ve)),
  variance_explained = ve,
  cumulative = cumsum(ve)
)

# ---- scores/sample coordinates ----
scores <- as.data.frame(pca$x)
scores$sample_id <- rownames(scores)
scores <- scores[, c("sample_id", setdiff(names(scores), "sample_id"))]

# ---- save tables + PCA object ----
fwrite(
  ve_df,
  file.path(out_dir, "pca_variance_explained.tsv"),
  sep = "\t"
)

fwrite(
  scores,
  file.path(out_dir, "pca_scores.tsv"),
  sep = "\t"
)

saveRDS(
  pca,
  file.path(out_dir, "pca_object.rds")
)

# Also save top variable genes used
top_gene_df <- data.frame(
  gene_id = top_genes,
  variance = gene_var2[top_genes]
)

fwrite(
  top_gene_df,
  file.path(out_dir, "pca_top_variable_genes.tsv"),
  sep = "\t"
)

# ---- plot 1: scree plot ----
k <- min(10, nrow(ve_df))

ve_plot <- ve_df[1:k, ]
ve_plot$PC <- factor(ve_plot$PC, levels = paste0("PC", seq_len(k)))

p_scree <- ggplot(ve_plot, aes(x = PC, y = variance_explained)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "PCA scree (top variable genes)",
    x = NULL,
    y = "Variance explained"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "pca_scree.png"),
  p_scree,
  width = 7,
  height = 4
)

# ---- plot 2: PC1 vs PC2 scatter ----
if (!all(c("PC1", "PC2") %in% names(scores))) {
  stop("PCA did not produce both PC1 and PC2. Check number of samples.")
}

p_scatter <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 2) +
  geom_text(aes(label = sample_id), vjust = -0.7, size = 3) +
  labs(
    title = "PCA (top variable genes): PC1 vs PC2",
    x = paste0("PC1 (", round(100 * ve[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * ve[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "pca_scatter_PC1_PC2.png"),
  p_scatter,
  width = 7,
  height = 5
)

# ---- optional metadata-annotated PCA plots ----
metadata_used <- FALSE

if (!is.null(metadata_path) && metadata_path != "") {
  
  if (!file.exists(metadata_path)) {
    warning("Metadata file was provided but not found: ", metadata_path)
  } else {
    
    meta <- fread(metadata_path)
    
    if (!"sample_id" %in% names(meta)) {
      warning("Metadata file does not contain a sample_id column. Skipping metadata PCA plots.")
    } else {
      
      meta$sample_id <- as.character(meta$sample_id)
      scores$sample_id <- as.character(scores$sample_id)
      
      common_meta_samples <- intersect(scores$sample_id, meta$sample_id)
      
      if (length(common_meta_samples) < 2) {
        warning(
          "Too few overlapping samples between PCA scores and metadata. ",
          "Found: ", length(common_meta_samples),
          ". Skipping metadata PCA plots."
        )
      } else {
        
        scores_meta <- scores %>%
          data.table::as.data.table() %>%
          merge(meta, by = "sample_id", all.x = TRUE) %>%
          as.data.frame()
        
        fwrite(
          scores_meta,
          file.path(out_dir, "pca_scores_with_metadata.tsv"),
          sep = "\t"
        )
        
        meta_cols <- setdiff(names(meta), "sample_id")
        
        if (length(meta_cols) == 0) {
          warning("Metadata file contains only sample_id and no annotation columns.")
        } else {
          
          metadata_plot_dir <- file.path(out_dir, "metadata_plots")
          dir.create(metadata_plot_dir, recursive = TRUE, showWarnings = FALSE)
          
          for (meta_col in meta_cols) {
            
            # Skip metadata columns that are completely missing after merge
            if (all(is.na(scores_meta[[meta_col]]))) {
              warning("Metadata column has all NA after merging with PCA scores: ", meta_col)
              next
            }
            
            # Use color for discrete/categorical or numeric columns.
            # ggplot handles both; numeric gets continuous color scale.
            p_meta <- ggplot(scores_meta, aes(x = PC1, y = PC2, color = .data[[meta_col]])) +
              geom_point(size = 3) +
              geom_text(aes(label = sample_id), vjust = -0.7, size = 3, color = "black") +
              labs(
                title = paste0("PCA (top variable genes): colored by ", meta_col),
                x = paste0("PC1 (", round(100 * ve[1], 1), "%)"),
                y = paste0("PC2 (", round(100 * ve[2], 1), "%)"),
                color = meta_col
              ) +
              theme_minimal()
            
            safe_col <- gsub("[^A-Za-z0-9_]+", "_", meta_col)
            
            ggsave(
              file.path(metadata_plot_dir, paste0("pca_scatter_PC1_PC2_by_", safe_col, ".png")),
              p_meta,
              width = 7,
              height = 5
            )
          }
          
          metadata_used <- TRUE
        }
      }
    }
  }
}

message("PCA complete.")
message("Input expression: ", expr_path)
message("Output directory: ", out_dir)
message("Samples: ", length(sample_cols))
message("Variable genes available: ", nrow(mat2))
message("Top variable genes used: ", n_use)

if (metadata_used) {
  message("Metadata PCA plots saved to: ", file.path(out_dir, "metadata_plots"))
} else {
  message("No metadata PCA plots generated.")
}