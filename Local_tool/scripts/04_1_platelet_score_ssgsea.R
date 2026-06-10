# 04_01_platelet_score_ssgsea.R
# Platelet ssGSEA score (GSVA new API) + plot vs PC1

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(GSVA)
})

# ---- inputs ----
expr_path   <- "./results/norm_log/counts_deseq2_log2norm_plus1.tsv"
pca_path    <- "./results/pca/pca_scores.tsv"
panglao_tsv <- "./resources/gene_sets/PanglaoDB_markers_27_Mar_2020.tsv"
ensembl_map <- "./resources/ensembl_cache/hsapiens_gene_biotypes.tsv"

out_dir <- "./results/platelet_qc/ssgsea"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load expression (gene_id + sample cols) ----
expr_dt <- fread(expr_path)
stopifnot("gene_id" %in% names(expr_dt))
sample_cols <- setdiff(names(expr_dt), "gene_id")
if (length(sample_cols) < 2) stop("Need >=2 samples to relate score to PC1.")

expr_mat <- as.matrix(expr_dt[, ..sample_cols])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- sub("\\..*$", "", expr_dt$gene_id)  # strip ENSG versions if present

# ---- load PCA scores ----
pca <- fread(pca_path)
stopifnot(all(c("sample_id", "PC1") %in% names(pca)))

# Ensure sample IDs match expression columns
common_samples <- intersect(pca$sample_id, colnames(expr_mat))
if (length(common_samples) < 2) stop("Sample IDs in PCA and expression do not overlap enough.")
expr_mat <- expr_mat[, common_samples, drop = FALSE]
pca <- pca %>% filter(sample_id %in% common_samples)

# ---- load Ensembl map (HGNC <-> ENSG) ----
map_dt <- fread(ensembl_map) %>%
  select(ensembl_gene_id, hgnc_symbol) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct()

# ---- load PanglaoDB and extract platelet gene symbols robustly ----
pg <- fread(panglao_tsv)

# Try to find a "cell type" column and a "gene symbol" column without hardcoding exact names
cell_col <- names(pg)[grepl("cell", names(pg), ignore.case = TRUE) & grepl("type", names(pg), ignore.case = TRUE)]
gene_col <- names(pg)[grepl("gene", names(pg), ignore.case = TRUE) | grepl("symbol", names(pg), ignore.case = TRUE)]

if (length(cell_col) == 0) stop("Could not find a cell-type column in PanglaoDB TSV.")
if (length(gene_col) == 0) stop("Could not find a gene/symbol column in PanglaoDB TSV.")

cell_col <- cell_col[1]
gene_col <- gene_col[1]

platelet_syms <- pg %>%
  filter(grepl("platelet", .data[[cell_col]], ignore.case = TRUE)) %>%
  pull(.data[[gene_col]]) %>%
  unique()

platelet_syms <- platelet_syms[!is.na(platelet_syms) & platelet_syms != ""]
if (length(platelet_syms) == 0) stop("No platelet markers found in PanglaoDB TSV (check file/columns).")

# Map HGNC symbols -> Ensembl IDs
platelet_ensg <- map_dt %>%
  filter(hgnc_symbol %in% platelet_syms) %>%
  pull(ensembl_gene_id) %>%
  unique()

# Intersect with expression matrix genes
genes_in_mat <- intersect(platelet_ensg, rownames(expr_mat))
if (length(genes_in_mat) < 5) stop("Too few platelet genes matched in expression matrix: ", length(genes_in_mat))

# ---- ssGSEA score (NEW GSVA API) ----
gene_sets <- list(Platelet = genes_in_mat)

# GSVA >= 2.x uses method-specific parameter objects
param <- GSVA::ssgseaParam(expr_mat, gene_sets, normalize = TRUE)
ss <- GSVA::gsva(param)  # returns a matrix: gene_set x samples

score_df <- data.frame(
  sample_id = colnames(ss),
  platelet_score_ssgsea = as.numeric(ss["Platelet", ])
)

# ---- save scores ----
fwrite(score_df, file.path(out_dir, "platelet_scores_ssgsea.tsv"), sep = "\t")

# ---- merge with PC1 + compute correlations ----
plot_df <- pca %>%
  select(sample_id, PC1) %>%
  left_join(score_df, by = "sample_id")

pear <- cor.test(plot_df$PC1, plot_df$platelet_score_ssgsea, method = "pearson")
spear <- cor.test(plot_df$PC1, plot_df$platelet_score_ssgsea, method = "spearman")

lab <- paste0(
  "R = ", round(unname(pear$estimate), 2), ", P = ", format.pval(pear$p.value, digits = 2), "\n",
  "\u03C1 = ", round(unname(spear$estimate), 2), ", P = ", format.pval(spear$p.value, digits = 2)
)

p <- ggplot(plot_df, aes(x = PC1, y = platelet_score_ssgsea)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.2, size = 4) +
  theme_minimal() +
  labs(x = "PC1", y = "Platelet score (ssGSEA)", title = "Platelet score vs PC1")

ggsave(file.path(out_dir, "platelet_score_vs_PC1_ssgsea.png"), p, width = 7, height = 4.5)

# ---- marker match stats ----
meta <- data.frame(
  platelet_markers_in_panglao = length(platelet_syms),
  platelet_markers_mapped_to_ensg = length(platelet_ensg),
  platelet_markers_found_in_matrix = length(genes_in_mat)
)
fwrite(meta, file.path(out_dir, "platelet_marker_match_stats.tsv"), sep = "\t")

