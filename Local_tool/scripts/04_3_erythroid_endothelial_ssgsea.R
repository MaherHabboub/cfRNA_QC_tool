# 04_3_erythroid_endothelial_ssgsea.R
# Erythroid + Endothelial ssGSEA scores (GSVA new API) + plots vs PC1

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

out_dir <- "./results/contam_qc/ssgsea"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load expression ----
expr_dt <- fread(expr_path)
stopifnot("gene_id" %in% names(expr_dt))
sample_cols <- setdiff(names(expr_dt), "gene_id")
if (length(sample_cols) < 2) stop("Need >=2 samples to relate scores to PC1.")

expr_mat <- as.matrix(expr_dt[, ..sample_cols])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- sub("\\..*$", "", expr_dt$gene_id)

# ---- load PCA ----
pca <- fread(pca_path)
stopifnot(all(c("sample_id", "PC1") %in% names(pca)))

common_samples <- intersect(pca$sample_id, colnames(expr_mat))
if (length(common_samples) < 2) stop("Sample IDs in PCA and expression do not overlap enough.")

expr_mat <- expr_mat[, common_samples, drop = FALSE]
pca <- pca %>% filter(sample_id %in% common_samples)

# ---- load Ensembl map ----
map_dt <- fread(ensembl_map) %>%
  select(ensembl_gene_id, hgnc_symbol) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct()

# ---- load PanglaoDB ----
pg <- fread(panglao_tsv)

cell_col <- names(pg)[grepl("cell", names(pg), ignore.case = TRUE) & grepl("type", names(pg), ignore.case = TRUE)]
gene_col <- names(pg)[grepl("gene", names(pg), ignore.case = TRUE) | grepl("symbol", names(pg), ignore.case = TRUE)]

if (length(cell_col) == 0) stop("Could not find a cell-type column in PanglaoDB TSV.")
if (length(gene_col) == 0) stop("Could not find a gene/symbol column in PanglaoDB TSV.")

cell_col <- cell_col[1]
gene_col <- gene_col[1]

# ---- helper to extract a gene set from PanglaoDB and map to ENSG ----
extract_gene_set_ensg <- function(pattern, set_name) {
  syms <- pg %>%
    filter(grepl(pattern, .data[[cell_col]], ignore.case = TRUE)) %>%
    pull(.data[[gene_col]]) %>%
    unique()
  
  syms <- syms[!is.na(syms) & syms != ""]
  if (length(syms) == 0) stop("No markers found for ", set_name, " using pattern: ", pattern)
  
  ensg <- map_dt %>%
    filter(hgnc_symbol %in% syms) %>%
    pull(ensembl_gene_id) %>%
    unique()
  
  genes_in_mat <- intersect(ensg, rownames(expr_mat))
  if (length(genes_in_mat) < 5) stop("Too few ", set_name, " genes matched in matrix: ", length(genes_in_mat))
  
  list(syms = syms, ensg = ensg, genes_in_mat = genes_in_mat)
}

# Patterns (adjust if PanglaoDB labels differ)
ery_pat <- "erythroid|hemoglobin|erythro"
end_pat <- "endothelial"

ery <- extract_gene_set_ensg(ery_pat, "Erythroid")
end <- extract_gene_set_ensg(end_pat, "Endothelial")

gene_sets <- list(
  Erythroid = ery$genes_in_mat,
  Endothelial = end$genes_in_mat
)

# ---- ssGSEA scores (NEW GSVA API) ----
param <- GSVA::ssgseaParam(expr_mat, gene_sets, normalize = TRUE)
ss <- GSVA::gsva(param)  # matrix: gene_set x samples

score_df <- data.frame(
  sample_id = colnames(ss),
  erythroid_score_ssgsea   = as.numeric(ss["Erythroid", ]),
  endothelial_score_ssgsea = as.numeric(ss["Endothelial", ])
)

fwrite(score_df, file.path(out_dir, "erythroid_endothelial_scores_ssgsea.tsv"), sep = "\t")

# ---- plotting function (paper-style) ----
plot_score_vs_pc1 <- function(df, score_col, title, ylab, out_png) {
  plot_df <- pca %>%
    select(sample_id, PC1) %>%
    left_join(df, by = "sample_id")
  
  pear <- cor.test(plot_df$PC1, plot_df[[score_col]], method = "pearson")
  spear <- cor.test(plot_df$PC1, plot_df[[score_col]], method = "spearman")
  
  lab <- paste0(
    "R = ", round(unname(pear$estimate), 2), ", P = ", format.pval(pear$p.value, digits = 2), "\n",
    "\u03C1 = ", round(unname(spear$estimate), 2), ", P = ", format.pval(spear$p.value, digits = 2)
  )
  
  p <- ggplot(plot_df, aes(x = PC1, y = .data[[score_col]])) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.2, size = 4) +
    theme_minimal() +
    labs(x = "PC1", y = ylab, title = title)
  
  ggsave(out_png, p, width = 7, height = 4.5)
}

plot_score_vs_pc1(score_df, "erythroid_score_ssgsea",
                  "Erythroid score vs PC1", "Erythroid score (ssGSEA)",
                  file.path(out_dir, "erythroid_score_vs_PC1_ssgsea.png"))

plot_score_vs_pc1(score_df, "endothelial_score_ssgsea",
                  "Endothelial score vs PC1", "Endothelial score (ssGSEA)",
                  file.path(out_dir, "endothelial_score_vs_PC1_ssgsea.png"))

# ---- marker match stats ----
meta <- data.frame(
  set = c("Erythroid", "Endothelial"),
  markers_in_panglao = c(length(ery$syms), length(end$syms)),
  markers_mapped_to_ensg = c(length(ery$ensg), length(end$ensg)),
  markers_found_in_matrix = c(length(ery$genes_in_mat), length(end$genes_in_mat))
)

fwrite(meta, file.path(out_dir, "marker_match_stats.tsv"), sep = "\t")