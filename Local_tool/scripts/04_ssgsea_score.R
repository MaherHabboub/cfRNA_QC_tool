# 04_ssgsea_score.R
# ssGSEA scores for platelet, erythroid/hemolysis, and endothelial contamination.
# Uses PanglaoDB marker sets + cached Ensembl HGNC-to-ENSG mapping.
# Produces score table, marker match stats, PC1 correlation plots, and score-PC heatmap.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(GSVA)
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
    c("--pca"),
    type = "character",
    default = "./results/pca/pca_scores.tsv",
    help = "PCA scores TSV produced by 03_pca_top_variable_genes.R."
  ),
  make_option(
    c("--panglao"),
    type = "character",
    default = "./resources/gene_sets/PanglaoDB_markers_27_Mar_2020.tsv",
    help = "PanglaoDB marker TSV."
  ),
  make_option(
    c("--ensembl"),
    type = "character",
    default = "./resources/ensembl_cache/hsapiens_gene_biotypes.tsv",
    help = "Cached Ensembl annotation TSV with ensembl_gene_id and hgnc_symbol."
  ),
  make_option(
    c("--out"),
    type = "character",
    default = "./results/gene_set_qc/ssgsea",
    help = "Output directory for ssGSEA scores and plots."
  ),
  make_option(
    c("--min_genes"),
    type = "integer",
    default = 5,
    help = "Minimum number of matched genes required for each gene set."
  ),
  make_option(
    c("--max_pcs"),
    type = "integer",
    default = 10,
    help = "Maximum number of PCs to test in the score-PC correlation heatmap."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

expr_path   <- opt$expr
pca_path    <- opt$pca
panglao_tsv <- opt$panglao
ensembl_map <- opt$ensembl
out_dir     <- opt$out
min_genes   <- opt$min_genes
max_pcs     <- opt$max_pcs

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- input checks ----
if (!file.exists(expr_path)) stop("Expression file not found: ", expr_path)
if (!file.exists(pca_path)) stop("PCA scores file not found: ", pca_path)
if (!file.exists(panglao_tsv)) stop("PanglaoDB file not found: ", panglao_tsv)
if (!file.exists(ensembl_map)) stop("Ensembl cache file not found: ", ensembl_map)
if (max_pcs < 1) stop("--max_pcs must be at least 1.")

# ---- load expression ----
expr_dt <- fread(expr_path)

if (!"gene_id" %in% names(expr_dt)) {
  stop("Expected a 'gene_id' column in expression file: ", expr_path)
}

sample_cols <- setdiff(names(expr_dt), "gene_id")

if (length(sample_cols) < 2) {
  stop("Need >=2 samples to relate scores to PCs. Found: ", length(sample_cols))
}

expr_mat <- as.matrix(expr_dt[, ..sample_cols])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- sub("\\..*$", "", expr_dt$gene_id)
colnames(expr_mat) <- sample_cols

if (anyNA(expr_mat)) {
  stop("NA values detected in expression matrix after numeric conversion.")
}

# ---- load PCA ----
pca <- fread(pca_path)

if (!"sample_id" %in% names(pca)) {
  stop("PCA file must contain a sample_id column.")
}

pc_cols_all <- names(pca)[grepl("^PC[0-9]+$", names(pca))]

if (length(pc_cols_all) < 1) {
  stop("No PC columns found in PCA file. Expected columns like PC1, PC2, ...")
}

pc_nums <- as.integer(sub("^PC", "", pc_cols_all))
pc_cols_all <- pc_cols_all[order(pc_nums)]
pc_cols <- pc_cols_all[seq_len(min(max_pcs, length(pc_cols_all)))]

if (!"PC1" %in% pc_cols_all) {
  stop("PCA file must contain PC1 for the individual score-vs-PC1 plots.")
}

common_samples <- intersect(pca$sample_id, colnames(expr_mat))

if (length(common_samples) < 2) {
  stop("Sample IDs in PCA and expression do not overlap enough.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
pca <- pca %>% filter(sample_id %in% common_samples)

# ---- load Ensembl map ----
map_dt <- fread(ensembl_map)

required_map_cols <- c("ensembl_gene_id", "hgnc_symbol")
missing_map_cols <- setdiff(required_map_cols, names(map_dt))

if (length(missing_map_cols) > 0) {
  stop("Ensembl map missing required columns: ", paste(missing_map_cols, collapse = ", "))
}

map_dt <- map_dt %>%
  select(ensembl_gene_id, hgnc_symbol) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct()

# ---- load PanglaoDB ----
pg <- fread(panglao_tsv)

cell_col <- names(pg)[
  grepl("cell", names(pg), ignore.case = TRUE) &
    grepl("type", names(pg), ignore.case = TRUE)
]

gene_col <- names(pg)[
  grepl("gene", names(pg), ignore.case = TRUE) |
    grepl("symbol", names(pg), ignore.case = TRUE)
]

if (length(cell_col) == 0) {
  stop("Could not find a cell-type column in PanglaoDB TSV.")
}

if (length(gene_col) == 0) {
  stop("Could not find a gene/symbol column in PanglaoDB TSV.")
}

cell_col <- cell_col[1]
gene_col <- gene_col[1]

# ---- gene set definitions ----
set_defs <- data.frame(
  set_name = c("Platelet", "Erythroid", "Endothelial"),
  pattern  = c("platelet", "erythroid|hemoglobin|erythro", "endothelial"),
  stringsAsFactors = FALSE
)

# ---- helper: extract PanglaoDB markers, map to ENSG, intersect with matrix ----
extract_gene_set_ensg <- function(pattern, set_name) {
  syms <- pg %>%
    filter(grepl(pattern, .data[[cell_col]], ignore.case = TRUE)) %>%
    pull(.data[[gene_col]]) %>%
    unique()
  
  syms <- syms[!is.na(syms) & syms != ""]
  
  if (length(syms) == 0) {
    stop("No markers found for ", set_name, " using pattern: ", pattern)
  }
  
  ensg <- map_dt %>%
    filter(hgnc_symbol %in% syms) %>%
    pull(ensembl_gene_id) %>%
    unique()
  
  genes_in_mat <- intersect(ensg, rownames(expr_mat))
  
  if (length(genes_in_mat) < min_genes) {
    stop(
      "Too few ", set_name, " genes matched in expression matrix: ",
      length(genes_in_mat), " found; minimum required is ", min_genes
    )
  }
  
  list(
    set_name = set_name,
    pattern = pattern,
    syms = syms,
    ensg = ensg,
    genes_in_mat = genes_in_mat
  )
}

sets <- lapply(seq_len(nrow(set_defs)), function(i) {
  extract_gene_set_ensg(
    pattern = set_defs$pattern[i],
    set_name = set_defs$set_name[i]
  )
})

names(sets) <- set_defs$set_name

gene_sets <- lapply(sets, function(x) x$genes_in_mat)

# ---- ssGSEA scores using new GSVA API ----
param <- GSVA::ssgseaParam(expr_mat, gene_sets, normalize = TRUE)
ss <- GSVA::gsva(param)

# ---- score table ----
score_df <- data.frame(sample_id = colnames(ss))

for (set_name in rownames(ss)) {
  col_name <- paste0(tolower(set_name), "_score_ssgsea")
  score_df[[col_name]] <- as.numeric(ss[set_name, ])
}

fwrite(
  score_df,
  file.path(out_dir, "ssgsea_scores.tsv"),
  sep = "\t"
)

# ---- marker match stats ----
meta <- data.frame(
  set = names(sets),
  panglao_pattern = sapply(sets, function(x) x$pattern),
  markers_in_panglao = sapply(sets, function(x) length(x$syms)),
  markers_mapped_to_ensg = sapply(sets, function(x) length(x$ensg)),
  markers_found_in_matrix = sapply(sets, function(x) length(x$genes_in_mat)),
  stringsAsFactors = FALSE
)

fwrite(
  meta,
  file.path(out_dir, "ssgsea_marker_match_stats.tsv"),
  sep = "\t"
)

# ---- plotting function: score vs PC1 ----
plot_score_vs_pc1 <- function(df, score_col, set_name, out_png) {
  plot_df <- pca %>%
    select(sample_id, PC1) %>%
    left_join(df, by = "sample_id")
  
  if (anyNA(plot_df[[score_col]])) {
    stop("Missing score values after joining with PCA for: ", score_col)
  }
  
  pear <- cor.test(plot_df$PC1, plot_df[[score_col]], method = "pearson")
  spear <- cor.test(plot_df$PC1, plot_df[[score_col]], method = "spearman")
  
  lab <- paste0(
    "R = ", round(unname(pear$estimate), 2),
    ", P = ", format.pval(pear$p.value, digits = 2), "\n",
    "\u03C1 = ", round(unname(spear$estimate), 2),
    ", P = ", format.pval(spear$p.value, digits = 2)
  )
  
  p <- ggplot(plot_df, aes(x = PC1, y = .data[[score_col]])) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = lab,
      hjust = 1.1,
      vjust = 1.2,
      size = 4
    ) +
    theme_minimal() +
    labs(
      x = "PC1",
      y = paste0(set_name, " score (ssGSEA)"),
      title = paste0(set_name, " score vs PC1")
    )
  
  ggsave(out_png, p, width = 7, height = 4.5)
}

# ---- make one PC1 plot per score ----
for (set_name in rownames(ss)) {
  score_col <- paste0(tolower(set_name), "_score_ssgsea")
  out_png <- file.path(
    out_dir,
    paste0(tolower(set_name), "_score_vs_PC1_ssgsea.png")
  )
  
  plot_score_vs_pc1(
    df = score_df,
    score_col = score_col,
    set_name = set_name,
    out_png = out_png
  )
}

# ---- score-PC correlation heatmap ----
score_cols <- setdiff(names(score_df), "sample_id")

score_pc_df <- score_df %>%
  inner_join(
    pca %>% select(sample_id, all_of(pc_cols)),
    by = "sample_id"
  )

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

cor_rows <- list()

for (score_col in score_cols) {
  for (pc in pc_cols) {
    pear <- safe_cor_test(score_pc_df[[score_col]], score_pc_df[[pc]], method = "pearson")
    spear <- safe_cor_test(score_pc_df[[score_col]], score_pc_df[[pc]], method = "spearman")
    
    pretty_score <- gsub("_score_ssgsea$", "", score_col)
    pretty_score <- tools::toTitleCase(pretty_score)
    
    cor_rows[[length(cor_rows) + 1]] <- data.frame(
      score = pretty_score,
      score_column = score_col,
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

score_cor_df <- rbindlist(cor_rows, fill = TRUE)

score_cor_df[, pearson_padj := p.adjust(pearson_p, method = "BH")]
score_cor_df[, spearman_padj := p.adjust(spearman_p, method = "BH")]

score_cor_df[, sig := fifelse(
  is.na(pearson_padj), "",
  fifelse(
    pearson_padj < 0.001, "***",
    fifelse(
      pearson_padj < 0.01, "**",
      fifelse(pearson_padj < 0.05, "*", "")
    )
  )
)]

score_cor_df[, label := ifelse(
  is.na(pearson_r),
  "",
  paste0(round(pearson_r, 2), sig)
)]

fwrite(
  score_cor_df,
  file.path(out_dir, "ssgsea_score_pc_correlations.tsv"),
  sep = "\t"
)

score_cor_df$score <- factor(
  score_cor_df$score,
  levels = c("Platelet", "Erythroid", "Endothelial")
)

score_cor_df$PC <- factor(score_cor_df$PC, levels = pc_cols)

# ---- heatmap with square tiles ----
p_heat <- ggplot(score_cor_df, aes(x = PC, y = score, fill = pearson_r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 4) +
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
    panel.grid = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9)
  ) +
  labs(
    caption = "* adjusted p < 0.05    ** adjusted p < 0.01    *** adjusted p < 0.001"
  )

# Dynamic dimensions so tiles remain square and readable.
# Width adapts to number of PCs; height adapts to number of score rows.
heat_width <- max(6, 0.8 * length(pc_cols) + 2)
heat_height <- max(3.5, 0.8 * length(unique(score_cor_df$score)) + 2)

ggsave(
  file.path(out_dir, "ssgsea_score_pc_correlation_heatmap.png"),
  p_heat,
  width = heat_width,
  height = heat_height
)

message("ssGSEA scoring complete.")
message("Input expression: ", expr_path)
message("Input PCA: ", pca_path)
message("Output directory: ", out_dir)
message("Gene sets scored: ", paste(rownames(ss), collapse = ", "))
message("PCs tested in score-PC heatmap: ", paste(pc_cols, collapse = ", "))
message("Heatmap size: ", round(heat_width, 2), " x ", round(heat_height, 2), " inches")