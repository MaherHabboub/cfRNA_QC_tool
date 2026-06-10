# run_all_local_qc.R
# Master runner for local count-based QC.
# Runs all R QC scripts in order and records pass/fail status.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---- command-line arguments ----
option_list <- list(
  make_option(
    c("--counts"),
    type = "character",
    default = "./data/cohortA_10samples_htseq_counts_combined_all.tsv",
    help = "Input HTSeq count matrix TSV. First column = gene IDs, remaining columns = samples."
  ),
  make_option(
    c("--metadata"),
    type = "character",
    default = "",
    help = "Optional sample metadata TSV. Must contain sample_id column matching count matrix sample names."
  ),
  make_option(
    c("--ensembl"),
    type = "character",
    default = "./resources/ensembl_cache/hsapiens_gene_biotypes.tsv",
    help = "Cached Ensembl annotation TSV."
  ),
  make_option(
    c("--panglao"),
    type = "character",
    default = "./resources/gene_sets/PanglaoDB_markers_27_Mar_2020.tsv",
    help = "PanglaoDB marker TSV."
  ),
  make_option(
    c("--sex_panel"),
    type = "character",
    default = "./resources/gene_sets/sex_qc_panel_XIST_plus_5Y.tsv",
    help = "Sex QC panel TSV with XIST + Y genes."
  ),
  make_option(
    c("--hpc_summary"),
    type = "character",
    default = "./data/cohort/hpc_qc_summary.tsv",
    help = "Aggregated HPC QC summary TSV."
  ),
  make_option(
    c("--out_root"),
    type = "character",
    default = "./results/Cohort_A/",
    help = "Root output directory."
  ),
  make_option(
    c("--top_n"),
    type = "integer",
    default = 1000,
    help = "Number of top variable genes for PCA."
  ),
  make_option(
    c("--max_pcs"),
    type = "integer",
    default = 5,
    help = "Maximum number of PCs to test against HPC QC metrics."
  ),
  make_option(
    c("--top_scatter"),
    type = "integer",
    default = 10,
    help = "Number of strongest HPC metric-PC correlations to plot."
  ),
  make_option(
    c("--report_title"),
    type = "character",
    default = "Local RNA-seq QC Report",
    help = "Title used in the generated HTML report."
  ),
  make_option(
    c("--stop_on_fail"),
    action = "store_true",
    default = FALSE,
    help = "Stop pipeline immediately if one script fails."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

counts_path   <- opt$counts
metadata_path <- opt$metadata
ensembl_map   <- opt$ensembl
panglao_tsv   <- opt$panglao
sex_panel     <- opt$sex_panel
hpc_summary   <- opt$hpc_summary
out_root      <- opt$out_root
top_n         <- opt$top_n
max_pcs       <- opt$max_pcs
top_scatter   <- opt$top_scatter
report_title  <- opt$report_title
stop_on_fail  <- opt$stop_on_fail

# ---- output folders ----
norm_dir        <- file.path(out_root, "norm_log")
biotype_dir     <- file.path(out_root, "biotype_qc")
pca_dir         <- file.path(out_root, "pca")
ssgsea_dir      <- file.path(out_root, "gene_set_qc", "ssgsea")
sex_dir         <- file.path(out_root, "sex_inference")
hpc_metric_dir  <- file.path(out_root, "hpc_metric_pc_qc")
report_dir      <- file.path(out_root, "report")
log_dir         <- file.path(out_root, "run_logs")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(log_dir, "local_qc_run_status.tsv")
log_path    <- file.path(log_dir, "local_qc_run_console.log")

# Clear old log if it exists
if (file.exists(log_path)) file.remove(log_path)

# ---- PCA args, with optional metadata ----
pca_args <- c(
  "--expr", file.path(norm_dir, "counts_deseq2_log2norm_plus1.tsv"),
  "--out", pca_dir,
  "--top_n", as.character(top_n)
)

if (!is.null(metadata_path) && metadata_path != "") {
  pca_args <- c(pca_args, "--metadata", metadata_path)
}

# ---- script paths ----
scripts <- list(
  list(
    name = "01_normalize_log",
    path = "./r_scripts/01_normalize_log.R",
    args = c(
      "--counts", counts_path,
      "--out", norm_dir
    )
  ),
  list(
    name = "02_biotype_distribution",
    path = "./r_scripts/02_biotype_distribution.R",
    args = c(
      "--counts", counts_path,
      "--cache", ensembl_map,
      "--out", biotype_dir
    )
  ),
  list(
    name = "03_pca_top_variable_genes",
    path = "./r_scripts/03_pca_top_variable_genes.R",
    args = pca_args
  ),
  list(
    name = "04_ssgsea_score",
    path = "./r_scripts/04_ssgsea_score.R",
    args = c(
      "--expr", file.path(norm_dir, "counts_deseq2_log2norm_plus1.tsv"),
      "--pca", file.path(pca_dir, "pca_scores.tsv"),
      "--panglao", panglao_tsv,
      "--ensembl", ensembl_map,
      "--out", ssgsea_dir,
      "--max_pcs", as.character(max_pcs)
    )
  ),
  list(
    name = "05_sex_inference_XIST_vs_Ypanel",
    path = "./r_scripts/05_sex_inference_XIST_vs_Ypanel.R",
    args = c(
      "--expr", file.path(norm_dir, "counts_deseq2_log2norm_plus1.tsv"),
      "--panel", sex_panel,
      "--out", sex_dir
    )
  ),
  list(
    name = "07_hpc_metric_pc_correlation",
    path = "./r_scripts/07_hpc_metric_pc_correlation.R",
    args = c(
      "--hpc_summary", hpc_summary,
      "--pca", file.path(pca_dir, "pca_scores.tsv"),
      "--out", hpc_metric_dir,
      "--max_pcs", as.character(max_pcs),
      "--top_scatter", as.character(top_scatter)
    )
  ),
  list(
    name = "08_generate_html_report",
    path = "./r_scripts/08_generate_html_report.R",
    args = c(
      "--results", out_root,
      "--template", "./r_scripts/08_report_template.Rmd",
      "--out", report_dir,
      "--title", report_title
    )
  )
)

# ---- helper to run one script ----
run_one_script <- function(step) {
  start_time <- Sys.time()
  
  script_name <- step$name
  script_path <- step$path
  script_args <- step$args
  
  cat(
    paste0(
      "\n\n==============================\n",
      script_name, "\n",
      script_path, "\n",
      "Started: ", start_time, "\n",
      "==============================\n"
    ),
    file = log_path,
    append = TRUE
  )
  
  if (!file.exists(script_path)) {
    end_time <- Sys.time()
    
    msg <- paste("Script not found:", script_path)
    
    cat(msg, "\n", file = log_path, append = TRUE)
    
    return(data.frame(
      step = script_name,
      script = script_path,
      status = "MISSING",
      start_time = as.character(start_time),
      end_time = as.character(end_time),
      error_message = msg,
      stringsAsFactors = FALSE
    ))
  }
  
  command_string <- paste("Rscript", paste(shQuote(c(script_path, script_args)), collapse = " "))
  cat("Command:\n", command_string, "\n\n", file = log_path, append = TRUE)
  
  message("Running: ", script_name)
  
  result <- system2(
    command = "Rscript",
    args = shQuote(c(script_path, script_args)),
    stdout = TRUE,
    stderr = TRUE
  )
  
  exit_status <- attr(result, "status")
  if (is.null(exit_status)) exit_status <- 0
  
  cat(paste(result, collapse = "\n"), "\n", file = log_path, append = TRUE)
  
  end_time <- Sys.time()
  
  if (exit_status == 0) {
    status <- "PASS"
    error_message <- ""
  } else {
    status <- "FAIL"
    error_message <- paste(result, collapse = " | ")
  }
  
  data.frame(
    step = script_name,
    script = script_path,
    status = status,
    start_time = as.character(start_time),
    end_time = as.character(end_time),
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

# ---- input summary ----
cat(
  paste0(
    "Local QC master run\n",
    "Counts: ", counts_path, "\n",
    "Metadata: ", ifelse(metadata_path == "", "not provided", metadata_path), "\n",
    "Ensembl cache: ", ensembl_map, "\n",
    "PanglaoDB: ", panglao_tsv, "\n",
    "Sex panel: ", sex_panel, "\n",
    "HPC summary: ", hpc_summary, "\n",
    "Output root: ", out_root, "\n",
    "Report output: ", report_dir, "\n",
    "Report title: ", report_title, "\n",
    "Top variable genes: ", top_n, "\n",
    "Max PCs for correlations: ", max_pcs, "\n",
    "Top scatterplots: ", top_scatter, "\n",
    "Stop on fail: ", stop_on_fail, "\n"
  ),
  file = log_path,
  append = TRUE
)

# ---- run scripts ----
status_list <- list()

for (step in scripts) {
  res <- run_one_script(step)
  status_list[[length(status_list) + 1]] <- res
  
  status_df_so_far <- rbindlist(status_list, fill = TRUE)
  fwrite(status_df_so_far, status_path, sep = "\t")
  
  if (res$status %in% c("FAIL", "MISSING") && stop_on_fail) {
    message("Stopping after failure in: ", step$name)
    break
  }
}

status_df <- rbindlist(status_list, fill = TRUE)
fwrite(status_df, status_path, sep = "\t")

# ---- final summary ----
message("Local QC run complete.")
message("Status table: ", status_path)
message("Console log: ", log_path)
message("HTML report should be here: ", file.path(report_dir, "local_qc_report.html"))

print(status_df)

failed <- status_df[status %in% c("FAIL", "MISSING")]

if (nrow(failed) > 0) {
  message("Some steps failed or were missing:")
  print(failed[, .(step, status, error_message)])
} else {
  message("All steps completed successfully.")
}