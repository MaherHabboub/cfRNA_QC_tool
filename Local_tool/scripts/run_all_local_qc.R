# run_all_local_qc.R
# Master runner for local count-based QC.
# Schedules existing analysis scripts, records pass/fail/skipped status, and
# renders the HTML report. Analysis logic remains in the individual scripts.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---- locate the bundled scripts and resources reliably ----
command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- command_args[grepl("^--file=", command_args)]

if (length(file_arg) != 1) {
  stop("Could not determine the location of run_all_local_qc.R.")
}

runner_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
scripts_dir <- dirname(runner_path)
local_tool_dir <- dirname(scripts_dir)

# ---- command-line arguments ----
option_list <- list(
  make_option(
    c("--config"),
    type = "character",
    default = "",
    help = "Path to an R config file. Start from Local_tool/config/local_qc_config.example.R."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$config) || opt$config == "") {
  stop("--config is required. Copy and edit Local_tool/config/local_qc_config.example.R.")
}

config_path <- normalizePath(opt$config, mustWork = FALSE)

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}

config_env <- new.env(parent = baseenv())
source(config_path, local = config_env)

config_value <- function(name, default = NULL, required = FALSE) {
  if (exists(name, envir = config_env, inherits = FALSE)) {
    return(get(name, envir = config_env, inherits = FALSE))
  }

  if (required) {
    stop("Config is missing required setting: ", name)
  }

  default
}

validate_scalar_string <- function(value, name, allow_empty = FALSE) {
  if (!is.character(value) || length(value) != 1 || is.na(value) || (!allow_empty && value == "")) {
    stop(name, " must be a single", if (allow_empty) " character string." else " non-empty character string.")
  }
  value
}

validate_scalar_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop(name, " must be TRUE or FALSE.")
  }
  value
}

validate_positive_integer <- function(value, name) {
  if (length(value) != 1 || is.na(value) || !is.numeric(value) || value < 1 || value != as.integer(value)) {
    stop(name, " must be a positive integer.")
  }
  as.integer(value)
}

# ---- config values ----
counts_path <- validate_scalar_string(config_value("COUNTS", required = TRUE), "COUNTS")
metadata_path <- validate_scalar_string(config_value("METADATA", ""), "METADATA", allow_empty = TRUE)
hpc_summary <- validate_scalar_string(config_value("HPC_SUMMARY", ""), "HPC_SUMMARY", allow_empty = TRUE)
out_root <- validate_scalar_string(config_value("OUT_ROOT", required = TRUE), "OUT_ROOT")
top_n <- validate_positive_integer(config_value("TOP_N", 1000L), "TOP_N")
max_pcs <- validate_positive_integer(config_value("MAX_PCS", 5L), "MAX_PCS")
top_scatter <- validate_positive_integer(config_value("TOP_SCATTER", 10L), "TOP_SCATTER")
report_title <- validate_scalar_string(config_value("REPORT_TITLE", "Local RNA-seq QC Report"), "REPORT_TITLE")
stop_on_fail <- validate_scalar_logical(config_value("STOP_ON_FAIL", FALSE), "STOP_ON_FAIL")

biotype_enabled <- validate_scalar_logical(config_value("BIOTYPE_ENABLED", TRUE), "BIOTYPE_ENABLED")
pca_ssgsea_enabled <- validate_scalar_logical(config_value("PCA_SSGSEA_ENABLED", TRUE), "PCA_SSGSEA_ENABLED")
sex_inference_enabled <- validate_scalar_logical(config_value("SEX_INFERENCE_ENABLED", TRUE), "SEX_INFERENCE_ENABLED")

# Bundled annotation and marker files are intentionally not user inputs.
ensembl_map <- file.path(local_tool_dir, "resources", "ensembl_cache", "hsapiens_gene_biotypes.tsv")
panglao_tsv <- file.path(local_tool_dir, "resources", "gene_sets", "PanglaoDB_markers_27_Mar_2020.tsv")
sex_panel <- file.path(local_tool_dir, "resources", "gene_sets", "sex_qc_panel_XIST_plus_5Y.tsv")

# ---- output folders ----
norm_dir <- file.path(out_root, "norm_log")
biotype_dir <- file.path(out_root, "biotype_qc")
pca_dir <- file.path(out_root, "pca")
ssgsea_dir <- file.path(out_root, "gene_set_qc", "ssgsea")
sex_dir <- file.path(out_root, "sex_inference")
hpc_metric_dir <- file.path(out_root, "hpc_metric_pc_qc")
report_dir <- file.path(out_root, "report")
log_dir <- file.path(out_root, "run_logs")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(log_dir, "local_qc_run_status.tsv")
manifest_path <- file.path(log_dir, "local_qc_module_manifest.tsv")
log_path <- file.path(log_dir, "local_qc_run_console.log")

if (file.exists(log_path)) file.remove(log_path)

hpc_metric_enabled <- hpc_summary != "" && pca_ssgsea_enabled
hpc_metric_reason <- if (hpc_summary == "") {
  "No HPC_SUMMARY path configured."
} else if (!pca_ssgsea_enabled) {
  "PCA_SSGSEA_ENABLED is FALSE; HPC metric correlation requires PCA scores."
} else {
  "Enabled because HPC_SUMMARY was supplied."
}

manifest <- data.table(
  module = c("normalization", "biotype", "pca_ssgsea", "sex_inference", "hpc_metric", "report"),
  enabled = c(TRUE, biotype_enabled, pca_ssgsea_enabled, sex_inference_enabled, hpc_metric_enabled, TRUE),
  reason = c(
    "Mandatory module.",
    if (biotype_enabled) "Enabled by BIOTYPE_ENABLED." else "Disabled by BIOTYPE_ENABLED.",
    if (pca_ssgsea_enabled) "Enabled by PCA_SSGSEA_ENABLED." else "Disabled by PCA_SSGSEA_ENABLED.",
    if (sex_inference_enabled) "Enabled by SEX_INFERENCE_ENABLED." else "Disabled by SEX_INFERENCE_ENABLED.",
    hpc_metric_reason,
    "Mandatory module."
  )
)
fwrite(manifest, manifest_path, sep = "\t")

# ---- PCA args, with optional metadata ----
pca_args <- c(
  "--expr", file.path(norm_dir, "counts_deseq2_log2norm_plus1.tsv"),
  "--out", pca_dir,
  "--top_n", as.character(top_n)
)

if (metadata_path != "") {
  pca_args <- c(pca_args, "--metadata", metadata_path)
}

make_step <- function(name, script, args) {
  list(name = name, path = file.path(scripts_dir, script), args = args, skip_reason = NULL)
}

make_skipped_step <- function(name, script, reason) {
  list(name = name, path = file.path(scripts_dir, script), args = character(), skip_reason = reason)
}

# Keep arguments passed to enabled analysis scripts unchanged from the prior runner.
pipeline_steps <- list(
  make_step(
    "01_normalize_log",
    "01_normalize_log.R",
    c("--counts", counts_path, "--out", norm_dir)
  ),
  if (biotype_enabled) {
    make_step(
      "02_biotype_distribution",
      "02_biotype_distribution.R",
      c("--counts", counts_path, "--cache", ensembl_map, "--out", biotype_dir)
    )
  } else {
    make_skipped_step("02_biotype_distribution", "02_biotype_distribution.R", "Disabled by BIOTYPE_ENABLED.")
  },
  if (pca_ssgsea_enabled) {
    make_step("03_pca_top_variable_genes", "03_pca_top_variable_genes.R", pca_args)
  } else {
    make_skipped_step("03_pca_top_variable_genes", "03_pca_top_variable_genes.R", "Disabled by PCA_SSGSEA_ENABLED.")
  },
  if (pca_ssgsea_enabled) {
    make_step(
      "04_ssgsea_score",
      "04_ssgsea_score.R",
      c(
        "--expr", file.path(norm_dir, "counts_deseq2_log2norm_plus1.tsv"),
        "--pca", file.path(pca_dir, "pca_scores.tsv"),
        "--panglao", panglao_tsv,
        "--ensembl", ensembl_map,
        "--out", ssgsea_dir,
        "--max_pcs", as.character(max_pcs)
      )
    )
  } else {
    make_skipped_step("04_ssgsea_score", "04_ssgsea_score.R", "Disabled by PCA_SSGSEA_ENABLED.")
  },
  if (sex_inference_enabled) {
    make_step(
      "05_sex_inference_XIST_vs_Ypanel",
      "05_sex_inference_XIST_vs_Ypanel.R",
      c("--expr", file.path(norm_dir, "counts_deseq2_log2norm_plus1.tsv"), "--panel", sex_panel, "--out", sex_dir)
    )
  } else {
    make_skipped_step("05_sex_inference_XIST_vs_Ypanel", "05_sex_inference_XIST_vs_Ypanel.R", "Disabled by SEX_INFERENCE_ENABLED.")
  },
  if (hpc_metric_enabled) {
    make_step(
      "07_hpc_metric_pc_correlation",
      "07_hpc_metric_pc_correlation.R",
      c(
        "--hpc_summary", hpc_summary,
        "--pca", file.path(pca_dir, "pca_scores.tsv"),
        "--out", hpc_metric_dir,
        "--max_pcs", as.character(max_pcs),
        "--top_scatter", as.character(top_scatter)
      )
    )
  } else {
    make_skipped_step("07_hpc_metric_pc_correlation", "07_hpc_metric_pc_correlation.R", hpc_metric_reason)
  }
)

report_step <- make_step(
  "08_generate_html_report",
  "08_generate_html_report.R",
  c(
    "--results", out_root,
    "--template", file.path(scripts_dir, "08_report_template.Rmd"),
    "--out", report_dir,
    "--title", report_title,
    "--manifest", manifest_path
  )
)

status_row <- function(step, status, error_message = "") {
  now <- as.character(Sys.time())
  data.frame(
    step = step$name,
    script = step$path,
    status = status,
    start_time = now,
    end_time = now,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

run_one_script <- function(step) {
  start_time <- Sys.time()

  cat(
    paste0(
      "\n\n==============================\n",
      step$name, "\n",
      step$path, "\n",
      "Started: ", start_time, "\n",
      "==============================\n"
    ),
    file = log_path,
    append = TRUE
  )

  if (!file.exists(step$path)) {
    msg <- paste("Script not found:", step$path)
    cat(msg, "\n", file = log_path, append = TRUE)
    return(status_row(step, "MISSING", msg))
  }

  command_string <- paste("Rscript", paste(shQuote(c(step$path, step$args)), collapse = " "))
  cat("Command:\n", command_string, "\n\n", file = log_path, append = TRUE)

  message("Running: ", step$name)
  result <- system2(
    command = "Rscript",
    args = shQuote(c(step$path, step$args)),
    stdout = TRUE,
    stderr = TRUE
  )

  exit_status <- attr(result, "status")
  if (is.null(exit_status)) exit_status <- 0

  cat(paste(result, collapse = "\n"), "\n", file = log_path, append = TRUE)

  end_time <- Sys.time()
  status <- if (exit_status == 0) "PASS" else "FAIL"
  error_message <- if (exit_status == 0) "" else paste(result, collapse = " | ")

  data.frame(
    step = step$name,
    script = step$path,
    status = status,
    start_time = as.character(start_time),
    end_time = as.character(end_time),
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

write_status <- function(status_list) {
  fwrite(rbindlist(status_list, fill = TRUE), status_path, sep = "\t")
}

# ---- input summary ----
cat(
  paste0(
    "Local QC master run\n",
    "Config: ", config_path, "\n",
    "Counts: ", counts_path, "\n",
    "Metadata: ", ifelse(metadata_path == "", "not provided", metadata_path), "\n",
    "Ensembl cache: ", ensembl_map, "\n",
    "PanglaoDB: ", panglao_tsv, "\n",
    "Sex panel: ", sex_panel, "\n",
    "HPC summary: ", ifelse(hpc_summary == "", "not provided", hpc_summary), "\n",
    "Output root: ", out_root, "\n",
    "Report output: ", report_dir, "\n",
    "Report title: ", report_title, "\n",
    "Top variable genes: ", top_n, "\n",
    "Max PCs for correlations: ", max_pcs, "\n",
    "Top scatterplots: ", top_scatter, "\n",
    "BIOTYPE_ENABLED: ", biotype_enabled, "\n",
    "PCA_SSGSEA_ENABLED: ", pca_ssgsea_enabled, "\n",
    "SEX_INFERENCE_ENABLED: ", sex_inference_enabled, "\n",
    "Stop on fail: ", stop_on_fail, "\n"
  ),
  file = log_path,
  append = TRUE
)

# ---- run configured analysis steps ----
status_list <- list()
stop_remaining_analysis <- FALSE

for (step in pipeline_steps) {
  if (!is.null(step$skip_reason)) {
    res <- status_row(step, "SKIPPED", step$skip_reason)
  } else if (stop_remaining_analysis) {
    res <- status_row(step, "SKIPPED", "Skipped because STOP_ON_FAIL stopped remaining analysis steps.")
  } else {
    res <- run_one_script(step)
  }

  status_list[[length(status_list) + 1]] <- res
  write_status(status_list)

  if (res$status %in% c("FAIL", "MISSING") && stop_on_fail) {
    message("Stopping remaining analysis steps after failure in: ", step$name)
    stop_remaining_analysis <- TRUE
  }
}

# The report is always attempted, including after an analysis failure.
write_status(status_list)
report_res <- run_one_script(report_step)
status_list[[length(status_list) + 1]] <- report_res
write_status(status_list)

# ---- final summary ----
status_df <- rbindlist(status_list, fill = TRUE)

message("Local QC run complete.")
message("Module manifest: ", manifest_path)
message("Status table: ", status_path)
message("Console log: ", log_path)
message("HTML report should be here: ", file.path(report_dir, "local_qc_report.html"))

print(status_df)

failed <- status_df[status %in% c("FAIL", "MISSING")]

if (nrow(failed) > 0) {
  message("Some steps failed or were missing:")
  print(failed[, .(step, status, error_message)])
} else {
  message("All scheduled steps completed successfully.")
}
