# 08_generate_html_report.R
# Renders the local QC HTML report from the R Markdown template.

suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
})

# ---- command-line arguments ----
option_list <- list(
  make_option(
    c("--results"),
    type = "character",
    default = "./results",
    help = "Root results directory produced by the QC pipeline."
  ),
  make_option(
    c("--template"),
    type = "character",
    default = "./r_scripts/08_report_template.Rmd",
    help = "Path to the R Markdown report template."
  ),
  make_option(
    c("--out"),
    type = "character",
    default = "./results/report",
    help = "Output directory for the HTML report."
  ),
  make_option(
    c("--title"),
    type = "character",
    default = "Local RNA-seq QC Report",
    help = "Report title."
  ),
  make_option(
    c("--manifest"),
    type = "character",
    default = "",
    help = "Optional module manifest TSV produced by run_all_local_qc.R."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

results_dir <- opt$results
template_path <- opt$template
out_dir <- opt$out
report_title <- opt$title
manifest_path <- opt$manifest

# ---- checks ----
if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

if (!file.exists(template_path)) {
  stop("Report template not found: ", template_path)
}

if (manifest_path != "" && !file.exists(manifest_path)) {
  stop("Module manifest not found: ", manifest_path)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Convert all important paths to absolute paths so the Rmd can find files
results_dir_abs <- normalizePath(results_dir, mustWork = TRUE)
template_path_abs <- normalizePath(template_path, mustWork = TRUE)
out_dir_abs <- normalizePath(out_dir, mustWork = TRUE)
manifest_path_abs <- if (manifest_path == "") "" else normalizePath(manifest_path, mustWork = TRUE)

output_file <- "local_qc_report.html"

message("Rendering report...")
message("Results directory: ", results_dir_abs)
message("Template: ", template_path_abs)
message("Output directory: ", out_dir_abs)
if (manifest_path_abs != "") message("Module manifest: ", manifest_path_abs)

rmarkdown::render(
  input = template_path_abs,
  output_file = output_file,
  output_dir = out_dir_abs,
  params = list(
    results_dir = results_dir_abs,
    report_title = report_title,
    run_manifest = manifest_path_abs
  ),
  envir = new.env(parent = globalenv())
)

message("HTML report generated:")
message(file.path(out_dir_abs, output_file))
