# RStudio launcher for the local cfRNA QC workflow.
#
# 1. Copy Local_tool/config/local_qc_config.example.R to a personal config file.
# 2. Edit the two paths below.
# 3. Open this file in RStudio and click Source.

# Path to this repository's Local_tool directory.
local_tool_dir <- "/path/to/cfRNA_QC_tool/Local_tool"

# Path to your edited copy of local_qc_config.example.R.
config_path <- file.path(local_tool_dir, "config", "local_qc_config.R")

# Do not edit below this line.
local_tool_dir <- normalizePath(local_tool_dir, mustWork = FALSE)
config_path <- normalizePath(config_path, mustWork = FALSE)
runner_path <- file.path(local_tool_dir, "scripts", "run_all_local_qc.R")
rscript_path <- file.path(R.home("bin"), "Rscript")

if (!dir.exists(local_tool_dir)) {
  stop("Local_tool directory not found: ", local_tool_dir)
}

if (!file.exists(config_path)) {
  stop(
    "Config file not found: ", config_path,
    "\nCopy Local_tool/config/local_qc_config.example.R and update config_path."
  )
}

if (!file.exists(runner_path)) {
  stop("Master runner not found: ", runner_path)
}

if (!file.exists(rscript_path)) {
  stop("Rscript executable not found for the current R installation: ", rscript_path)
}

message("Launching local QC workflow...")
message("Config: ", config_path)

exit_status <- system2(
  command = rscript_path,
  args = c(shQuote(runner_path), "--config", shQuote(config_path)),
  stdout = "",
  stderr = ""
)

if (!identical(exit_status, 0L)) {
  stop("The local QC workflow exited with status ", exit_status, ". See the RStudio console and run logs for details.")
}

message("Local QC workflow completed successfully.")
