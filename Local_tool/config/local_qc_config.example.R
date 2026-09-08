# Copy this file outside the repository, update the paths, then run:
# Rscript Local_tool/scripts/run_all_local_qc.R --config /path/to/local_qc_config.R

# Required user paths
COUNTS <- "/Users/maherhabboub/Desktop/University/Thesis/Project/rnaqc_tool/data/silverseq_htseq_counts_combined.tsv"
OUT_ROOT <- "/Users/maherhabboub/Desktop/University/Thesis/Project/rnaqc_tool/results/Silver_seq/"

# Optional user paths. Leave empty to omit the associated input.
METADATA <- ""
HPC_SUMMARY <- ""

# Module selection. Normalization always runs and deliberately has no switch.
BIOTYPE_ENABLED <- TRUE
PCA_SSGSEA_ENABLED <- TRUE
SEX_INFERENCE_ENABLED <- TRUE

# Existing runner settings
TOP_N <- 1000L
MAX_PCS <- 5L
TOP_SCATTER <- 10L
REPORT_TITLE <- "Local RNA-seq QC Report"
STOP_ON_FAIL <- FALSE
