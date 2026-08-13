#!/bin/bash

# ============================================================
# Cohort A 10-sample QC config
# ============================================================

# Optional cluster setup for the submission host. Leave these empty when the
# desired Slurm cluster has already been selected before submitting the workflow.
CLUSTER_MODULE="cluster/doduo"
CLUSTER_ENV_MODULE="env/software/doduo"

# Input manifest
SAMPLESHEET="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master/inputs/cohortA_10samples_samplesheet.tsv"

# Reference files
GTF="/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf"
EXON_BED="/data/gent/vo/000/gvo00027/resources/Ensembl_bedregions/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S_exons_sorted_merged.bed"

# Output directory
OUTDIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master"
