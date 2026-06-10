#!/bin/bash

# ============================================================
# Cohort A 10-sample QC config
# ============================================================

# Input files/directories
SAMPLESHEET="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master/inputs/cohortA_10samples_samplesheet.tsv"
FASTQ_DIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/fastq"
BAM_DIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/map"

# Reference files
GTF="/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf"
EXON_BED="/data/gent/vo/000/gvo00027/resources/Ensembl_bedregions/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S_exons_sorted_merged.bed"

# Output directory
OUTDIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master"
BED12="${OUTDIR}/annotation/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.bed12.bed"