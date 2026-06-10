#!/bin/bash

# ============================================================
# Cohort A 10-sample QC config
# ============================================================

SAMPLESHEET="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master/inputs/cohortA_10samples_samplesheet.tsv"

FASTQ_DIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/fastq"
BAM_DIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/map"
COUNT_DIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/counts"

GTF="/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf"
OUTDIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master"

BED12="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master/annotation/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.bed12.bed"
EXON_BED="/data/gent/vo/000/gvo00027/resources/Ensembl_bedregions/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S_exons_sorted_merged.bed"

THREADS="4"

MULTIQC_MODULE="MultiQC/1.28-foss-2024a"
