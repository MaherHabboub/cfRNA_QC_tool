#!/bin/bash

# ============================================================
# Cohort A 10-sample QC config
# ============================================================

# Optional cluster setup for the submission host. Leave these empty when the
# desired Slurm cluster has already been selected before submitting the workflow.
CLUSTER_MODULE="cluster/doduo"
CLUSTER_ENV_MODULE="env/software/doduo"

# FastQC uses this many threads per sample job.
FASTQC_THREADS=2

# Downsample BAMs before duplication and gene-body QC. Set to "no" to keep
# those two modules on the full original BAMs.
DOWNSAMPLE_ENABLED="yes"
DOWNSAMPLE_TARGET_ALIGNMENTS=1000000
DOWNSAMPLE_SEED=42
DOWNSAMPLE_THREADS=4

# Core QC modules. Each runs by default; set a switch to "no" to omit that
# module from submission and from the current run's MultiQC/aggregate results.
FASTQC_ENABLED="yes"
MAPPING_ENABLED="yes"
DUPLICATION_ENABLED="yes"
INSERT_SIZE_ENABLED="yes"
GENEBODY_ENABLED="yes"
READ_DISTRIBUTION_ENABLED="yes"
SPLICE_JUNCTION_ENABLED="yes"
STRANDEDNESS_ENABLED="yes"
DROPOFF_ENABLED="yes"
# Kraken2 microbial screening of STAR-unmapped reads. Both SE and PE samples
# are supported; PE samples require STAR's Unmapped.out.mate1 and mate2 files.
KRAKEN_ENABLED="yes"

# The following retain the tested HPC Kraken2 installation and cfRNA database.
# Override them when using the tool on another HPC system.
KRAKEN_CONFIG="/scratch/gent/vo/000/gvo00027/tools/kraken/config.sh"
KRAKEN_DB="/scratch/gent/vo/000/gvo00027/tools/kraken/database/cfrna_k2_bacteria_archaea_viral_human_fungi_20260828_full"
KRAKEN_TOP_TAXA=20
KRAKEN_TOP_GENERA=15

# Input manifest. Its optional ninth column, transcriptome_bam, selects the
# transcriptome insert-size method for that sample.
SAMPLESHEET="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master/inputs/cohortA_10samples_samplesheet.tsv"

# Reference files
GTF="/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf"
EXON_BED="/data/gent/vo/000/gvo00027/resources/Ensembl_bedregions/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S_exons_sorted_merged.bed"

# Output directory
OUTDIR="/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/QC_master"
