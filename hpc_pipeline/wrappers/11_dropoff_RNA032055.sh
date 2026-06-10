#!/bin/bash

#SBATCH -J dropoff_RNA032055
#SBATCH -D /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool
#SBATCH -t 12:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/11_dropoff_RNA032055_%j.out
#SBATCH -e /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/11_dropoff_RNA032055_%j.err

set -euo pipefail

echo "============================================================"
echo "Running step: 11_dropoff_RNA032055"
echo "Module: /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Dropoff.sh"
echo "Config: /kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/configs/cohortA_10samples_config.sh"
echo "Extra args: RNA032055"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

bash "/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Dropoff.sh" "/kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/configs/cohortA_10samples_config.sh" RNA032055

echo "============================================================"
echo "Finished step: 11_dropoff_RNA032055"
echo "Finished: $(date)"
echo "============================================================"
