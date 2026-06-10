#!/bin/bash

#SBATCH -J dropoff
#SBATCH -D /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool
#SBATCH -t 24:00:00
#SBATCH --mem=175G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/11_dropoff_%j.out
#SBATCH -e /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/11_dropoff_%j.err

set -euo pipefail

echo "============================================================"
echo "Running step: 11_dropoff"
echo "Module: /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Dropoff.sh"
echo "Config: /kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/configs/cohortA_10samples_config.sh"
echo "Extra args: "
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

bash "/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Dropoff.sh" "/kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/Cohort_A/configs/cohortA_10samples_config.sh" 

echo "============================================================"
echo "Finished step: 11_dropoff"
echo "Finished: $(date)"
echo "============================================================"
