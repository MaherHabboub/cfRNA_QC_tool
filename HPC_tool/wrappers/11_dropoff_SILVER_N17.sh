#!/bin/bash

#SBATCH -J dropoff_SILVER_N17
#SBATCH -D /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool
#SBATCH -t 12:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/11_dropoff_SILVER_N17_%j.out
#SBATCH -e /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/11_dropoff_SILVER_N17_%j.err

set -euo pipefail

echo "============================================================"
echo "Running step: 11_dropoff_SILVER_N17"
echo "Module: /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Dropoff.sh"
echo "Config: /kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh"
echo "Extra args: SILVER_N17"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

bash "/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Dropoff.sh" "/kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh" SILVER_N17

echo "============================================================"
echo "Finished step: 11_dropoff_SILVER_N17"
echo "Finished: $(date)"
echo "============================================================"
