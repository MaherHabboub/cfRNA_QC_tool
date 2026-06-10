#!/bin/bash

#SBATCH -J make_bins_only
#SBATCH -D /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool
#SBATCH -t 03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/02_make_bins_only_%j.out
#SBATCH -e /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/02_make_bins_only_%j.err

set -euo pipefail

echo "============================================================"
echo "Running Make_Dropoff_Bins.sh only"
echo "Config: /kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

bash "/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Make_Dropoff_Bins.sh" "/kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh"

echo "============================================================"
echo "Finished Make_Dropoff_Bins.sh"
echo "Finished: $(date)"
echo "============================================================"
