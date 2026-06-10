#!/bin/bash

#SBATCH -J fastqc
#SBATCH -D /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/03_fastqc_%j.out
#SBATCH -e /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/03_fastqc_%j.err

set -euo pipefail

echo "============================================================"
echo "Running step: 03_fastqc"
echo "Module: /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Fastqc.sh"
echo "Config: /kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh"
echo "Extra args: "
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

bash "/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Fastqc.sh" "/kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh" 

echo "============================================================"
echo "Finished step: 03_fastqc"
echo "Finished: $(date)"
echo "============================================================"
