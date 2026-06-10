#!/bin/bash

#SBATCH -J strandedness
#SBATCH -D /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/10_strandedness_%j.out
#SBATCH -e /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/10_strandedness_%j.err

set -euo pipefail

echo "============================================================"
echo "Running step: 10_strandedness"
echo "Module: /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Strandedness.sh"
echo "Config: /kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh"
echo "Extra args: "
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

bash "/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Strandedness.sh" "/kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh" 

echo "============================================================"
echo "Finished step: 10_strandedness"
echo "Finished: $(date)"
echo "============================================================"
