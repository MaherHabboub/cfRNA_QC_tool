#!/bin/bash

#SBATCH -J splice_junction
#SBATCH -D /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/09_splice_junction_%j.out
#SBATCH -e /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/09_splice_junction_%j.err

set -euo pipefail

echo "============================================================"
echo "Running step: 09_splice_junction"
echo "Module: /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Splice_Junction.sh"
echo "Config: /kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh"
echo "Extra args: "
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

bash "/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Splice_Junction.sh" "/kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh" 

echo "============================================================"
echo "Finished step: 09_splice_junction"
echo "Finished: $(date)"
echo "============================================================"
