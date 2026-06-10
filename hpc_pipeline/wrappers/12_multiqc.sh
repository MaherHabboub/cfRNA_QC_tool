#!/bin/bash

#SBATCH -J multiqc
#SBATCH -D /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/12_multiqc_%j.out
#SBATCH -e /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/logs/12_multiqc_%j.err

set -euo pipefail

echo "============================================================"
echo "Running step: 12_multiqc"
echo "Module: /scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Multiqc.sh"
echo "Config: /kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh"
echo "Extra args: "
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

bash "/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/modules/Multiqc.sh" "/kyukon/scratch/gent/vo/000/gvo00027/projects/MHB/HPC_QC_Tool/configs/silverseq_config.sh" 

echo "============================================================"
echo "Finished step: 12_multiqc"
echo "Finished: $(date)"
echo "============================================================"
