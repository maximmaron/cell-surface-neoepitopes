#!/bin/bash
#SBATCH --job-name=drimseq_dtu
#SBATCH --output=logs/05_drimseq_%A_%a.log
#SBATCH --array=1-8
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00

set -euo pipefail

PIPELINE_DIR="/data1/abdelwao/maxim/splicing_pipeline/kallisto"
mkdir -p "${PIPELINE_DIR}/logs"
cd "${PIPELINE_DIR}"

source /admin/software/anaconda/1.11.1/bin/activate
conda activate kallisto

echo "[$(date)] Starting DRIMSeq DTU — task ${SLURM_ARRAY_TASK_ID}"
echo "[$(date)] R version: $(R --version | head -1)"
echo "[$(date)] Conda env: ${CONDA_DEFAULT_ENV}"
 
Rscript scripts/05_drimseq_dtu.R "${SLURM_ARRAY_TASK_ID}"
 
echo "[$(date)] Done — task ${SLURM_ARRAY_TASK_ID}"