#!/bin/bash
# =============================================================================
# 02_cluster_introns.sh
# Non-array SLURM job: collect all junc files and run leafcutter_cluster.py
#
# Depends on: 01_bam_to_junc.sh completing successfully for ALL samples.
# Usage:  sbatch --export=COMPARISON=SF3B1_vs_normal 02_cluster_introns.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=lc_cluster
#SBATCH --output=logs/lc_cluster_%j.out
#SBATCH --error=logs/lc_cluster_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

source /admin/software/anaconda/1.11.1/bin/activate
conda activate leafcutter

COMPARISON="${COMPARISON:?COMPARISON env var must be set (e.g. --export=COMPARISON=SF3B1_vs_normal)}"

WORK_DIR="/data1/abdelwao/maxim/splicing_pipeline/leafcutter/${COMPARISON}"
JUNC_DIR="${WORK_DIR}/junc_files"
CLUSTER_DIR="${WORK_DIR}/clustering"
mkdir -p "${CLUSTER_DIR}" logs

# ── Collect all junc files from both cohorts ──────────────────────────────────
JUNC_LIST="${WORK_DIR}/all_junc_files.txt"
ls "${JUNC_DIR}"/*.junc > "${JUNC_LIST}"
N=$(wc -l < "${JUNC_LIST}")
echo "Total junc files found: ${N}"

if [[ "${N}" -eq 0 ]]; then
    echo "ERROR: No junc files in ${JUNC_DIR}" >&2
    exit 1
fi

# ── Run intron clustering ─────────────────────────────────────────────────────
# -j  : file listing junc files (one per line)
# -m  : minimum reads per cluster (50 is standard)
# -l  : maximum intron length (match regtools -M 500000)
# -o  : output prefix
# -r  : output directory
leafcutter-cluster \
    -j "${JUNC_LIST}" \
    -m 50 \
    -l 500000 \
    -o "${COMPARISON}" \
    -r "${CLUSTER_DIR}"

echo ""
echo "Clustering complete."
echo "Output: ${CLUSTER_DIR}/${COMPARISON}_perind_numers.counts.gz"
echo ""
echo "Next: build groups file and run 03_differential_splicing.sh"

conda deactivate