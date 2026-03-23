#!/bin/bash
# =============================================================================
# 02_cluster_introns.sh
# Non-array SLURM job: collect all junc files and run leafcutter-cluster
#
# Depends on: 01_bam_to_junc.sh completing successfully for ALL samples.
# Usage:  sbatch 02_cluster_introns.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=lc_cluster
#SBATCH --output=logs/lc_cluster_%j.out
#SBATCH --error=logs/lc_cluster_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

source /admin/software/anaconda/1.11.1/bin/activate
conda activate leafcutter

WORK_DIR="/scratch/abdelwao/leafcutter_srsf2"
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
# -m  : minimum reads per cluster (50 is appropriate for 191 samples)
# -l  : maximum intron length
# -o  : output prefix
# -r  : run directory (output goes here)
# --checkchrom: remove non-standard chroms (important for hg38 with alt contigs)
leafcutter-cluster \
    -j "${JUNC_LIST}" \
    -m 50 \
    -l 500000 \
    -o "srsf2_vs_normal" \
    -r "${CLUSTER_DIR}" 
    
echo ""
echo "Clustering complete."
echo "Output: ${CLUSTER_DIR}/srsf2_vs_normal_perind_numers.counts.gz"
echo ""
echo "Next: build groups file and run 03_differential_splicing.sh"