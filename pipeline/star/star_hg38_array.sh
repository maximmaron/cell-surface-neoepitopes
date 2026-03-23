#!/bin/bash
# =============================================================================
# 00_prep_bam_lists.sh
# Convert comma-separated BAM list files (as used by rMATS) to one-per-line
# format required by this LeafCutter pipeline, and compute SLURM array bounds.
#
# Run this ONCE from your working directory before submitting SLURM jobs.
# Usage:  bash 00_prep_bam_lists.sh
# =============================================================================

WORK_DIR="/data1/abdelwao/shared/splicing_analysis/leafcutter/SRSF2_vs_normal"
LIST_DIR="${WORK_DIR}/bam_lists"
mkdir -p "${LIST_DIR}"

# Source: your existing rMATS-format comma-delimited lists
SRSF2_CSV="/data1/abdelwao/shared/splicing_analysis/leafcutter/SRSF2_vs_normal/srsf2_file_list.txt"
NORMAL_CSV="/data1/abdelwao/shared/splicing_analysis/leafcutter/SRSF2_vs_normal/normal_control_file_list.txt"

# ── Convert to one-per-line ───────────────────────────────────────────────────
tr ',' '\n' < "${SRSF2_CSV}"  | sed '/^$/d' > "${LIST_DIR}/srsf2_bams.txt"
tr ',' '\n' < "${NORMAL_CSV}" | sed '/^$/d' > "${LIST_DIR}/normal_bams.txt"

N_SRSF2=$(wc -l < "${LIST_DIR}/srsf2_bams.txt")
N_NORMAL=$(wc -l < "${LIST_DIR}/normal_bams.txt")

echo "SRSF2  BAMs: ${N_SRSF2}  -> ${LIST_DIR}/srsf2_bams.txt"
echo "Normal BAMs: ${N_NORMAL} -> ${LIST_DIR}/normal_bams.txt"

echo ""
echo "Submit jobs with:"
echo "  sbatch --export=COHORT=srsf2  --array=1-${N_SRSF2}  01_bam_to_junc.sh"
echo "  sbatch --export=COHORT=normal --array=1-${N_NORMAL} 01_bam_to_junc.sh"