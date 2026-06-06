#!/bin/bash
# =============================================================================
# 01_bam_to_junc.sh
# SLURM array job: convert BAM -> junction files using regtools
#
# Strandedness: all BAMs were aligned with STAR --outSAMstrandField intronMotif,
# which writes XS tags from the intron splice motif (GT-AG, GC-AG, AT-AC)
# for every spliced read, regardless of library strandedness. regtools -s XS
# therefore works correctly for all cohorts (BEAT, TCGA, Leucegene, HPA,
# Bodymap2, Maiga, Madan, Pellagatti) without per-sample lookup.
#
# Usage (COHORT matches the label used in 00_prep_bam_lists.sh):
#   sbatch --export=COHORT=cases  --array=1-N  01_bam_to_junc.sh
#   sbatch --export=COHORT=controls --array=1-N  01_bam_to_junc.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=lc_bam2junc
#SBATCH --output=logs/bam2junc_%A_%a.out
#SBATCH --error=logs/bam2junc_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

set -euo pipefail

source /admin/software/anaconda/1.11.1/bin/activate
conda activate leafcutter

# ── Paths ─────────────────────────────────────────────────────────────────────
WORK_DIR="/data1/abdelwao/maxim/splicing_pipeline/leafcutter/ZRSR2_vs_normal"
JUNC_DIR="${WORK_DIR}/junc_files"
LIST_DIR="${WORK_DIR}/bam_lists"
mkdir -p "${JUNC_DIR}" logs

# ── Resolve BAM list from COHORT ──────────────────────────────────────────────
COHORT="${COHORT:?COHORT env var must be set (e.g. --export=COHORT=cases)}"
BAM_LIST_FILE="${LIST_DIR}/${COHORT}_bams.txt"

if [[ ! -f "${BAM_LIST_FILE}" ]]; then
    echo "ERROR: BAM list not found: ${BAM_LIST_FILE}" >&2
    echo "Run 00_prep_bam_lists.sh first to generate BAM lists." >&2
    exit 1
fi

# ── Select BAM for this array task ───────────────────────────────────────────
mapfile -t BAMS < "${BAM_LIST_FILE}"
NBAMS=${#BAMS[@]}
echo "Cohort: ${COHORT} | Total BAMs: ${NBAMS}"

IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))
BAM="${BAMS[$IDX]}"

if [[ -z "${BAM}" ]]; then
    echo "ERROR: No BAM at index ${IDX} (SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID})" >&2
    exit 1
fi

if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found on disk: ${BAM}" >&2
    exit 1
fi

SAMPLE=$(basename "${BAM}" _Aligned.sortedByCoord.out.bam)
JUNC="${JUNC_DIR}/${COHORT}_${SAMPLE}.junc"

echo "Sample   : ${SAMPLE}"
echo "BAM      : ${BAM}"
echo "JUNC     : ${JUNC}"

# ── Skip if already done ──────────────────────────────────────────────────────
if [[ -f "${JUNC}" && -s "${JUNC}" ]]; then
    echo "Skipping — junction file already exists: ${JUNC}"
    exit 0
fi

# ── Index BAM if needed ───────────────────────────────────────────────────────
if [[ ! -f "${BAM}.bai" ]]; then
    echo "Indexing ${BAM}..."
    samtools index -@ "${SLURM_CPUS_PER_TASK}" "${BAM}"
fi

# ── Extract junctions ─────────────────────────────────────────────────────────
# -a 8      : minimum anchor length (bp overhang on each side of junction)
# -m 50     : minimum intron length
# -M 500000 : maximum intron length
# -s XS     : use XS tag written by STAR --outSAMstrandField intronMotif;
#             works for all cohorts regardless of original library strandedness
regtools junctions extract \
    -a 8 \
    -m 50 \
    -M 500000 \
    -s XS \
    "${BAM}" \
    -o "${JUNC}"

echo "Done: ${JUNC}"

conda deactivate