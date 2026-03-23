#!/bin/bash
# =============================================================================
# 01_bam_to_junc.sh
# SLURM array job: convert BAM files -> junction files using regtools
# Per-sample strandedness is looked up from bam_summary.tsv (infer_experiment
# results) rather than hardcoded, since cohorts differ:
#   BEAT AML / TCGA-AB : reverse  -> regtools -s RF
#   HPA / Bodymap2 / Maiga / Madan / Leucegene : unstranded -> regtools -s XS
#
# Usage:
#   sbatch --export=COHORT=srsf2  --array=1-N  01_bam_to_junc.sh
#   sbatch --export=COHORT=normal --array=1-N  01_bam_to_junc.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=lc_bam2junc
#SBATCH --output=logs/bam2junc_%A_%a.out
#SBATCH --error=logs/bam2junc_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

source /admin/software/anaconda/1.11.1/bin/activate
conda activate leafcutter

# ── Paths ─────────────────────────────────────────────────────────────────────
WORK_DIR="/scratch/abdelwao/leafcutter_srsf2"
JUNC_DIR="${WORK_DIR}/junc_files"
mkdir -p "${JUNC_DIR}" logs

# infer_experiment.py results — used to determine per-sample regtools -s flag
BAM_SUMMARY="/data1/abdelwao/shared/splicing_analysis/rmats/bam_qc_results/bam_summary.tsv"

# ── Select cohort ─────────────────────────────────────────────────────────────
COHORT="${COHORT:-srsf2}"

if [[ "${COHORT}" == "srsf2" ]]; then
    BAM_LIST_FILE="/data1/abdelwao/shared/splicing_analysis/leafcutter/SRSF2_vs_normal/bam_lists/srsf2_bams.txt"
elif [[ "${COHORT}" == "normal" ]]; then
    BAM_LIST_FILE="/data1/abdelwao/shared/splicing_analysis/leafcutter/SRSF2_vs_normal/bam_lists/normal_bams.txt"
else
    echo "Unknown COHORT=${COHORT}. Use srsf2 or normal." >&2
    exit 1
fi

# ── Select BAM for this array task ───────────────────────────────────────────
mapfile -t BAMS < "${BAM_LIST_FILE}"
NBAMS=${#BAMS[@]}
echo "Cohort: ${COHORT} | Total BAMs: ${NBAMS}"

IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))
BAM="${BAMS[$IDX]}"

if [[ -z "${BAM}" ]]; then
    echo "No BAM at index ${IDX}" >&2
    exit 1
fi

SAMPLE=$(basename "${BAM}" _Aligned.sortedByCoord.out.bam)
JUNC="${JUNC_DIR}/${COHORT}_${SAMPLE}.junc"

echo "Processing : ${SAMPLE}"
echo "  BAM      : ${BAM}"
echo "  JUNC     : ${JUNC}"

# ── Look up strandedness from bam_summary.tsv ─────────────────────────────────
# Column 12 is stranded_call: "unstranded", "forward", or "reverse"
# awk matches on bam_path column (col 2) then prints stranded_call (col 12)
STRAND_CALL=$(awk -F'\t' -v bam="${BAM}" 'NR>1 && $2==bam {print $12}' "${BAM_SUMMARY}")

if [[ -z "${STRAND_CALL}" ]]; then
    echo "WARNING: ${BAM} not found in bam_summary.tsv — defaulting to XS (unstranded)" >&2
    STRAND_CALL="unstranded"
fi

# Map infer_experiment stranded_call -> regtools -s flag (regtools 1.0.0 syntax)
case "${STRAND_CALL}" in
    unstranded) REGTOOLS_STRAND="XS" ;;
    forward)    REGTOOLS_STRAND="FR" ;;
    reverse)    REGTOOLS_STRAND="RF" ;;
    *)
        echo "WARNING: unrecognized stranded_call '${STRAND_CALL}' — defaulting to XS" >&2
        REGTOOLS_STRAND="XS"
        ;;
esac

echo "  Strand call : ${STRAND_CALL} -> regtools -s ${REGTOOLS_STRAND}"

# ── Index BAM if needed ───────────────────────────────────────────────────────
if [[ ! -f "${BAM}.bai" ]]; then
    echo "Indexing ${BAM}..."
    samtools index -@ "${SLURM_CPUS_PER_TASK}" "${BAM}"
fi

# ── Extract junctions (regtools 1.0.0 syntax) ─────────────────────────────────
# -a 8      : minimum anchor length
# -m 50     : minimum intron length
# -M 500000 : maximum intron length
# -s        : strandedness (XS/RF/FR — looked up per sample above)
regtools junctions extract \
    -a 8 \
    -m 50 \
    -M 500000 \
    -s "${REGTOOLS_STRAND}" \
    "${BAM}" \
    -o "${JUNC}"

if [[ $? -eq 0 ]]; then
    echo "Done: ${JUNC}"
else
    echo "regtools failed for ${SAMPLE}" >&2
    exit 1
fi

conda deactivate