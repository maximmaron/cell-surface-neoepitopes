#!/bin/bash
# =============================================================================
# 05_stringtie_per_sample.sh
# SLURM array job: per-sample StringTie transcript assembly.
#
# Skip logic: ATOMIC — skips only if BOTH the GTF and strand file exist
# and are non-empty. If either is missing, both are regenerated.
#
# Strandedness:
#   - Primary call: fraction > 0.75 in one direction
#   - Fallback: ratio > 4:1 between directions (handles low-quality samples
#     where high "failed to determine" deflates absolute fractions)
#   - Genuinely unstranded libraries (~0.5/0.5) are correctly left unstranded
#
# Manifest columns used:
#   sample_id(1), cohort(2), tier(3), disease(4), seq_type(5), fq1(6), fq2(7)
#   BAM is derived from dirname(fq1)/${sample_id}_Aligned.sortedByCoord.out.bam
#
# Usage:
#   bash 05_stringtie_per_sample.sh   ← prints submission command
#   sbatch --array=1-N%20 05_stringtie_per_sample.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=stringtie
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
# NOTE: --output and --error set in submission command to land in OUTDIR/logs/

# ── Paths ─────────────────────────────────────────────────────────────────────
MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"
BED_FILE="/data1/abdelwao/maxim/annotations/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed12"
OUTDIR="/data1/abdelwao/maxim/splicing_pipeline/stringtie"
N_READS=200000

mkdir -p "${OUTDIR}/logs" "${OUTDIR}/gtf" "${OUTDIR}/abundance" "${OUTDIR}/strand"

# ── Setup mode (no SLURM_ARRAY_TASK_ID) ──────────────────────────────────────
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    N=$(( $(wc -l < "${MANIFEST}") - 1 ))
    echo "Manifest : ${MANIFEST}"
    echo "Samples  : ${N}"
    echo "Output   : ${OUTDIR}"
    echo ""
    echo "Submit with:"
    echo "  sbatch \\"
    echo "    --array=1-${N}%20 \\"
    echo "    --output=${OUTDIR}/logs/stringtie_%A_%a.out \\"
    echo "    --error=${OUTDIR}/logs/stringtie_%A_%a.err \\"
    echo "    $(realpath $0)"
    echo ""
    echo "To test a single sample:"
    echo "  sbatch \\"
    echo "    --array=1 \\"
    echo "    --output=${OUTDIR}/logs/stringtie_%A_%a.out \\"
    echo "    --error=${OUTDIR}/logs/stringtie_%A_%a.err \\"
    echo "    $(realpath $0)"
    exit 0
fi

# ── Load tools ────────────────────────────────────────────────────────────────
source /admin/software/anaconda/1.11.1/bin/activate
conda activate gatk

# ── Read sample from manifest ─────────────────────────────────────────────────
LINE=$(awk -F'\t' -v task="${SLURM_ARRAY_TASK_ID}" 'NR==task+1' "${MANIFEST}")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: No entry for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

SAMPLE_ID=$(echo "${LINE}" | cut -f1)
COHORT=$(echo "${LINE}"    | cut -f2)
TIER=$(echo "${LINE}"      | cut -f3)
FQ1=$(echo "${LINE}"       | cut -f6)

BAM="$(dirname ${FQ1})/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"
OUT_GTF="${OUTDIR}/gtf/${SAMPLE_ID}.gtf"
OUT_ABUND="${OUTDIR}/abundance/${SAMPLE_ID}.abundance.tsv"
STRAND_FILE="${OUTDIR}/strand/${SAMPLE_ID}.strand.txt"
STRAND_RAW="${OUTDIR}/strand/${SAMPLE_ID}_infer_experiment.txt"

echo "======================================="
echo "Job    : ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Sample : ${SAMPLE_ID}"
echo "Cohort : ${COHORT} (Tier ${TIER})"
echo "BAM    : ${BAM}"
echo "Started: $(date)"
echo "======================================="

# ── Validate BAM ──────────────────────────────────────────────────────────────
if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}" >&2
    exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
    echo "ERROR: BAM index not found: ${BAM}.bai" >&2
    exit 1
fi

# ── Skip if already complete (ATOMIC: both GTF and strand file must exist) ────
if [[ -f "${OUT_GTF}"     && -s "${OUT_GTF}" ]] && \
   [[ -f "${STRAND_FILE}" && -s "${STRAND_FILE}" ]]; then
    echo "Skipping ${SAMPLE_ID} — GTF and strand file already exist"
    exit 0
fi

# ── Detect strandedness with RSeQC infer_experiment.py ───────────────────────
echo ""
echo "── Detecting strandedness (RSeQC) ──────────────────────"

infer_experiment.py \
    -r "${BED_FILE}" \
    -i "${BAM}" \
    -s "${N_READS}" \
    > "${STRAND_RAW}" 2>/dev/null

echo "  RSeQC output:"
cat "${STRAND_RAW}" | sed 's/^/    /'

# Parse fractions — use explicit -z guard (grep -oP exits 0 on no match
# in some versions, making || echo "0" unreliable)
FORWARD_PCT=$(grep "1++,1--,2+-,2-+" "${STRAND_RAW}" | grep -oP '(?<=: )[0-9.]+')
REVERSE_PCT=$(grep "1+-,1-+,2++,2--" "${STRAND_RAW}" | grep -oP '(?<=: )[0-9.]+')
[[ -z "${FORWARD_PCT}" ]] && FORWARD_PCT="0"
[[ -z "${REVERSE_PCT}" ]] && REVERSE_PCT="0"

# Primary call: absolute threshold > 0.75
FWD_CALL=$(echo "${FORWARD_PCT} > 0.75" | bc -l 2>/dev/null || echo 0)
REV_CALL=$(echo "${REVERSE_PCT} > 0.75" | bc -l 2>/dev/null || echo 0)

# Ratio-based fallback: if neither exceeds 0.75, check if one direction
# dominates >4:1. Handles samples with high "failed to determine" fraction
# (e.g. low-quality BEAT samples) where absolute values are deflated but
# the direction is unambiguous.
if [[ "${FWD_CALL}" -eq 0 && "${REV_CALL}" -eq 0 ]]; then
    if (( $(echo "${FORWARD_PCT} > 0 && ${REVERSE_PCT} > 0" | bc -l) )); then
        FWD_CALL=$(echo "${FORWARD_PCT} / ${REVERSE_PCT} > 4" | bc -l 2>/dev/null || echo 0)
        REV_CALL=$(echo "${REVERSE_PCT} / ${FORWARD_PCT} > 4" | bc -l 2>/dev/null || echo 0)
        [[ "${FWD_CALL}" -eq 1 ]] && echo "  Ratio fallback: forward dominant (${FORWARD_PCT}/${REVERSE_PCT})"
        [[ "${REV_CALL}" -eq 1 ]] && echo "  Ratio fallback: reverse dominant (${REVERSE_PCT}/${FORWARD_PCT})"
    fi
fi

if [[ "${FWD_CALL}" -eq 1 ]]; then
    STRAND_CALL="forward"
    STRAND_FLAG="--fr"
    RMATS_LIBTYPE="fr-secondstrand"
elif [[ "${REV_CALL}" -eq 1 ]]; then
    STRAND_CALL="reverse"
    STRAND_FLAG="--rf"
    RMATS_LIBTYPE="fr-firststrand"
else
    STRAND_CALL="unstranded"
    STRAND_FLAG=""
    RMATS_LIBTYPE="fr-unstranded"
fi

echo ""
echo "  Forward pct : ${FORWARD_PCT}"
echo "  Reverse pct : ${REVERSE_PCT}"
echo "  Call        : ${STRAND_CALL} ${STRAND_FLAG:-(no flag)}"
echo "  rMATS type  : ${RMATS_LIBTYPE}"

echo -e "${SAMPLE_ID}\t${COHORT}\t${STRAND_CALL}\t${STRAND_FLAG}\t${FORWARD_PCT}\t${REVERSE_PCT}\t${RMATS_LIBTYPE}" \
    > "${STRAND_FILE}"

# ── Switch to StringTie conda env ─────────────────────────────────────────────
conda deactivate
conda activate stringtie

# ── Run StringTie ─────────────────────────────────────────────────────────────
echo ""
echo "── StringTie assembly ──────────────────────────────────"

stringtie \
    "${BAM}" \
    -G "${GTF}" \
    -o "${OUT_GTF}" \
    -A "${OUT_ABUND}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    ${STRAND_FLAG} \
    -c 1 \
    -f 0.01 \
    -m 200 \
    -j 1 \
    -a 10 \
    -l "${SAMPLE_ID}"

exit_code=$?

if [[ ${exit_code} -ne 0 ]]; then
    echo "ERROR: StringTie failed for ${SAMPLE_ID}" >&2
    # Remove partial outputs so the atomic skip check doesn't leave stale files
    rm -f "${OUT_GTF}" "${OUT_ABUND}" "${STRAND_FILE}"
    exit 1
fi

# Verify GTF is non-empty
if [[ ! -s "${OUT_GTF}" ]]; then
    echo "ERROR: StringTie produced empty GTF for ${SAMPLE_ID}" >&2
    rm -f "${OUT_GTF}" "${OUT_ABUND}" "${STRAND_FILE}"
    exit 1
fi

# ── Summary ───────────────────────────────────────────────────────────────────
N_TRANSCRIPTS=$(awk '$3=="transcript"' "${OUT_GTF}" | wc -l)
N_NOVEL=$(awk '$3=="transcript"' "${OUT_GTF}" | grep -vc 'reference_id' || true)

echo ""
echo "── Complete: ${SAMPLE_ID} ──────────────────────────────"
echo "  Transcripts : ${N_TRANSCRIPTS}"
echo "  Novel       : ${N_NOVEL}"
echo "  Strandedness: ${STRAND_CALL}"
echo "  Finished    : $(date)"

conda deactivate