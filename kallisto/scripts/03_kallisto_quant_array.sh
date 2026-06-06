#!/bin/bash
#SBATCH --job-name=kallisto_quant
#SBATCH --output=logs/03_kallisto_quant_%A_%a.log
#SBATCH --array=1-1166%40
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# ==============================================================================
# 03_kallisto_quant_array.sh
#
# Runs kallisto quant directly from FASTQ files.
#
# Manifest columns (tab-separated, header row):
#   sample_id | fastq_r1 | fastq_r2 | strandedness | sf_group | cohort | condition
#
# For single-end samples fastq_r2 should be "NA".
#
# Strandedness → kallisto flags:
#   fr-firststrand  → --rf-stranded   (dUTP / TruSeq stranded)
#   fr-secondstrand → --fr-stranded   (ligation-based)
#   unstranded      → (no flag)
#
# Dependencies: kallisto
# ==============================================================================


set -euo pipefail

PIPELINE_DIR="/data1/abdelwao/maxim/splicing_pipeline/kallisto"
MANIFEST="${PIPELINE_DIR}/reference/sample_manifest.tsv"
INDEX="${PIPELINE_DIR}/reference/surface_kallisto.idx"
OUT_DIR="${PIPELINE_DIR}/kallisto_quant"
THREADS=16

mkdir -p "${OUT_DIR}"

source /admin/software/anaconda/1.11.1/bin/activate
conda activate kallisto

# ── Parse manifest row ────────────────────────────────────────────────────────
LINE=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n+1' "${MANIFEST}")
SAMPLE_ID=$(echo "${LINE}" | cut -f1)
FQ1=$(echo       "${LINE}" | cut -f2)
FQ2=$(echo       "${LINE}" | cut -f3)
STRAND=$(echo    "${LINE}" | cut -f4)

echo "[$(date)] Task:     ${SLURM_ARRAY_TASK_ID}"
echo "[$(date)] Sample:   ${SAMPLE_ID}"
echo "[$(date)] FASTQ R1: ${FQ1}"
echo "[$(date)] FASTQ R2: ${FQ2}"
echo "[$(date)] Strand:   ${STRAND}"

SAMPLE_OUT="${OUT_DIR}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_OUT}"

# Skip if already complete
if [[ -f "${SAMPLE_OUT}/abundance.h5" ]]; then
    echo "[$(date)] abundance.h5 exists — skipping."
    exit 0
fi

# Validate R1
if [[ ! -f "${FQ1}" ]]; then
    echo "ERROR: FASTQ R1 not found: ${FQ1}" >&2
    exit 1
fi

# ── Strandedness flag ─────────────────────────────────────────────────────────
STRAND_FLAG=""
case "${STRAND}" in
    fr-firststrand)  STRAND_FLAG="--rf-stranded" ;;
    fr-secondstrand) STRAND_FLAG="--fr-stranded"  ;;
    unstranded)      STRAND_FLAG=""               ;;
    *)               STRAND_FLAG=""               ;;
esac

# ── Paired-end vs single-end ──────────────────────────────────────────────────
if [[ "${FQ2}" != "NA" ]] && [[ -f "${FQ2}" ]]; then
    echo "[$(date)] Mode: paired-end"

    kallisto quant \
        -i "${INDEX}" \
        -o "${SAMPLE_OUT}" \
        -t "${THREADS}" \
        --bootstrap-samples 100 \
        ${STRAND_FLAG} \
        "${FQ1}" "${FQ2}"

else
    echo "[$(date)] Mode: single-end"
    # Use read_length_median from bam_summary if present in manifest (col 8)
    # otherwise fall back to 200
    READ_LEN=$(echo "${LINE}" | cut -f8)
    if [[ "${READ_LEN}" =~ ^[0-9]+$ ]]; then
        FRAG_LEN="${READ_LEN}"
    else
        FRAG_LEN=200
    fi
    FRAG_SD=30

    kallisto quant \
        -i "${INDEX}" \
        -o "${SAMPLE_OUT}" \
        -t "${THREADS}" \
        --bootstrap-samples 100 \
        --single \
        -l "${FRAG_LEN}" \
        -s "${FRAG_SD}" \
        ${STRAND_FLAG} \
        "${FQ1}"
fi

echo "[$(date)] Done: ${SAMPLE_ID}"

# Sanity check
if [[ ! -f "${SAMPLE_OUT}/abundance.h5" ]]; then
    echo "ERROR: abundance.h5 not produced for ${SAMPLE_ID}" >&2
    exit 1
fi