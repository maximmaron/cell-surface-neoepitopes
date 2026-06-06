#!/bin/bash
# =============================================================================
# 01_star_pass1.sh
# SLURM array job: fastp trimming + STAR pass 1 alignment
#
# Goal: generate SJ.out.tab junction files for pooled two-pass.
# No BAM output — pass 2 produces the final BAMs.
#
# Index selection is automatic based on detected read length:
#   ≤ 50bp → gencode_v49_rl50   (TCGA, Bodymap2)
#   > 50bp → gencode_v49_rl100  (all other cohorts)
#
# Outputs per sample (written alongside FASTQs):
#   {SAMPLE}_SJ.out.tab       ← junction file for pooled two-pass
#   {SAMPLE}_Log.final.out    ← alignment QC stats
#
# fastp QC reports:
#   /data1/abdelwao/maxim/splicing_pipeline/qc/fastp/
#
# Usage:
#   bash 01_star_pass1.sh          ← prints submission command
#   sbatch --array=1-N%20 01_star_pass1.sh  ← run jobs
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=star_pass1
#SBATCH --output=logs/pass1_%A_%a.out
#SBATCH --error=logs/pass1_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16

# ── Containers ────────────────────────────────────────────────────────────────
STAR_SIF="/data1/abdelwao/maxim/containers/star_2.7.11b.sif"
FASTP_SIF="/data1/abdelwao/maxim/containers/fastp_latest.sif"

# ── Reference ─────────────────────────────────────────────────────────────────
INDEX_BASE="/data1/abdelwao/maxim/annotations/Homo_sapiens/STAR"
GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"

# ── Manifest + QC ─────────────────────────────────────────────────────────────
MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
QC_DIR="/data1/abdelwao/maxim/splicing_pipeline/qc/fastp"
SCRATCH_BASE="/scratch/abdelwao"

mkdir -p logs "${QC_DIR}"

# ── Setup mode: print submission command if not running under SLURM ───────────
if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    N=$(( $(wc -l < "${MANIFEST}") - 1 ))
    echo "Manifest : ${MANIFEST}"
    echo "Samples  : ${N}"
    echo ""
    echo "Submit with:"
    echo "  sbatch --array=1-${N}%10 $(basename $0)"
    echo ""
    echo "To test a single sample first:"
    echo "  sbatch --array=1 $(basename $0)"
    exit 0
fi

# ── Read sample from manifest ─────────────────────────────────────────────────
LINE=$(awk -F'\t' -v task="${SLURM_ARRAY_TASK_ID}" 'NR==task+1' "${MANIFEST}")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: No entry for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

SAMPLE_ID=$(echo "${LINE}" | awk -F'\t' '{print $1}')
COHORT=$(echo "${LINE}"    | awk -F'\t' '{print $2}')
TIER=$(echo "${LINE}"      | awk -F'\t' '{print $3}')
SEQ_TYPE=$(echo "${LINE}"  | awk -F'\t' '{print $5}')
FQ1=$(echo "${LINE}"       | awk -F'\t' '{print $6}')
FQ2=$(echo "${LINE}"       | awk -F'\t' '{print $7}')

echo "======================================="
echo "Job   : ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE_ID}"
echo "Cohort: ${COHORT} (Tier ${TIER})"
echo "Type  : ${SEQ_TYPE}"
echo "FQ1   : ${FQ1}"
[[ "${SEQ_TYPE}" == "paired" ]] && echo "FQ2   : ${FQ2}"
echo "======================================="

# ── Validate inputs ───────────────────────────────────────────────────────────
[[ ! -f "${FQ1}" ]] && { echo "ERROR: FQ1 not found: ${FQ1}" >&2; exit 1; }
[[ "${SEQ_TYPE}" == "paired" && ! -f "${FQ2}" ]] && \
    { echo "ERROR: FQ2 not found: ${FQ2}" >&2; exit 1; }

# ── Output directory = same directory as FASTQs ───────────────────────────────
OUT_DIR=$(dirname "${FQ1}")
OUTPUT_PREFIX="${OUT_DIR}/${SAMPLE_ID}_"

# ── Skip if already completed ─────────────────────────────────────────────────
if [[ -f "${OUT_DIR}/${SAMPLE_ID}_SJ.out.tab" && \
      -f "${OUT_DIR}/${SAMPLE_ID}_Log.final.out" ]]; then
    echo "Skipping ${SAMPLE_ID} — pass 1 already complete"
    exit 0
fi

# ── Scratch setup ─────────────────────────────────────────────────────────────
SCRATCH="${SCRATCH_BASE}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${SCRATCH}"
trap "echo 'Cleaning up scratch: ${SCRATCH}'; rm -rf '${SCRATCH}'" EXIT

# ── Detect read length → select index ────────────────────────────────────────
READ_LENGTH=$(zcat "${FQ1}" | \
    awk 'NR%4==2{print length($0); count++; if(count>=100) exit}' | \
    sort -n | tail -1)

if [[ "${READ_LENGTH}" -le 50 ]]; then
    INDEX_DIR="${INDEX_BASE}/gencode_v49_rl50"
else
    INDEX_DIR="${INDEX_BASE}/gencode_v49_rl100"
fi

echo "Read length : ${READ_LENGTH}bp → index: $(basename ${INDEX_DIR})"

# ── Validate index exists ─────────────────────────────────────────────────────
if [[ ! -f "${INDEX_DIR}/SA" ]]; then
    echo "ERROR: STAR index not found: ${INDEX_DIR}" >&2
    echo "Run: sbatch --export=READ_LENGTH=${READ_LENGTH} build_star_index.sh" >&2
    exit 1
fi

# ── fastp trimming ────────────────────────────────────────────────────────────
TRIMMED_R1="${SCRATCH}/${SAMPLE_ID}_R1.fastq"
TRIMMED_R2="${SCRATCH}/${SAMPLE_ID}_R2.fastq"

echo ""
echo "── fastp ───────────────────────────────────────────────"

if [[ "${SEQ_TYPE}" == "paired" ]]; then
    singularity exec \
        --bind /data1/abdelwao:/data1/abdelwao,${SCRATCH_BASE}:${SCRATCH_BASE} \
        "${FASTP_SIF}" fastp \
        --in1 "${FQ1}" \
        --in2 "${FQ2}" \
        --out1 "${TRIMMED_R1}" \
        --out2 "${TRIMMED_R2}" \
        --thread "${SLURM_CPUS_PER_TASK}" \
        --qualified_quality_phred 20 \
        --length_required 25 \
        --detect_adapter_for_pe \
        --json "${QC_DIR}/${SAMPLE_ID}_fastp.json" \
        --html "${QC_DIR}/${SAMPLE_ID}_fastp.html"
else
    singularity exec \
        --bind /data1/abdelwao:/data1/abdelwao,${SCRATCH_BASE}:${SCRATCH_BASE} \
        "${FASTP_SIF}" fastp \
        --in1 "${FQ1}" \
        --out1 "${TRIMMED_R1}" \
        --thread "${SLURM_CPUS_PER_TASK}" \
        --qualified_quality_phred 20 \
        --length_required 25 \
        --json "${QC_DIR}/${SAMPLE_ID}_fastp.json" \
        --html "${QC_DIR}/${SAMPLE_ID}_fastp.html"
fi

[[ $? -ne 0 ]] && { echo "ERROR: fastp failed for ${SAMPLE_ID}" >&2; exit 1; }

# ── STAR pass 1 — junction discovery only, no BAM ────────────────────────────
echo ""
echo "── STAR pass 1 (junction discovery) ───────────────────"

[[ "${SEQ_TYPE}" == "paired" ]] && \
    READ_FILES="${TRIMMED_R1} ${TRIMMED_R2}" || \
    READ_FILES="${TRIMMED_R1}"

singularity exec \
    --bind /data1/abdelwao:/data1/abdelwao,${SCRATCH_BASE}:${SCRATCH_BASE} \
    "${STAR_SIF}" STAR \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "${INDEX_DIR}" \
    --readFilesIn ${READ_FILES} \
    --outFileNamePrefix "${OUTPUT_PREFIX}" \
    --outTmpDir "${SCRATCH}/STAR_tmp" \
    --sjdbGTFfile "${GTF}" \
    --twopassMode None \
    --outSAMtype None \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs None \
    --alignSJoverhangMin 6 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterMismatchNmax 999 


[[ $? -ne 0 ]] && { echo "ERROR: STAR pass 1 failed for ${SAMPLE_ID}" >&2; exit 1; }

# ── Verify outputs ────────────────────────────────────────────────────────────
SJ="${OUT_DIR}/${SAMPLE_ID}_SJ.out.tab"
[[ ! -f "${SJ}" ]] && \
    { echo "ERROR: SJ.out.tab not found" >&2; exit 1; }

SJ_COUNT=$(wc -l < "${SJ}")
[[ "${SJ_COUNT}" -eq 0 ]] && \
    { echo "ERROR: SJ.out.tab is empty" >&2; exit 1; }

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "── Complete ────────────────────────────────────────────"
echo "SJ.out.tab : ${SJ} (${SJ_COUNT} junctions)"
echo "Log        : ${OUT_DIR}/${SAMPLE_ID}_Log.final.out"
echo "fastp QC   : ${QC_DIR}/${SAMPLE_ID}_fastp.json"
echo ""
grep -E "Uniquely mapped reads %|Number of splices: Total|GT/AG|AT/AC|Non-canonical" \
    "${OUT_DIR}/${SAMPLE_ID}_Log.final.out"
echo ""
echo "Completed: ${SAMPLE_ID}"