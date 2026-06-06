#!/bin/bash
# =============================================================================
# 04_star_pass2.sh
# SLURM array job: fastp trimming + STAR pass 2 alignment
#
# Uses pooled pass 2 index containing all novel junctions discovered
# across all 1166 samples in pass 1. Produces final BAMs for downstream
# analysis (rMATS, LeafCutter, StringTie, neoepitope pipeline).
#
# fastp writes trimmed FASTQs to scratch, STAR reads from scratch.
# No named pipes — more reliable than streaming approach.
#
# Index selection automatic based on detected read length:
#   <= 50bp -> gencode_v49_pass2_rl50   (TCGA, Bodymap2)
#   >  50bp -> gencode_v49_pass2_rl100  (all other cohorts)
#
# Outputs per sample (written alongside FASTQs):
#   {SAMPLE}_Aligned.sortedByCoord.out.bam
#   {SAMPLE}_Aligned.sortedByCoord.out.bam.bai
#   {SAMPLE}_Aligned.toTranscriptome.out.bam   <- for StringTie/RSEM
#   {SAMPLE}_SJ.out.tab                         <- pass 2 junction file
#   {SAMPLE}_Log.final.out
#
# fastp QC reports:
#   /data1/abdelwao/maxim/splicing_pipeline/qc/fastp_pass2/
#
# Usage:
#   bash 04_star_pass2.sh          <- prints submission command
#   sbatch --array=1-N%20 04_star_pass2.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=star_pass2
#SBATCH --output=logs/pass2_%A_%a.out
#SBATCH --error=logs/pass2_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

# ── Containers ────────────────────────────────────────────────────────────────
STAR_SIF="/data1/abdelwao/maxim/containers/star_2.7.11b.sif"
FASTP_SIF="/data1/abdelwao/maxim/containers/fastp_latest.sif"
SAMTOOLS_SIF="/data1/abdelwao/maxim/containers/samtools_latest.sif"

# ── Reference ─────────────────────────────────────────────────────────────────
INDEX_BASE="/data1/abdelwao/maxim/annotations/Homo_sapiens/STAR"
GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"

# ── Manifest + QC ─────────────────────────────────────────────────────────────
MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
SCRATCH_BASE="/scratch/abdelwao/"

# ── Setup mode ────────────────────────────────────────────────────────────────
if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    N=$(( $(wc -l < "${MANIFEST}") - 1 ))
    echo "Manifest : ${MANIFEST}"
    echo "Samples  : ${N}"
    echo ""
    echo "Submit with:"
    echo "  sbatch --array=1-${N}%10 $(basename $0)"
    echo ""
    echo "To submit specific tasks only:"
    echo "  sbatch --array=TASK1,TASK2,...%10 $(basename $0)"
    echo ""
    echo "To test a single sample:"
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
if [[ -f "${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam" && \
      -f "${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam.bai" && \
      -f "${OUT_DIR}/${SAMPLE_ID}_Log.final.out" ]]; then
    echo "Skipping ${SAMPLE_ID} — pass 2 already complete"
    exit 0
fi

# ── Scratch setup ─────────────────────────────────────────────────────────────
SCRATCH="${SCRATCH_BASE}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${SCRATCH}"
trap "echo 'Cleaning up scratch: ${SCRATCH}'; rm -rf '${SCRATCH}'" EXIT

# ── Check scratch space before starting ───────────────────────────────────────
SCRATCH_FREE=$(df "${SCRATCH_BASE}" | awk 'NR==2{print $4}')
SCRATCH_FREE_GB=$(( SCRATCH_FREE / 1024 / 1024 ))
echo "Scratch free: ${SCRATCH_FREE_GB}GB"
if [[ "${SCRATCH_FREE_GB}" -lt 80 ]]; then
    echo "ERROR: Less than 80GB free on scratch — aborting to prevent failures" >&2
    echo "       Free up scratch space and resubmit" >&2
    exit 1
fi

# ── Detect read length → select pass 2 index ─────────────────────────────────
READ_LENGTH=$(zcat "${FQ1}" | \
    awk 'NR%4==2{print length($0); count++; if(count>=100) exit}' | \
    sort -n | tail -1)

if [[ "${READ_LENGTH}" -le 50 ]]; then
    INDEX_DIR="${INDEX_BASE}/gencode_v49_pass2_rl50"
else
    INDEX_DIR="${INDEX_BASE}/gencode_v49_pass2_rl100"
fi

echo "Read length : ${READ_LENGTH}bp -> index: $(basename ${INDEX_DIR})"

# ── Validate index exists ─────────────────────────────────────────────────────
if [[ ! -f "${INDEX_DIR}/SA" ]]; then
    echo "ERROR: Pass 2 index not found: ${INDEX_DIR}" >&2
    echo "Run: sbatch --export=READ_LENGTH=${READ_LENGTH} 03_build_pass2_index.sh" >&2
    exit 1
fi

# ── fastp trimming — write to scratch ────────────────────────────────────────
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
        --detect_adapter_for_pe 

else
    singularity exec \
        --bind /data1/abdelwao:/data1/abdelwao,${SCRATCH_BASE}:${SCRATCH_BASE} \
        "${FASTP_SIF}" fastp \
        --in1 "${FQ1}" \
        --out1 "${TRIMMED_R1}" \
        --thread "${SLURM_CPUS_PER_TASK}" \
        --qualified_quality_phred 20 \
        --length_required 25 

fi

if [[ $? -ne 0 ]]; then
    echo "ERROR: fastp failed for ${SAMPLE_ID}" >&2
    exit 1
fi

# ── STAR pass 2 — final BAM output ───────────────────────────────────────────
echo ""
echo "── STAR pass 2 ─────────────────────────────────────────"

if [[ "${SEQ_TYPE}" == "paired" ]]; then
    READ_FILES="${TRIMMED_R1} ${TRIMMED_R2}"
else
    READ_FILES="${TRIMMED_R1}"
fi

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
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI AS NM MD \
    --outFilterIntronMotifs None \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterMismatchNmax 999 \
    --quantMode TranscriptomeSAM \
    --outBAMsortingBinsN 200 \
    --limitBAMsortRAM 120000000000 \
    --limitSjdbInsertNsj 5000000

if [[ $? -ne 0 ]]; then
    echo "ERROR: STAR pass 2 failed for ${SAMPLE_ID}" >&2
    exit 1
fi

# ── Index BAM ─────────────────────────────────────────────────────────────────
echo ""
echo "── Indexing BAM ────────────────────────────────────────"
singularity exec \
    --bind /data1/abdelwao:/data1/abdelwao \
    "${SAMTOOLS_SIF}" samtools index \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"

if [[ $? -ne 0 ]]; then
    echo "ERROR: samtools index failed for ${SAMPLE_ID}" >&2
    exit 1
fi

# ── Verify outputs ────────────────────────────────────────────────────────────
BAM="${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"
BAI="${BAM}.bai"
SJ="${OUT_DIR}/${SAMPLE_ID}_SJ.out.tab"

[[ ! -f "${BAM}" ]] && { echo "ERROR: BAM not found" >&2; exit 1; }
[[ ! -f "${BAI}" ]] && { echo "ERROR: BAI not found" >&2; exit 1; }
[[ ! -f "${SJ}" ]]  && { echo "ERROR: SJ.out.tab not found" >&2; exit 1; }

BAM_SIZE=$(du -sh "${BAM}" | cut -f1)
SJ_COUNT=$(wc -l < "${SJ}")

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "── Pass 2 complete ─────────────────────────────────────"
echo "BAM        : ${BAM} (${BAM_SIZE})"
echo "BAI        : ${BAI}"
echo "SJ.out.tab : ${SJ} (${SJ_COUNT} junctions)"
echo "Log        : ${OUT_DIR}/${SAMPLE_ID}_Log.final.out"
echo ""
grep -E "Uniquely mapped reads %|Number of splices: Total|GT/AG|AT/AC|Non-canonical" \
    "${OUT_DIR}/${SAMPLE_ID}_Log.final.out"
echo ""
echo "Completed: ${SAMPLE_ID}"