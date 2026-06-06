#!/bin/bash
# =============================================================================
# 03_build_pass2_index.sh
# Build STAR pass 2 genome index incorporating pooled novel junctions.
#
# Mirrors the same READ_LENGTH parameter pattern as 00b_build_star_index.sh.
#
# Usage:
#   sbatch --export=READ_LENGTH=100 03_build_pass2_index.sh
#   sbatch --export=READ_LENGTH=50  03_build_pass2_index.sh
#
# Optionally specify junction directory (versioned):
#   sbatch --export=READ_LENGTH=100,JUNC_DIR=/path/to/junctions/v1_20260325 \
#       03_build_pass2_index.sh
#
# If JUNC_DIR not set, auto-detects most recent v1_* directory.
# Junction file is always: {JUNC_DIR}/final_filtered_junctions.tab
#
# Read length groups:
#   READ_LENGTH=100 -> gencode_v49_pass2_rl100 (BEAT, Leucegene, Pellagatti,
#                                               Madan, Maiga, HPA)
#   READ_LENGTH=50  -> gencode_v49_pass2_rl50  (TCGA, Bodymap2)
#
# Requirements: ~100G RAM, ~2 hours per index
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=star_pass2_index
#SBATCH --output=logs/pass2_index_%j.out
#SBATCH --error=logs/pass2_index_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

# ── Containers ────────────────────────────────────────────────────────────────
STAR_SIF="/data1/abdelwao/maxim/containers/star_2.7.11b.sif"

# ── Reference ─────────────────────────────────────────────────────────────────
GENOME_FA="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/GRCh38.primary_assembly.genome.fa"
GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"
INDEX_BASE="/data1/abdelwao/maxim/annotations/Homo_sapiens/STAR"

# ── Junction directory ────────────────────────────────────────────────────────
# Points to versioned directory produced by 02_pool_junctions.sh
# Auto-detects most recent v1_* directory if JUNC_DIR not set
if [[ -z "${JUNC_DIR}" ]]; then
    JUNC_DIR=$(ls -td \
        /data1/abdelwao/maxim/splicing_pipeline/junctions/v1_* \
        2>/dev/null | head -1)
    echo "Auto-detected junction directory: ${JUNC_DIR}"
fi

JUNCTION_FILE="${JUNC_DIR}/final_filtered_junctions.tab"

mkdir -p logs

# ── Validate READ_LENGTH ──────────────────────────────────────────────────────
if [[ -z "${READ_LENGTH}" ]]; then
    echo "ERROR: READ_LENGTH not set." >&2
    echo "Usage: sbatch --export=READ_LENGTH=100 03_build_pass2_index.sh" >&2
    echo "       sbatch --export=READ_LENGTH=50  03_build_pass2_index.sh" >&2
    exit 1
fi

if ! [[ "${READ_LENGTH}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: READ_LENGTH must be an integer, got: ${READ_LENGTH}" >&2
    exit 1
fi

# ── Validate junction directory and file ──────────────────────────────────────
if [[ -z "${JUNC_DIR}" || ! -d "${JUNC_DIR}" ]]; then
    echo "ERROR: Junction directory not found: ${JUNC_DIR}" >&2
    echo "Run 02_pool_junctions.sh first, or set JUNC_DIR explicitly:" >&2
    echo "  sbatch --export=READ_LENGTH=${READ_LENGTH},JUNC_DIR=/path/to/v1_DATE \\" >&2
    echo "      03_build_pass2_index.sh" >&2
    exit 1
fi

if [[ ! -f "${JUNCTION_FILE}" ]]; then
    echo "ERROR: Junction file not found: ${JUNCTION_FILE}" >&2
    echo "Expected: ${JUNC_DIR}/final_filtered_junctions.tab" >&2
    exit 1
fi

OVERHANG=$((READ_LENGTH - 1))
INDEX_DIR="${INDEX_BASE}/gencode_v49_pass2_rl${READ_LENGTH}"
JUNC_COUNT=$(wc -l < "${JUNCTION_FILE}")

echo "======================================="
echo "Building STAR pass 2 index"
echo "  Read length   : ${READ_LENGTH}bp"
echo "  sjdbOverhang  : ${OVERHANG}"
echo "  Output        : ${INDEX_DIR}"
echo "  Junction dir  : ${JUNC_DIR}"
echo "  Junction file : ${JUNCTION_FILE}"
echo "  Junctions     : ${JUNC_COUNT}"
echo "  Genome FASTA  : ${GENOME_FA}"
echo "  GTF           : ${GTF}"
echo "  Threads       : ${SLURM_CPUS_PER_TASK}"
echo "======================================="
echo ""

# ── Skip if already built ─────────────────────────────────────────────────────
if [[ -f "${INDEX_DIR}/SA" && -f "${INDEX_DIR}/Genome" ]]; then
    echo "Index already exists — skipping build."
    echo "Delete ${INDEX_DIR} to force rebuild."
    exit 0
fi

mkdir -p "${INDEX_DIR}"

# ── Build index ───────────────────────────────────────────────────────────────
singularity exec \
    --bind /data1/abdelwao:/data1/abdelwao \
    "${STAR_SIF}" STAR \
    --runMode genomeGenerate \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "${INDEX_DIR}" \
    --genomeFastaFiles "${GENOME_FA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang "${OVERHANG}" \
    --sjdbFileChrStartEnd "${JUNCTION_FILE}" \
    --limitSjdbInsertNsj 5000000

if [[ $? -eq 0 ]]; then
    echo ""
    echo "Index built successfully: ${INDEX_DIR}"
    echo ""
    ls -lh "${INDEX_DIR}/"
    echo ""
    echo "Next steps:"
    echo "  Submit other read length if needed:"
    echo "    sbatch --export=READ_LENGTH=50 03_build_pass2_index.sh"
    echo "  Then run pass 2:"
    echo "    bash 04_star_pass2.sh"
else
    echo "ERROR: Index build failed for READ_LENGTH=${READ_LENGTH}" >&2
    exit 1
fi

echo "======================================="
echo "Completed: $(date)"