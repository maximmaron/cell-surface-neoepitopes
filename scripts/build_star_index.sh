#!/bin/bash
# =============================================================================
# 00b_build_star_index.sh
# Build STAR genome index for a specific read length.
# Run once per read length group before pass 1 alignment.
#
# Two indexes needed for your cohorts:
#   READ_LENGTH=50  → gencode_v49_rl50  (TCGA, Bodymap2)
#   READ_LENGTH=100 → gencode_v49_rl100 (BEAT, Leucegene, Pellagatti,
#                                        Madan, Maiga, HPA)
#
# Usage:
#   sbatch --export=READ_LENGTH=100 00b_build_star_index.sh
#   sbatch --export=READ_LENGTH=50  00b_build_star_index.sh
#
# Requirements: ~100G RAM, ~2 hours per index
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=star_index
#SBATCH --output=logs/star_index_%j.out
#SBATCH --error=logs/star_index_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16

STAR_SIF="/data1/abdelwao/maxim/containers/star_2.7.11b.sif"
GENOME_FA="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/GRCh38.primary_assembly.genome.fa"
GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"

mkdir -p logs

# ── Validate READ_LENGTH ──────────────────────────────────────────────────────
if [[ -z "${READ_LENGTH}" ]]; then
    echo "ERROR: READ_LENGTH not set." >&2
    echo "Usage: sbatch --export=READ_LENGTH=100 00b_build_star_index.sh" >&2
    exit 1
fi

if ! [[ "${READ_LENGTH}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: READ_LENGTH must be an integer, got: ${READ_LENGTH}" >&2
    exit 1
fi

OVERHANG=$((READ_LENGTH - 1))
INDEX_DIR="/data1/abdelwao/maxim/annotations/Homo_sapiens/STAR/gencode_v49_rl${READ_LENGTH}"

echo "======================================="
echo "Building STAR index"
echo "  Read length   : ${READ_LENGTH}bp"
echo "  sjdbOverhang  : ${OVERHANG}"
echo "  Output        : ${INDEX_DIR}"
echo "  FASTA         : ${GENOME_FA}"
echo "  GTF           : ${GTF}"
echo "  Threads       : ${SLURM_CPUS_PER_TASK}"
echo "======================================="

# ── Skip if already built ─────────────────────────────────────────────────────
if [[ -f "${INDEX_DIR}/SA" && -f "${INDEX_DIR}/Genome" ]]; then
    echo "Index already exists — skipping build."
    echo "Delete ${INDEX_DIR} to force rebuild."
    exit 0
fi

mkdir -p "${INDEX_DIR}"

singularity exec \
    --bind /data1/abdelwao:/data1/abdelwao \
    "${STAR_SIF}" STAR \
    --runMode genomeGenerate \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "${INDEX_DIR}" \
    --genomeFastaFiles "${GENOME_FA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang "${OVERHANG}"

if [[ $? -eq 0 ]]; then
    echo ""
    echo "Index built successfully: ${INDEX_DIR}"
    ls -lh "${INDEX_DIR}/"
else
    echo "ERROR: Index build failed" >&2
    exit 1
fi