#!/bin/bash
# =============================================================================
# submit_bam_qc.sh
#
# Builds the BAM manifest from sf_final_metadata_with_del.tsv (col 10 = bam_path),
# then submits the QC array job and a dependent aggregation job.
#
# This runs once across all samples before any rMATS comparison.
# The resulting bam_summary.tsv is the QC lookup table used by
# submit_rmats_persample.sh for all downstream comparisons.
#
# Usage:
#   bash submit_bam_qc.sh -m /path/to/sf_final_metadata_with_del.tsv
#
# Note on BED12 generation from GENCODE v49 GTF (run once):
#   gtfToGenePred gencode.v49.annotation.gtf /dev/stdout \
#       | genePredToBed /dev/stdin gencode.v49.annotation.bed12
# =============================================================================

set -euo pipefail

# --- Fixed paths ---
BED_FILE=/data1/abdelwao/maxim/annotations/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed12
OUT_DIR=/data1/abdelwao/maxim/splicing_pipeline/rmats/bam_qc
ARRAY_SCRIPT=/data1/abdelwao/maxim/splicing_pipeline/rmats/scripts/extract_bam_info.sh
AGGREGATE_SCRIPT=/data1/abdelwao/maxim/splicing_pipeline/rmats/scripts/aggregate_bam_qc.sh

# --- Argument parsing ---
metadata_tsv=""

usage() {
    echo "Usage: $0 -m /path/to/sf_final_metadata_with_del.tsv"
    exit 1
}

while getopts "m:h" opt; do
    case $opt in
        m) metadata_tsv="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -z "$metadata_tsv" || ! -f "$metadata_tsv" ]] && { echo "ERROR: metadata TSV not found: '$metadata_tsv'"; usage; }

# --- Setup ---
mkdir -p "${OUT_DIR}"/{strand,flagstat,per_sample,logs}

# --- Build BAM manifest from metadata TSV ---
# bam_path is col 10; skip header, skip blank/dot entries, deduplicate
MANIFEST=${OUT_DIR}/bam_manifest.txt

awk -F'\t' 'NR>1 && $10 != "" && $10 != "." {print $10}' "$metadata_tsv" \
    | sort -u \
    > "$MANIFEST"

N_BAMS=$(wc -l < "$MANIFEST")

if [[ "$N_BAMS" -eq 0 ]]; then
    echo "ERROR: No BAM paths found in column 10 of $metadata_tsv"
    echo "Check that bam_path column is populated and paths are not '.'"
    exit 1
fi

# Warn about any BAMs listed in the TSV that don't exist on disk
n_missing=0
while IFS= read -r bam; do
    if [[ ! -f "$bam" ]]; then
        echo "WARNING: BAM not found on disk: $bam"
        n_missing=$((n_missing + 1))
    fi
done < "$MANIFEST"

if [[ "$n_missing" -gt 0 ]]; then
    echo ""
    echo "WARNING: $n_missing BAM(s) listed in metadata are missing on disk."
    echo "These will cause array tasks to fail. Fix paths in metadata TSV before proceeding."
    echo ""
fi

echo "Found $N_BAMS BAMs in metadata TSV"
echo "Manifest written to: $MANIFEST"
echo ""

# Print breakdown by sample_group
echo "Sample group breakdown:"
awk -F'\t' 'NR>1 && $10 != "" && $10 != "." {print $5}' "$metadata_tsv" \
    | sort | uniq -c | sed 's/^/  /'
echo ""

MAX_INDEX=$((N_BAMS - 1))

# --- Submit array job ---
# Pass explicit --output/--error so logs land in OUT_DIR/logs/, not relative to
# wherever this script is called from.
ARRAY_JOB_ID=$(sbatch \
    --array=0-${MAX_INDEX} \
    --output="${OUT_DIR}/logs/bam_qc_%A_%a.out" \
    --error="${OUT_DIR}/logs/bam_qc_%A_%a.err" \
    --export=ALL,BED_FILE="${BED_FILE}",OUT_DIR="${OUT_DIR}" \
    "$ARRAY_SCRIPT" \
    | awk '{print $NF}')

echo "Submitted array job: $ARRAY_JOB_ID (${N_BAMS} tasks)"

# --- Submit aggregation job dependent on array completion ---
AGG_JOB_ID=$(sbatch \
    --dependency=afterok:${ARRAY_JOB_ID} \
    --job-name=bam_qc_aggregate \
    --cpus-per-task=1 \
    --mem=4G \
    --time=00:30:00 \
    --output="${OUT_DIR}/logs/aggregate_%j.out" \
    --error="${OUT_DIR}/logs/aggregate_%j.err" \
    --export=ALL,OUT_DIR="${OUT_DIR}" \
    "$AGGREGATE_SCRIPT" \
    | awk '{print $NF}')

echo "Submitted aggregation job: $AGG_JOB_ID (runs after array completes)"
echo ""
echo "Monitor with: squeue -u $USER"
echo ""
echo "Output: ${OUT_DIR}/bam_summary.tsv"