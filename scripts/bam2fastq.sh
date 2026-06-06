#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=bam2fastq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb
#SBATCH --time=48:00:00
#SBATCH --output=logs/bam2fastq_%A_%a.log

###############################################################################
# BAM to FASTQ Conversion with Sample ID Lookup
# Usage: sbatch --array=0-N bam2fastq.sh ##N is number of BAM files - 1
###############################################################################

source /admin/software/anaconda/1.11.1/bin/activate
conda activate samtools

set -euo pipefail

# Configuration - UPDATE THESE
BAM_DIR="/data1/abdelwao/shared/splicing_analysis/dbGaP/beat/bam_files"
OUTPUT_DIR="fastq_output"
SAMPLESHEET="gdc_sample_sheet_merged.txt"  # Your samplesheet file
THREADS=$SLURM_CPUS_PER_TASK

# Get all BAM files
BAM_FILES=($BAM_DIR/*.bam)
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# Get BAM filename (without path)
BAM_FILENAME=$(basename "$BAM_FILE")

echo "============================================"
echo "Processing: $BAM_FILENAME"
echo "============================================"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Look up Sample ID from samplesheet
# Search for the BAM filename in File.Name column and extract Sample.ID
if [ -f "$SAMPLESHEET" ]; then
    SAMPLE_ID=$(awk -v bam="$BAM_FILENAME" 'BEGIN{FS="\t"} NR>1 && $2==bam {print $7; exit}' "$SAMPLESHEET")
    
    if [ -z "$SAMPLE_ID" ]; then
        echo "WARNING: BAM file not found in samplesheet: $BAM_FILENAME"
        echo "Using BAM filename as sample ID"
        SAMPLE_ID=$(basename "$BAM_FILE" .bam)
    else
        echo "Found Sample ID: $SAMPLE_ID"
    fi
else
    echo "WARNING: Samplesheet not found: $SAMPLESHEET"
    echo "Using BAM filename as sample ID"
    SAMPLE_ID=$(basename "$BAM_FILE" .bam)
fi

echo "Final Sample ID: $SAMPLE_ID"
echo ""

# Check for pigz
if command -v pigz &> /dev/null; then
    GZIP_CMD="pigz -p $THREADS"
    echo "Using pigz for compression"
else
    GZIP_CMD="gzip"
    echo "Using gzip for compression"
fi

# Detect paired-end
echo "Detecting read type..."
PAIRED_COUNT=$(samtools view -c -f 1 "$BAM_FILE" | head -1)

if [ "$PAIRED_COUNT" -gt 0 ]; then
    # Paired-end
    echo "→ Paired-end detected"
    echo ""
    
    # Sort by name
    echo "Sorting BAM..."
    SORTED_BAM="${OUTPUT_DIR}/.${SAMPLE_ID}.sorted.bam"
    samtools sort -n -@ "$THREADS" -o "$SORTED_BAM" "$BAM_FILE"
    
    # Extract reads
    echo "Extracting reads..."
    samtools fastq -@ "$THREADS" \
        -1 >($GZIP_CMD > "${OUTPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz") \
        -2 >($GZIP_CMD > "${OUTPUT_DIR}/${SAMPLE_ID}_R2.fastq.gz") \
        -s >($GZIP_CMD > "${OUTPUT_DIR}/${SAMPLE_ID}_singleton.fastq.gz") \
        "$SORTED_BAM"
    
    # Remove sorted BAM
    rm -f "$SORTED_BAM"
    
    # Count reads
    R1_COUNT=$(zcat "${OUTPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz" | wc -l)
    R1_READS=$((R1_COUNT / 4))
    
    echo "✓ R1 reads: $R1_READS"
    echo "✓ R2 reads: $R1_READS"
    
else
    # Single-end
    echo "→ Single-end detected"
    echo ""
    
    echo "Extracting reads..."
    samtools fastq -@ "$THREADS" "$BAM_FILE" | $GZIP_CMD > "${OUTPUT_DIR}/${SAMPLE_ID}.fastq.gz"
    
    # Count reads
    READ_COUNT=$(zcat "${OUTPUT_DIR}/${SAMPLE_ID}.fastq.gz" | wc -l)
    READS=$((READ_COUNT / 4))
    
    echo "✓ Reads: $READS"
fi

echo ""
echo "============================================"
echo "Completed: $SAMPLE_ID"
echo "============================================"
echo "Output files:"
if [ "$PAIRED_COUNT" -gt 0 ]; then
    echo "  ${OUTPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz"
    echo "  ${OUTPUT_DIR}/${SAMPLE_ID}_R2.fastq.gz"
else
    echo "  ${OUTPUT_DIR}/${SAMPLE_ID}.fastq.gz"
fi


conda deactivate
