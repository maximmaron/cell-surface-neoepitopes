#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_download_%A_%a.out
#SBATCH --error=logs/sra_download_%A_%a.err
#SBATCH --array=1-72  # Adjust based on number of SRR accessions
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu

# Create logs directory if it doesn't exist
mkdir -p logs

# Configuration
SRR_LIST="$PWD/srr_list.txt"  # File with one SRR accession per line
OUTPUT_DIR="$PWD/srr_downloads"  # Change to your output directory
TEMP_DIR="$PWD/tmp"  # Temporary directory for large files
NGC_FILE="/data1/abdelwao/shared/splicing_analysis/dbGaP/prj_41432.ngc"  # Your dbGaP repository key

# Create output and temp directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

# Load SRA Toolkit module (adjust based on your cluster)
source /admin/software/anaconda/1.11.1/bin/activate
conda activate sratools

# Get the SRR accession for this array task
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SRR_LIST})

echo "================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing: ${SRR}"
echo "Started at: $(date)"
echo "================================================"

# Verify NGC file exists
if [ ! -f "${NGC_FILE}" ]; then
    echo "ERROR: NGC file not found at ${NGC_FILE}"
    exit 1
fi

export TMPDIR=${TEMP_DIR}

# Step 1: Prefetch the SRA file using NGC
echo "Starting prefetch for ${SRR}..."
prefetch --ngc ${NGC_FILE} ${SRR} \
    --max-size 100G \
    --output-directory ${TEMP_DIR}

# Check if prefetch was successful
if [ $? -ne 0 ]; then
    echo "ERROR: Prefetch failed for ${SRR}"
    exit 1
fi

echo "Prefetch completed for ${SRR}"

# Step 2: Rename the downloaded SRA file if it has dbGaP suffix
echo "Checking for dbGaP-renamed SRA file..."
DBGAP_SRA=$(find ${TEMP_DIR}/${SRR} -name "${SRR}_dbgap_*.sra" 2>/dev/null | head -n 1)
if [ -n "$DBGAP_SRA" ]; then
    echo "Found dbGaP SRA file: $DBGAP_SRA"
    echo "Renaming to ${SRR}.sra..."
    mv "$DBGAP_SRA" "${TEMP_DIR}/${SRR}/${SRR}.sra"
fi

# Step 3: Convert to FASTQ using fasterq-dump with NGC
echo "Converting ${SRR} to FASTQ..."
fasterq-dump --ngc ${NGC_FILE} ${TEMP_DIR}/${SRR}/${SRR}.sra \
    --split-files \
    --threads ${SLURM_CPUS_PER_TASK} \
    --outdir ${OUTPUT_DIR} \
    --temp ${TEMP_DIR} \
    --progress

# Check if conversion was successful
if [ $? -ne 0 ]; then
    echo "ERROR: FASTQ conversion failed for ${SRR}"
    exit 1
fi

echo "FASTQ conversion completed for ${SRR}"

# Step 3: Compress FASTQ files
echo "Compressing FASTQ files with pigz (using ${SLURM_CPUS_PER_TASK} threads)..."
for fastq in ${OUTPUT_DIR}/${SRR}*.fastq; do
    if [ -f "$fastq" ]; then
        echo "Compressing $fastq..."
        pigz -p ${SLURM_CPUS_PER_TASK} -9 $fastq
    fi
done

# Step 4: Validate the download (optional but recommended)
echo "Validating ${SRR}..."
vdb-validate ${TEMP_DIR}/${SRR}/${SRR}.sra

if [ $? -eq 0 ]; then
    echo "Validation successful for ${SRR}"
    # Clean up prefetched SRA file after successful validation
    echo "Cleaning up prefetched SRA file..."
    rm -rf ${TEMP_DIR}/${SRR}
else
    echo "WARNING: Validation failed for ${SRR} - keeping SRA file for inspection"
fi

# Step 5: Generate MD5 checksums
echo "Generating MD5 checksums..."
cd ${OUTPUT_DIR}
for file in ${SRR}*.fastq.gz; do
    if [ -f "$file" ]; then
        md5sum $file > ${file}.md5
    fi
done

echo "================================================"
echo "Completed processing: ${SRR}"
echo "Finished at: $(date)"
echo "================================================"

# Print file sizes for verification
echo "Output files:"
ls -lh ${OUTPUT_DIR}/${SRR}*

exit 0
