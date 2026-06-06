#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --array=0-1165
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=cpu
#SBATCH --output=/data1/abdelwao/maxim/splicing_pipeline/gatk/bamcoverage/slurm_logs/bamcoverage_%A_%a.out
#SBATCH --error=/data1/abdelwao/maxim/splicing_pipeline/gatk/bamcoverage/slurm_logs/bamcoverage_%A_%a.err

source /admin/software/anaconda/1.11.1/bin/activate
conda activate deeptools

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "/data1/abdelwao/maxim/splicing_pipeline/gatk/bamcoverage/bam_list.tsv")
SAMPLE_ID=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE"       | cut -f2)

OUT_BW="$(dirname $BAM)/${SAMPLE_ID}.bw"

echo "Task $SLURM_ARRAY_TASK_ID | $SAMPLE_ID"
echo "  BAM: $BAM"
echo "  BW:  $OUT_BW"

bamCoverage \
    --bam                "$BAM" \
    --outFileName        "$OUT_BW" \
    --outFileFormat      bigwig \
    --binSize            1 \
    --normalizeUsing     CPM \
    --minMappingQuality  20 \
    --numberOfProcessors 16
