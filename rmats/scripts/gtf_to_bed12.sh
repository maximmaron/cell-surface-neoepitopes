#!/bin/bash
#SBATCH --job-name=gtf_to_bed12
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/gtf_to_bed12_%j.out
#SBATCH --error=logs/gtf_to_bed12_%j.err

# =============================================================================
# gtf_to_bed12.sh
# Converts GENCODE GTF to BED12 format for RSeQC infer_experiment.py
# =============================================================================

# --- Paths ---
GTF=/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf
OUT_DIR=/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/
GENEPRED=${OUT_DIR}/genes.genePred
BED12=${OUT_DIR}/genes.bed12

mkdir -p $OUT_DIR
mkdir -p logs

# --- Load/activate tools ---
source /admin/software/anaconda/1.11.1/bin/activate
conda activate ucsc-tools


# --- Check inputs ---
if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF not found: $GTF"
    exit 1
fi

echo "Input GTF : $GTF"
echo "Output BED: $BED12"
echo ""

# --- Step 1: GTF to genePred ---
echo "[$(date)] Running gtfToGenePred..."
gtfToGenePred \
    -genePredExt \
    -geneNameAsName2 \
    -ignoreGroupsWithoutExons \
    "$GTF" \
    "$GENEPRED"

if [[ $? -ne 0 || ! -s "$GENEPRED" ]]; then
    echo "ERROR: gtfToGenePred failed or produced empty output"
    exit 1
fi
echo "genePred records: $(wc -l < $GENEPRED)"

# --- Step 2: genePred to BED12 ---
echo "[$(date)] Running genePredToBed..."
genePredToBed \
    "$GENEPRED" \
    "$BED12"

if [[ $? -ne 0 || ! -s "$BED12" ]]; then
    echo "ERROR: genePredToBed failed or produced empty output"
    exit 1
fi
echo "BED12 records: $(wc -l < $BED12)"

# --- Step 3: Verify chromosome naming matches BAM files ---
echo ""
echo "Chromosome names in BED12 (first 5):"
cut -f1 "$BED12" | sort -u | head -5

echo ""
echo "Done. BED12 written to: $BED12"
echo ""
echo "IMPORTANT: Verify chromosome naming matches your BAM files before running RSeQC:"
echo "  cut -f1 $BED12 | sort -u | head -5"
echo "  samtools view -H your_sample.bam | grep '^@SQ' | cut -f2 | head -5"

conda deactivate