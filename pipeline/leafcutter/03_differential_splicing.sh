#!/bin/bash
# =============================================================================
# 03_differential_splicing.sh
# SLURM job: run leafcutter-ds (Python implementation)
# Differential intron excision: SRSF2 vs Normal
#
# Usage:  sbatch 03_differential_splicing.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=lc_ds
#SBATCH --output=logs/lc_ds_%j.out
#SBATCH --error=logs/lc_ds_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

source /admin/software/anaconda/1.11.1/bin/activate
conda activate leafcutter

WORK_DIR="/scratch/abdelwao/leafcutter_srsf2"
CLUST_DIR="${WORK_DIR}/clustering"
DS_DIR="${WORK_DIR}/ds_output"
mkdir -p "${DS_DIR}" logs

# ── Inputs ────────────────────────────────────────────────────────────────────
COUNTS="${CLUST_DIR}/srsf2_vs_normal_perind_numers.counts.gz"
LC_GROUPS="${WORK_DIR}/ds_input/srsf2_vs_normal_groups.txt"

# Exon annotation file for cluster-to-gene labeling (optional but recommended)
# Columns: chr, start, end, strand, gene_name
# Generate from GTF with: 
#   python -c "
#   import gzip, re
#   with gzip.open('gencode.v49.gtf.gz','rt') as f, open('gencode_v49_exons.txt','w') as o:
#       for line in f:
#           if '\texon\t' not in line: continue
#           f9 = line.split('\t')[8]
#           gene = re.search('gene_name \"([^\"]+)\"', f9)
#           if not gene: continue
#           c = line.split('\t')
#           o.write(f'{c[0]}\t{c[3]}\t{c[4]}\t{c[6]}\t{gene.group(1)}\n')
#   "
EXONS="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode_v49_exons.txt.gz"

# ── Run leafcutter-ds ─────────────────────────────────────────────────────────
# -p NUM_THREADS    : parallel threads; match cpus-per-task
# -c MIN_COVERAGE   : minimum total reads per cluster (default 20)
# -i                : minimum samples with >=1 read per intron (default 5)
# -g                : minimum samples per group with min_coverage (default 3)
# -0 SRSF2         : baseline group for effect size direction
#                     (deltaPSI = SRSF2 - Normal, positive = gained in mutant)
# -e EXON_FILE      : optional gene annotation for cluster labeling
# -o OUTPUT_PREFIX  : prefix for output files

cd "${DS_DIR}"

if [[ -f "${EXONS}" ]]; then
    echo "Using exon annotation: ${EXONS}"
    leafcutter-ds \
        --num_threads 16 \
        -c 20 \
        -i 5 \
        -g 3 \
        -0 SRSF2 \
        -e "${EXONS}" \
        -o "srsf2_vs_normal" \
        "${COUNTS}" \
        "${LC_GROUPS}"
else
    echo "Exon annotation not found — running without gene labeling"
    leafcutter-ds \
        --num_threads 16 \
        -c 20 \
        -i 5 \
        -g 3 \
        -0 SRSF2 \
        -o "srsf2_vs_normal" \
        "${COUNTS}" \
        "${LC_GROUPS}"
fi

echo ""
echo "Differential splicing complete."
ls -lh "${DS_DIR}"/srsf2_vs_normal*
echo ""
echo "Next: run 04_leafcutter_analysis.R"