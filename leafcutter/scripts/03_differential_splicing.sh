#!/bin/bash
# =============================================================================
# 03_differential_splicing.sh
# SLURM job: run leafcutter-ds (Python implementation)
# Differential intron excision between case and control groups.
#
# Usage:  sbatch --export=COMPARISON=SF3B1_vs_normal 03_differential_splicing.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=lc_ds
#SBATCH --output=logs/lc_ds_%j.out
#SBATCH --error=logs/lc_ds_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

set -euo pipefail

source /admin/software/anaconda/1.11.1/bin/activate
conda activate leafcutter

COMPARISON="${COMPARISON:?COMPARISON env var must be set (e.g. --export=COMPARISON=SF3B1_vs_normal)}"

# Derive case group label for -0 flag (deltaPSI = case - control, positive = gained in mutant)
CASE_GROUP="${COMPARISON%%_vs_*}"

WORK_DIR="/data1/abdelwao/maxim/splicing_pipeline/leafcutter/${COMPARISON}"
CLUST_DIR="${WORK_DIR}/clustering"
DS_DIR="${WORK_DIR}/ds_output"
mkdir -p "${DS_DIR}" logs

# ── Inputs ────────────────────────────────────────────────────────────────────
COUNTS="${CLUST_DIR}/${COMPARISON}_perind_numers.counts.gz"
LC_GROUPS="${WORK_DIR}/ds_input/${COMPARISON}_groups.txt"

# Exon annotation file for cluster-to-gene labeling (optional but recommended)
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
EXONS="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode_v49_exons_header.txt.gz"

for f in "${COUNTS}" "${LC_GROUPS}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input not found: ${f}" >&2
        exit 1
    fi
done

echo "Comparison : ${COMPARISON}"
echo "Case group : ${CASE_GROUP} (baseline for deltaPSI direction)"
echo "Counts     : ${COUNTS}"
echo "Groups     : ${LC_GROUPS}"
echo ""

# ── Run leafcutter-ds ─────────────────────────────────────────────────────────
# --num_threads     : parallel threads; match cpus-per-task
# -c MIN_COVERAGE   : minimum total reads per cluster (default 20)
# -i                : minimum samples with >=1 read per intron (default 5)
# -g                : minimum samples per group with min_coverage (default 3)
# -0 CASE_GROUP     : baseline group for effect size direction
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
        -0 "${CASE_GROUP}" \
        -e "${EXONS}" \
        -o "${COMPARISON}" \
        "${COUNTS}" \
        "${LC_GROUPS}"
else
    echo "Exon annotation not found — running without gene labeling"
    leafcutter-ds \
        --num_threads 16 \
        -c 20 \
        -i 5 \
        -g 3 \
        -0 "${CASE_GROUP}" \
        -o "${COMPARISON}" \
        "${COUNTS}" \
        "${LC_GROUPS}"
fi

echo ""
echo "Differential splicing complete."
ls -lh "${DS_DIR}/${COMPARISON}"*
echo ""
echo "Next: run 04_leafcutter_analysis.R"

conda deactivate