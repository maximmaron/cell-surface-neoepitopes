#!/bin/bash
# =============================================================================
# 02b_build_groups_file.sh
# Build the LeafCutter groups file from BAM lists.
#
# Output: two-column TSV (sample_id <TAB> group)
#   col1: sample ID matching counts matrix column names (cases_sample / controls_sample)
#   col2: group label derived from COMPARISON (e.g. SF3B1_vs_normal -> SF3B1 / normal)
#         Control written first = baseline in LeafCutter DS
#
# Run AFTER clustering so the verification step can check against the
# actual counts matrix header.
# Usage: COMPARISON=SF3B1_vs_normal bash 02b_build_groups_file.sh
# =============================================================================

set -euo pipefail

COMPARISON="${COMPARISON:?COMPARISON env var must be set (e.g. SF3B1_vs_normal)}"

# Derive human-readable group labels by splitting COMPARISON on _vs_
# SF3B1_vs_normal -> CASE_GROUP=SF3B1, CTRL_GROUP=normal
CASE_GROUP="${COMPARISON%%_vs_*}"
CTRL_GROUP="${COMPARISON##*_vs_}"

# BAM list filenames are always cases/controls as written by 00_prep_bam_lists.sh
CASE_LABEL="cases"
CTRL_LABEL="controls"

WORK_DIR="/data1/abdelwao/maxim/splicing_pipeline/leafcutter/${COMPARISON}"
LIST_DIR="${WORK_DIR}/bam_lists"
CLUSTER_DIR="${WORK_DIR}/clustering"
DS_DIR="${WORK_DIR}/ds_input"
mkdir -p "${DS_DIR}"

CASE_BAMS="${LIST_DIR}/${CASE_LABEL}_bams.txt"
CTRL_BAMS="${LIST_DIR}/${CTRL_LABEL}_bams.txt"
COUNTS_GZ="${CLUSTER_DIR}/${COMPARISON}_perind_numers.counts.gz"
LC_GROUPS_FILE="${DS_DIR}/${COMPARISON}_groups.txt"

echo "Comparison : ${COMPARISON}"
echo "Case group : ${CASE_GROUP} (from ${CASE_BAMS})"
echo "Ctrl group : ${CTRL_GROUP} (from ${CTRL_BAMS})"
echo ""

for f in "${CASE_BAMS}" "${CTRL_BAMS}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: BAM list not found: ${f}" >&2
        echo "Run 00_prep_bam_lists.sh first." >&2
        exit 1
    fi
done

# ── Build groups file — control first (baseline) ──────────────────────────────
> "${LC_GROUPS_FILE}"

while IFS= read -r bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    printf "%s_%s\t%s\n" "${CTRL_LABEL}" "${sample}" "${CTRL_GROUP}"
done < "${CTRL_BAMS}" >> "${LC_GROUPS_FILE}"

while IFS= read -r bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    printf "%s_%s\t%s\n" "${CASE_LABEL}" "${sample}" "${CASE_GROUP}"
done < "${CASE_BAMS}" >> "${LC_GROUPS_FILE}"

N=$(wc -l < "${LC_GROUPS_FILE}")
echo "Groups file written: ${LC_GROUPS_FILE} (${N} samples)"
echo ""
echo "Preview:"
head -3 "${LC_GROUPS_FILE}"
echo "..."
tail -3 "${LC_GROUPS_FILE}"

# ── Verify against counts matrix header ───────────────────────────────────────
if [[ ! -f "${COUNTS_GZ}" ]]; then
    echo ""
    echo "WARNING: counts matrix not found at ${COUNTS_GZ} — skipping verification"
    exit 0
fi

echo ""
echo "Verifying sample IDs against counts matrix header..."

mapfile -t MAT_COLS < <(zcat "${COUNTS_GZ}" | head -1 | tr ' ' '\n' | grep -v '^$')
N_MAT=${#MAT_COLS[@]}
echo "Samples in counts matrix: ${N_MAT}"

N_MATCHED=0
N_UNMATCHED=0
UNMATCHED_LIST=()

while IFS=$'\t' read -r sample_id group; do
    if printf '%s\n' "${MAT_COLS[@]}" | grep -qx "${sample_id}"; then
        (( N_MATCHED++ ))
    else
        (( N_UNMATCHED++ ))
        UNMATCHED_LIST+=("${sample_id}")
    fi
done < "${LC_GROUPS_FILE}"

echo "Matched  : ${N_MATCHED}"
echo "Unmatched: ${N_UNMATCHED}"

if [[ ${N_UNMATCHED} -gt 0 ]]; then
    echo ""
    echo "WARNING — unmatched sample IDs (first 10):"
    printf '  %s\n' "${UNMATCHED_LIST[@]:0:10}"
    echo ""
    echo "First 5 matrix column names for comparison:"
    printf '  %s\n' "${MAT_COLS[@]:0:5}"
    exit 1
else
    echo ""
    echo "All ${N_MATCHED} sample IDs matched successfully."
    echo "Ready for: sbatch --export=COMPARISON=${COMPARISON} 03_differential_splicing.sh"
fi