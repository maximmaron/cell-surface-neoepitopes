#!/bin/bash
# =============================================================================
# 02b_build_groups_file.sh
# Build the LeafCutter groups file from BAM lists.
#
# Output: two-column TSV (sample_id <TAB> group)
#   col1: sample ID matching counts matrix column names (cohort_sample)
#   col2: group label (Normal or SRSF2) — Normal first = baseline
#
# Run AFTER clustering so the verification step can check against the
# actual counts matrix header.
# Usage: bash 02b_build_groups_file.sh
# =============================================================================

SRSF2_BAMS="/data1/abdelwao/shared/splicing_analysis/leafcutter/SRSF2_vs_normal/bam_lists/srsf2_bams.txt"
NORMAL_BAMS="/data1/abdelwao/shared/splicing_analysis/leafcutter/SRSF2_vs_normal/bam_lists/normal_bams.txt"
COUNTS_GZ="/scratch/abdelwao/leafcutter_srsf2/clustering/srsf2_vs_normal_perind_numers.counts.gz"
OUT_DIR="/scratch/abdelwao/leafcutter_srsf2/ds_input"
mkdir -p "${OUT_DIR}"

LC_GROUPS_FILE="${OUT_DIR}/srsf2_vs_normal_groups.txt"

# ── Build groups file — Normal first (baseline) ───────────────────────────────
> "${LC_GROUPS_FILE}"   # truncate/create

while IFS= read -r bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    printf "normal_%s\tNormal\n" "${sample}"
done < "${NORMAL_BAMS}" >> "${LC_GROUPS_FILE}"

while IFS= read -r bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    printf "srsf2_%s\tSRSF2\n" "${sample}"
done < "${SRSF2_BAMS}" >> "${LC_GROUPS_FILE}"

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

# Extract column names from header (skip first empty field from split)
mapfile -t MAT_COLS < <(zcat "${COUNTS_GZ}" | head -1 | tr ' ' '\n' | grep -v '^$')
N_MAT=${#MAT_COLS[@]}
echo "Samples in counts matrix: ${N_MAT}"

# Check each groups file ID against matrix columns
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
    echo "Ready for differential splicing."
fi