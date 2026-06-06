#!/bin/bash
#SBATCH --job-name=bam_qc_aggregate
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# NOTE: --output and --error are passed explicitly by submit_bam_qc.sh.

# =============================================================================
# aggregate_bam_qc.sh
# Runs after the array job completes.
# Combines per-sample TSVs into bam_summary.tsv.
#
# Passed via --export from submit_bam_qc.sh:
#   OUT_DIR
# =============================================================================

OUT_DIR=${OUT_DIR:-/data1/abdelwao/maxim/splicing_pipeline/rmats/bam_qc}
SUMMARY=${OUT_DIR}/bam_summary.tsv

echo "Aggregating results from: ${OUT_DIR}/per_sample/"
echo ""

# =============================================================================
# 1. Combine all per-sample TSVs
# =============================================================================

# Guard: fail clearly if no per-sample TSVs exist
shopt -s nullglob
tsv_files=("${OUT_DIR}/per_sample/"*.tsv)
shopt -u nullglob

if [[ "${#tsv_files[@]}" -eq 0 ]]; then
    echo "ERROR: No per-sample TSV files found in ${OUT_DIR}/per_sample/"
    echo "Check that the array job completed successfully before running aggregation."
    exit 1
fi

# Write header from first file, then data rows from all files
head -1 "${tsv_files[0]}" > "$SUMMARY"
for f in "${tsv_files[@]}"; do
    tail -n +2 "$f" >> "$SUMMARY"
done

echo "Combined summary: $SUMMARY"
echo "Total samples: $(( $(wc -l < "$SUMMARY") - 1 ))"
echo ""

# =============================================================================
# 2. Print summary table
# =============================================================================
echo "=============================================================================================="
echo "SAMPLE SUMMARY"
echo "=============================================================================================="
column -t -s $'\t' "$SUMMARY"
echo ""

# =============================================================================
# 3. Flag samples needing review
# =============================================================================
echo "=============================================================================================="
echo "SAMPLES REQUIRING REVIEW"
echo "=============================================================================================="

any_flagged=0
while IFS=$'\t' read -r sample bam maxlen medlen paired total mapped rate fwd rev undet strand libtype lib_flag; do
    [[ "$sample" == "sample" ]] && continue

    if [[ "$fwd" != "NA" && "$rev" != "NA" ]]; then
        is_ambiguous=$(echo "$fwd > 0.6 && $fwd < 0.75" | bc -l 2>/dev/null || echo 0)
        if [[ "$is_ambiguous" -eq 1 ]]; then
            echo "  [AMBIGUOUS STRAND] $sample  fwd=${fwd}  rev=${rev}"
            any_flagged=1
        fi
    fi

    if [[ "$rate" != "NA" ]]; then
        low_map=$(echo "$rate < 70" | bc -l 2>/dev/null || echo 0)
        if [[ "$low_map" -eq 1 ]]; then
            echo "  [LOW MAPPING RATE] $sample  rate=${rate}%"
            any_flagged=1
        fi
    fi

    if [[ "$paired" == "no" ]]; then
        echo "  [SINGLE-END]       $sample"
        any_flagged=1
    fi

done < "$SUMMARY"

[[ "$any_flagged" -eq 0 ]] && echo "  None"
echo ""

# =============================================================================
# 4. Library type breakdown — informational only
# submit_rmats_persample.sh handles per-sample libtype lookup directly
# =============================================================================
echo "=============================================================================================="
echo "LIBRARY TYPE BREAKDOWN"
echo "=============================================================================================="

declare -A GROUPS

while IFS=$'\t' read -r sample bam maxlen medlen paired total mapped rate fwd rev undet strand libtype lib_flag; do
    [[ "$sample" == "sample" ]] && continue
    group_key="${libtype}|${lib_flag}"
    if [[ -z "${GROUPS[$group_key]+_}" ]]; then
        GROUPS[$group_key]="$sample"
    else
        GROUPS[$group_key]="${GROUPS[$group_key]},$sample"
    fi
done < "$SUMMARY"

for key in "${!GROUPS[@]}"; do
    libtype=$(echo "$key" | cut -d'|' -f1)
    lib_flag=$(echo "$key" | cut -d'|' -f2)
    samples="${GROUPS[$key]}"
    n_samples=$(echo "$samples" | tr ',' '\n' | wc -l)
    echo "  libType=${libtype}  ${lib_flag}  (${n_samples} samples)"
    echo "  Samples: $(echo "$samples" | tr ',' '\n' | head -5 | tr '\n' ' ')$([ "$n_samples" -gt 5 ] && echo '...')"
    echo ""
done

echo "Done. Review ${SUMMARY} before running rMATS comparisons."