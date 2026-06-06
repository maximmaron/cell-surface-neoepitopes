#!/bin/bash
# =============================================================================
# pass2_sanity_check.sh
# Sanity check and QC summary for pass 2 BAM outputs.
#
# Checks per sample:
#   1. BAM exists and is non-empty
#   2. BAI exists and is non-empty
#   3. Log.final.out exists
#   4. Mapping rate above cohort-specific threshold
#   5. Junction count above minimum threshold
#   6. AT/AC splice count (ZRSR2 validation)
#   7. BAM file size reasonable
#
# Compares pass 2 vs pass 1 metrics where available:
#   - Mapping rate
#   - Junction count
#   - AT/AC splice count
#
# Outputs:
#   - Summary counts (complete/missing/flagged)
#   - Per-cohort mapping rate table
#   - Flagged samples with specific issues
#   - Resubmission command for missing samples
#
# Usage: bash pass2_sanity_check.sh
# =============================================================================

MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
OUT_DIR="/data1/abdelwao/maxim/splicing_pipeline/qc"
PASS1_QC="${OUT_DIR}/pass1_qc_summary_v3.tsv"

mkdir -p "${OUT_DIR}"

SUMMARY_FILE="${OUT_DIR}/pass2_qc_summary.tsv"
FLAGGED_FILE="${OUT_DIR}/pass2_flagged_samples.tsv"
MISSING_FILE="${OUT_DIR}/pass2_missing_samples.tsv"
REPORT="${OUT_DIR}/pass2_sanity_check.txt"

MIN_BAM_SIZE=104857600  # 100MB

# ── Thresholds ────────────────────────────────────────────────────────────────
declare -A MIN_MAPPING
MIN_MAPPING["beat"]=80
MIN_MAPPING["tcga"]=65
MIN_MAPPING["leucegene_SRSF2"]=80
MIN_MAPPING["leucegene_SF3B1"]=80
MIN_MAPPING["leucegene_U2AF1"]=80
MIN_MAPPING["leucegene_ZRSR2"]=80
MIN_MAPPING["pellagatti_SF3B1"]=75
MIN_MAPPING["pellagatti_SRSF2"]=75
MIN_MAPPING["pellagatti_U2AF1"]=75
MIN_MAPPING["pellagatti_WT"]=75
MIN_MAPPING["madan_2015_ZRSR2"]=75
MIN_MAPPING["madan_2015_WT"]=75
MIN_MAPPING["maiga_normal"]=75
MIN_MAPPING["madan_normal"]=75
MIN_MAPPING["hpa_2014"]=75
MIN_MAPPING["bodymap2"]=65

declare -A MIN_JUNCTIONS
MIN_JUNCTIONS["beat"]=500000
MIN_JUNCTIONS["tcga"]=100000
MIN_JUNCTIONS["leucegene_SRSF2"]=500000
MIN_JUNCTIONS["leucegene_SF3B1"]=500000
MIN_JUNCTIONS["leucegene_U2AF1"]=500000
MIN_JUNCTIONS["leucegene_ZRSR2"]=500000
MIN_JUNCTIONS["pellagatti_SF3B1"]=200000
MIN_JUNCTIONS["pellagatti_SRSF2"]=200000
MIN_JUNCTIONS["pellagatti_U2AF1"]=200000
MIN_JUNCTIONS["pellagatti_WT"]=200000
MIN_JUNCTIONS["madan_2015_ZRSR2"]=200000
MIN_JUNCTIONS["madan_2015_WT"]=200000
MIN_JUNCTIONS["maiga_normal"]=200000
MIN_JUNCTIONS["madan_normal"]=200000
MIN_JUNCTIONS["hpa_2014"]=200000
MIN_JUNCTIONS["bodymap2"]=100000

# ── Initialize output files ───────────────────────────────────────────────────
echo -e "task_num\tsample_id\tcohort\treason" > "${MISSING_FILE}"
echo -e "sample_id\tcohort\tmapping_pct\tjunctions\tatac_splices\tbam_size_mb\tflags" \
    > "${FLAGGED_FILE}"
echo -e "sample_id\tcohort\ttier\tmapping_pct\tjunctions\tatac_splices\tgt_ag_splices\tinput_reads\tbam_size_mb" \
    > "${SUMMARY_FILE}"

# ── Counters ──────────────────────────────────────────────────────────────────
TOTAL=0
COMPLETE=0
MISSING=0
FLAGGED=0
task_num=1
missing_tasks=()
missing_large=()
missing_other=()

# ── Per-cohort accumulators ───────────────────────────────────────────────────
declare -A COHORT_COUNT
declare -A COHORT_MAP_SUM
declare -A COHORT_MAP_MIN
declare -A COHORT_MAP_MAX
declare -A COHORT_JUNC_SUM
declare -A COHORT_ATAC_SUM
declare -A COHORT_BAM_SUM

# Redirect all output to both terminal and report file
exec > >(tee "${REPORT}") 2>&1

echo "======================================="
echo "Pass 2 sanity check — $(date)"
echo "Manifest: ${MANIFEST}"
echo "Report  : ${REPORT}"
echo "======================================="
echo ""

# ── Main loop ─────────────────────────────────────────────────────────────────
while IFS=$'\t' read -r sample_id cohort tier disease seq_type fq1 fq2; do
    ((TOTAL++))
    out_dir=$(dirname "${fq1}")
    bam="${out_dir}/${sample_id}_Aligned.sortedByCoord.out.bam"
    bai="${out_dir}/${sample_id}_Aligned.sortedByCoord.out.bam.bai"
    log="${out_dir}/${sample_id}_Log.final.out"

    # ── Check existence ───────────────────────────────────────────────────────
    if [[ ! -f "${bam}" || ! -f "${bai}" || ! -f "${log}" ]]; then
        reason="Missing"
        [[ ! -f "${bam}" ]] && reason="Missing BAM"
        [[ ! -f "${bai}" ]] && reason="Missing BAI"
        [[ ! -f "${log}" ]] && reason="Missing Log"
        echo -e "${task_num}\t${sample_id}\t${cohort}\t${reason}" >> "${MISSING_FILE}"
        missing_tasks+=("${task_num}")
        # Determine concurrency tier by compressed FASTQ size
        # Large: >3GB compressed (~30GB+ uncompressed) → needs low concurrency
        fq_size=$(stat -c%s "${fq1}" 2>/dev/null || echo 0)
        fq_gb=$(echo "scale=1; ${fq_size}/1024/1024/1024" | bc)
        if (( $(echo "${fq_gb} > 3" | bc -l) )); then
            missing_large+=("${task_num}")
        else
            missing_other+=("${task_num}")
        fi
        ((MISSING++))
        ((task_num++))
        continue
    fi

    # ── Check BAM size ────────────────────────────────────────────────────────
    bam_size=$(stat -c%s "${bam}" 2>/dev/null || echo 0)
    bam_mb=$(( bam_size / 1024 / 1024 ))

    if [[ "${bam_size}" -eq 0 ]]; then
        echo -e "${task_num}\t${sample_id}\t${cohort}\tEmpty BAM (0 bytes)" \
            >> "${MISSING_FILE}"
        missing_tasks+=("${task_num}")
        # Determine concurrency tier by compressed FASTQ size
        # Large: >3GB compressed (~30GB+ uncompressed) → needs low concurrency
        fq_size=$(stat -c%s "${fq1}" 2>/dev/null || echo 0)
        fq_gb=$(echo "scale=1; ${fq_size}/1024/1024/1024" | bc)
        if (( $(echo "${fq_gb} > 3" | bc -l) )); then
            missing_large+=("${task_num}")
        else
            missing_other+=("${task_num}")
        fi
        ((MISSING++))
        ((task_num++))
        continue
    fi

    # ── Extract QC metrics from Log.final.out ─────────────────────────────────
    mapping_pct=$(awk -F'|' '/Uniquely mapped reads %/{
        gsub(/ |%/,"",$2); print $2+0}' "${log}")
    input_reads=$(awk -F'|' '/Number of input reads/{
        gsub(/ /,"",$2); print $2+0}' "${log}")
    gt_ag=$(awk -F'|' '/Number of splices: GT\/AG/{
        gsub(/ /,"",$2); print $2+0}' "${log}")
    at_ac=$(awk -F'|' '/Number of splices: AT\/AC/{
        gsub(/ /,"",$2); print $2+0}' "${log}")
    sj_count=$(wc -l < "${out_dir}/${sample_id}_SJ.out.tab" 2>/dev/null || echo 0)

    # ── Write to summary ──────────────────────────────────────────────────────
    echo -e "${sample_id}\t${cohort}\t${tier}\t${mapping_pct}\t${sj_count}\t${at_ac}\t${gt_ag}\t${input_reads}\t${bam_mb}" \
        >> "${SUMMARY_FILE}"

    # ── Quality flags ─────────────────────────────────────────────────────────
    flags=""

    # Mapping rate
    min_map=${MIN_MAPPING[$cohort]:-75}
    if (( $(echo "${mapping_pct} < ${min_map}" | bc -l) )); then
        flags="${flags}LOW_MAPPING(${mapping_pct}%<${min_map}%);"
    fi

    # Junction count
    min_junc=${MIN_JUNCTIONS[$cohort]:-200000}
    if [[ "${sj_count}" -lt "${min_junc}" ]]; then
        flags="${flags}LOW_JUNCTIONS(${sj_count}<${min_junc});"
    fi

    # Small BAM
    if [[ "${bam_size}" -lt "${MIN_BAM_SIZE}" ]]; then
        flags="${flags}SMALL_BAM(${bam_mb}MB);"
    fi

    # AT/AC — warn if very low for ZRSR2 cohort
    if [[ "${cohort}" =~ zrsr2|ZRSR2 ]] && [[ "${at_ac}" -lt 10000 ]]; then
        flags="${flags}LOW_ATAC(${at_ac});"
    fi

    if [[ -n "${flags}" ]]; then
        echo -e "${sample_id}\t${cohort}\t${mapping_pct}\t${sj_count}\t${at_ac}\t${bam_mb}\t${flags}" \
            >> "${FLAGGED_FILE}"
        ((FLAGGED++))
    fi

    # ── Accumulate cohort stats ───────────────────────────────────────────────
    COHORT_COUNT[$cohort]=$(( ${COHORT_COUNT[$cohort]:-0} + 1 ))
    COHORT_MAP_SUM[$cohort]=$(echo "${COHORT_MAP_SUM[$cohort]:-0} + ${mapping_pct}" | bc)
    COHORT_JUNC_SUM[$cohort]=$(( ${COHORT_JUNC_SUM[$cohort]:-0} + sj_count ))
    COHORT_ATAC_SUM[$cohort]=$(( ${COHORT_ATAC_SUM[$cohort]:-0} + at_ac ))
    COHORT_BAM_SUM[$cohort]=$(( ${COHORT_BAM_SUM[$cohort]:-0} + bam_mb ))

    if [[ -z "${COHORT_MAP_MIN[$cohort]}" ]] || \
       (( $(echo "${mapping_pct} < ${COHORT_MAP_MIN[$cohort]}" | bc -l) )); then
        COHORT_MAP_MIN[$cohort]="${mapping_pct}"
    fi
    if [[ -z "${COHORT_MAP_MAX[$cohort]}" ]] || \
       (( $(echo "${mapping_pct} > ${COHORT_MAP_MAX[$cohort]}" | bc -l) )); then
        COHORT_MAP_MAX[$cohort]="${mapping_pct}"
    fi

    ((COMPLETE++))
    ((task_num++))

done < <(tail -n +2 "${MANIFEST}")

# ── Overall summary ───────────────────────────────────────────────────────────
echo "── Overall ─────────────────────────────────────────────"
echo "  Total samples : ${TOTAL}"
echo "  Complete      : ${COMPLETE}"
echo "  Missing       : ${MISSING}"
echo "  Flagged       : ${FLAGGED}"
echo ""

# ── Per-cohort summary ────────────────────────────────────────────────────────
echo "── Per-cohort summary ──────────────────────────────────"
printf "%-25s %5s %8s %8s %8s %10s %10s %8s\n" \
    "Cohort" "N" "Mean%" "Min%" "Max%" "MeanJunc" "MeanAT/AC" "MeanBAM"
printf "%-25s %5s %8s %8s %8s %10s %10s %8s\n" \
    "------" "-" "-----" "----" "----" "--------" "---------" "-------"

for cohort in $(echo "${!COHORT_COUNT[@]}" | tr ' ' '\n' | sort); do
    n=${COHORT_COUNT[$cohort]}
    mean_map=$(echo "scale=1; ${COHORT_MAP_SUM[$cohort]} / ${n}" | bc)
    mean_junc=$(( ${COHORT_JUNC_SUM[$cohort]} / n ))
    mean_atac=$(( ${COHORT_ATAC_SUM[$cohort]} / n ))
    mean_bam=$(( ${COHORT_BAM_SUM[$cohort]} / n ))
    printf "%-25s %5d %7s%% %7s%% %7s%% %10d %10d %6dMB\n" \
        "${cohort}" "${n}" "${mean_map}" \
        "${COHORT_MAP_MIN[$cohort]}" "${COHORT_MAP_MAX[$cohort]}" \
        "${mean_junc}" "${mean_atac}" "${mean_bam}"
done

# ── Pass 1 vs Pass 2 comparison ───────────────────────────────────────────────
echo ""
echo "── Pass 1 vs Pass 2 comparison ────────────────────────"
if [[ -f "${PASS1_QC}" ]]; then
    echo "  Metric              Pass 1          Pass 2"
    echo "  ------              ------          ------"

    P1_MAP=$(awk -F'\t' 'NR>1 && $4!=""{sum+=$4; count++} \
        END{printf "%.1f", sum/count}' "${PASS1_QC}" 2>/dev/null)
    P2_MAP=$(awk -F'\t' 'NR>1 && $4!=""{sum+=$4; count++} \
        END{printf "%.1f", sum/count}' "${SUMMARY_FILE}" 2>/dev/null)

    P1_JUNC=$(awk -F'\t' 'NR>1 && $6!=""{sum+=$6; count++} \
        END{printf "%d", sum/count}' "${PASS1_QC}" 2>/dev/null)
    P2_JUNC=$(awk -F'\t' 'NR>1 && $5!=""{sum+=$5; count++} \
        END{printf "%d", sum/count}' "${SUMMARY_FILE}" 2>/dev/null)

    P1_ATAC=$(awk -F'\t' 'NR>1 && $6!=""{sum+=$6; count++} \
        END{printf "%d", sum/count}' "${PASS1_QC}" 2>/dev/null)
    P2_ATAC=$(awk -F'\t' 'NR>1 && $6!=""{sum+=$6; count++} \
        END{printf "%d", sum/count}' "${SUMMARY_FILE}" 2>/dev/null)

    printf "  %-20s %-15s %-15s\n" "Mean mapping %" "${P1_MAP}%" "${P2_MAP}%"
    printf "  %-20s %-15s %-15s\n" "Mean junctions" "${P1_JUNC}" "${P2_JUNC}"
    printf "  %-20s %-15s %-15s\n" "Mean AT/AC splices" "${P1_ATAC}" "${P2_ATAC}"
else
    echo "  Pass 1 QC summary not found — skipping comparison"
    echo "  Expected: ${PASS1_QC}"
fi

# ── Flagged samples ───────────────────────────────────────────────────────────
echo ""
echo "── Flagged samples (${FLAGGED}) ─────────────────────────"
if [[ "${FLAGGED}" -gt 0 ]]; then
    column -t "${FLAGGED_FILE}" | head -20
    [[ "${FLAGGED}" -gt 20 ]] && echo "  ... and $(( FLAGGED - 20 )) more — see ${FLAGGED_FILE}"
else
    echo "  None — all completed samples passed QC"
fi

# ── Missing samples ───────────────────────────────────────────────────────────
echo ""
echo "── Missing samples (${MISSING}) ─────────────────────────"
if [[ "${MISSING}" -gt 0 ]]; then
    echo ""
    N_LARGE=${#missing_large[@]}
    N_OTHER=${#missing_other[@]}

    if [[ "${N_LARGE}" -gt 0 ]]; then
        LARGE_SPEC=$(IFS=','; echo "${missing_large[*]}")
        echo "  Large FASTQs (${N_LARGE} samples, ~200GB scratch each):"
        echo "    sbatch --array=${LARGE_SPEC}%4 04_star_pass2.sh"
        echo ""
    fi

    if [[ "${N_OTHER}" -gt 0 ]]; then
        OTHER_SPEC=$(IFS=','; echo "${missing_other[*]}")
        echo "  Other cohorts (${N_OTHER} samples, ~50GB scratch each):"
        echo "    sbatch --array=${OTHER_SPEC}%15 04_star_pass2.sh"
        echo ""
    fi

    ALL_SPEC=$(IFS=','; echo "${missing_tasks[*]}")
    echo "  Or all together (conservative):"
    echo "    sbatch --array=${ALL_SPEC}%4 04_star_pass2.sh"
else
    echo "  None — all samples complete"
    echo ""
    echo "  Ready for next step: bash 05_stringtie_per_sample.sh"
fi

echo ""
echo "── Output files ────────────────────────────────────────"
echo "  QC summary : ${SUMMARY_FILE}"
echo "  Flagged    : ${FLAGGED_FILE}"
echo "  Missing    : ${MISSING_FILE}"
echo "  Report     : ${REPORT}"
echo ""
echo "Completed: $(date)"
echo "======================================="