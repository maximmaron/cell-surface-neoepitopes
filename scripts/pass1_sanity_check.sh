#!/bin/bash
# =============================================================================
# 01b_pass1_sanity_check.sh
# Sanity check for pass 1 completion across all samples in manifest.
#
# Checks per sample:
#   1. SJ.out.tab exists and is non-empty
#   2. Log.final.out exists
#   3. No FATAL errors in Log.out
#   4. Mapping rate above cohort-specific threshold
#   5. Junction count above minimum threshold
#   6. fastp QC JSON exists
#
# Outputs:
#   - Summary counts (complete/missing/flagged)
#   - List of missing samples with task numbers for resubmission
#   - List of flagged samples with specific issues
#   - Per-cohort mapping rate summary
#
# Usage: bash pass1_sanity_check.sh
# =============================================================================

MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
QC_DIR="/data1/abdelwao/maxim/splicing_pipeline/qc/fastp"
LOG_DIR="logs"
OUT_DIR="/data1/abdelwao/maxim/splicing_pipeline/qc"

mkdir -p "${OUT_DIR}"

MISSING_FILE="${OUT_DIR}/pass1_missing_samples.tsv"
FLAGGED_FILE="${OUT_DIR}/pass1_flagged_samples.tsv"
SUMMARY_FILE="${OUT_DIR}/pass1_qc_summary.tsv"

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
MIN_JUNCTIONS["beat"]=100000
MIN_JUNCTIONS["tcga"]=50000
MIN_JUNCTIONS["leucegene_SRSF2"]=100000
MIN_JUNCTIONS["leucegene_SF3B1"]=100000
MIN_JUNCTIONS["leucegene_U2AF1"]=100000
MIN_JUNCTIONS["leucegene_ZRSR2"]=100000
MIN_JUNCTIONS["pellagatti_SF3B1"]=50000
MIN_JUNCTIONS["pellagatti_SRSF2"]=50000
MIN_JUNCTIONS["pellagatti_U2AF1"]=50000
MIN_JUNCTIONS["pellagatti_WT"]=50000
MIN_JUNCTIONS["madan_2015_ZRSR2"]=50000
MIN_JUNCTIONS["madan_2015_WT"]=50000
MIN_JUNCTIONS["maiga_normal"]=50000
MIN_JUNCTIONS["madan_normal"]=50000
MIN_JUNCTIONS["hpa_2014"]=50000
MIN_JUNCTIONS["bodymap2"]=50000

echo "======================================="
echo "Pass 1 sanity check — $(date)"
echo "Manifest: ${MANIFEST}"
echo "======================================="
echo ""

# ── Initialize output files ───────────────────────────────────────────────────
echo -e "task_num\tsample_id\tcohort\treason" > "${MISSING_FILE}"
echo -e "sample_id\tcohort\tmapping_pct\tjunctions\tflags" > "${FLAGGED_FILE}"
echo -e "sample_id\tcohort\ttier\tmapping_pct\tjunctions\tat_ac_splices\tgt_ag_splices\tinput_reads" \
    > "${SUMMARY_FILE}"

# ── Counters ──────────────────────────────────────────────────────────────────
TOTAL=0
COMPLETE=0
MISSING=0
FLAGGED=0
task_num=1

# ── Per-cohort accumulators ───────────────────────────────────────────────────
declare -A COHORT_COUNT
declare -A COHORT_MAP_SUM
declare -A COHORT_MAP_MIN
declare -A COHORT_MAP_MAX
declare -A COHORT_JUNC_SUM

# ── Main loop ─────────────────────────────────────────────────────────────────
while IFS=$'\t' read -r sample_id cohort tier disease seq_type fq1 fq2; do
    ((TOTAL++))
    out_dir=$(dirname "${fq1}")
    sj="${out_dir}/${sample_id}_SJ.out.tab"
    log="${out_dir}/${sample_id}_Log.final.out"
    logout="${out_dir}/${sample_id}_Log.out"
    fastp_json="${QC_DIR}/${sample_id}_fastp.json"

    # ── Check 1: SJ.out.tab exists ────────────────────────────────────────────
    if [[ ! -f "${sj}" ]]; then
        echo -e "${task_num}\t${sample_id}\t${cohort}\tMissing SJ.out.tab" \
            >> "${MISSING_FILE}"
        ((MISSING++))
        ((task_num++))
        continue
    fi

    # ── Check 2: Log.final.out exists ────────────────────────────────────────
    if [[ ! -f "${log}" ]]; then
        echo -e "${task_num}\t${sample_id}\t${cohort}\tMissing Log.final.out" \
            >> "${MISSING_FILE}"
        ((MISSING++))
        ((task_num++))
        continue
    fi

    # ── Check 3: SJ.out.tab non-empty ────────────────────────────────────────
    sj_count=$(wc -l < "${sj}")
    if [[ "${sj_count}" -eq 0 ]]; then
        echo -e "${task_num}\t${sample_id}\t${cohort}\tEmpty SJ.out.tab" \
            >> "${MISSING_FILE}"
        ((MISSING++))
        ((task_num++))
        continue
    fi

    # ── Check 4: No FATAL errors in Log.out ──────────────────────────────────
    if [[ -f "${logout}" ]]; then
        if grep -q "FATAL\|EXITING because of fatal" "${logout}" 2>/dev/null; then
            echo -e "${task_num}\t${sample_id}\t${cohort}\tFATAL error in Log.out" \
                >> "${MISSING_FILE}"
            ((MISSING++))
            ((task_num++))
            continue
        fi
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

    # ── Write to summary ──────────────────────────────────────────────────────
    echo -e "${sample_id}\t${cohort}\t${tier}\t${mapping_pct}\t${sj_count}\t${at_ac}\t${gt_ag}\t${input_reads}" \
        >> "${SUMMARY_FILE}"

    # ── Check 5: Mapping rate ─────────────────────────────────────────────────
    flags=""
    min_map=${MIN_MAPPING[$cohort]:-75}
    if (( $(echo "${mapping_pct} < ${min_map}" | bc -l) )); then
        flags="${flags}LOW_MAPPING(${mapping_pct}%<${min_map}%);"
    fi

    # ── Check 6: Junction count ───────────────────────────────────────────────
    min_junc=${MIN_JUNCTIONS[$cohort]:-50000}
    if [[ "${sj_count}" -lt "${min_junc}" ]]; then
        flags="${flags}LOW_JUNCTIONS(${sj_count}<${min_junc});"
    fi

    # ── Check 7: fastp JSON exists ────────────────────────────────────────────
    if [[ ! -f "${fastp_json}" ]]; then
        flags="${flags}MISSING_FASTP_QC;"
    fi

    # ── Check 8: Suspiciously low input reads ─────────────────────────────────
    if [[ "${input_reads}" -lt 1000000 ]]; then
        flags="${flags}LOW_INPUT_READS(${input_reads});"
    fi

    # ── Flag if issues found ──────────────────────────────────────────────────
    if [[ -n "${flags}" ]]; then
        echo -e "${sample_id}\t${cohort}\t${mapping_pct}\t${sj_count}\t${flags}" \
            >> "${FLAGGED_FILE}"
        ((FLAGGED++))
    fi

    # ── Accumulate cohort stats ───────────────────────────────────────────────
    COHORT_COUNT[$cohort]=$(( ${COHORT_COUNT[$cohort]:-0} + 1 ))
    COHORT_MAP_SUM[$cohort]=$(echo "${COHORT_MAP_SUM[$cohort]:-0} + ${mapping_pct}" | bc)
    COHORT_JUNC_SUM[$cohort]=$(( ${COHORT_JUNC_SUM[$cohort]:-0} + sj_count ))

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

# ── Print overall summary ─────────────────────────────────────────────────────
echo "── Overall ─────────────────────────────────────────────"
echo "  Total samples : ${TOTAL}"
echo "  Complete      : ${COMPLETE}"
echo "  Missing       : ${MISSING}"
echo "  Flagged       : ${FLAGGED}"
echo ""

# ── Print per-cohort summary ──────────────────────────────────────────────────
echo "── Per-cohort mapping rate summary ────────────────────"
printf "%-25s %6s %8s %8s %8s %10s\n" \
    "Cohort" "N" "Mean%" "Min%" "Max%" "MeanJunc"
printf "%-25s %6s %8s %8s %8s %10s\n" \
    "------" "-" "-----" "----" "----" "--------"

for cohort in $(echo "${!COHORT_COUNT[@]}" | tr ' ' '\n' | sort); do
    n=${COHORT_COUNT[$cohort]}
    mean_map=$(echo "scale=1; ${COHORT_MAP_SUM[$cohort]} / ${n}" | bc)
    mean_junc=$(( ${COHORT_JUNC_SUM[$cohort]} / n ))
    printf "%-25s %6d %7s%% %7s%% %7s%% %10d\n" \
        "${cohort}" "${n}" "${mean_map}" \
        "${COHORT_MAP_MIN[$cohort]}" "${COHORT_MAP_MAX[$cohort]}" \
        "${mean_junc}"
done

# ── Print missing samples with task numbers ───────────────────────────────────
echo ""
echo "── Missing samples (${MISSING}) ─────────────────────────"
if [[ "${MISSING}" -gt 0 ]]; then
    cat "${MISSING_FILE}"
    echo ""

    # Generate resubmission command
    MISSING_TASKS=$(awk -F'\t' 'NR>1{print $1}' "${MISSING_FILE}" | \
        tr '\n' ',' | sed 's/,$//')
    echo "── Resubmission command ────────────────────────────────"
    echo "  sbatch --array=${MISSING_TASKS}%20 01_star_pass1.sh"
else
    echo "  None — all samples complete"
fi

# ── Print flagged samples ─────────────────────────────────────────────────────
echo ""
echo "── Flagged samples (${FLAGGED}) ─────────────────────────"
if [[ "${FLAGGED}" -gt 0 ]]; then
    column -t "${FLAGGED_FILE}"
else
    echo "  None — all completed samples passed QC thresholds"
fi

echo ""
echo "── Output files ────────────────────────────────────────"
echo "  QC summary : ${SUMMARY_FILE}"
echo "  Missing    : ${MISSING_FILE}"
echo "  Flagged    : ${FLAGGED_FILE}"
echo ""
echo "Completed: $(date)"
echo "======================================="