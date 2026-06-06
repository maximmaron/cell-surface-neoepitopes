#!/bin/bash
# =============================================================================
# find_missing_pass2.sh
# Find samples missing or incomplete pass 2 BAM outputs.
#
# Checks per sample:
#   - BAM exists AND size > 0
#   - BAI exists AND size > 0
#   - Log.final.out exists AND size > 0
#   - BAM is not suspiciously small (< 100MB for a real RNA-seq BAM)
#
# Leucegene samples are separated into their own submission command
# with lower concurrency (%4) due to large FASTQ sizes (~200GB scratch).
#
# Usage: bash find_missing_pass2.sh
# =============================================================================

MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
MIN_BAM_SIZE=104857600  # 100MB minimum

missing_tasks=()
missing_large=()
missing_other=()
complete=0
empty_bam=0
task_num=1

while IFS=$'\t' read -r sample_id cohort tier disease seq_type fq1 fq2; do
    out_dir=$(dirname "${fq1}")
    bam="${out_dir}/${sample_id}_Aligned.sortedByCoord.out.bam"
    bai="${out_dir}/${sample_id}_Aligned.sortedByCoord.out.bam.bai"
    log="${out_dir}/${sample_id}_Log.final.out"

    is_missing=0

    # Check existence
    if [[ ! -f "${bam}" || ! -f "${bai}" || ! -f "${log}" ]]; then
        is_missing=1
    fi

    # Check BAM is non-empty
    if [[ "${is_missing}" -eq 0 ]]; then
        bam_size=$(stat -c%s "${bam}" 2>/dev/null || echo 0)
        if [[ "${bam_size}" -eq 0 ]]; then
            echo "EMPTY BAM: ${sample_id} (0 bytes) — removing and will rerun"
            rm -f "${bam}" "${bai}"
            ((empty_bam++))
            is_missing=1
        elif [[ "${bam_size}" -lt "${MIN_BAM_SIZE}" ]]; then
            bam_mb=$(( bam_size / 1024 / 1024 ))
            echo "SMALL BAM: ${sample_id} (${bam_mb}MB) — flagged for rerun"
            is_missing=1
        fi
    fi

    # Check BAI is non-empty
    if [[ "${is_missing}" -eq 0 ]]; then
        bai_size=$(stat -c%s "${bai}" 2>/dev/null || echo 0)
        if [[ "${bai_size}" -eq 0 ]]; then
            echo "EMPTY BAI: ${sample_id} — removing and will rerun"
            rm -f "${bai}"
            is_missing=1
        fi
    fi

    if [[ "${is_missing}" -eq 1 ]]; then
        missing_tasks+=("${task_num}")
        # Determine concurrency tier by compressed FASTQ size
        # Large: >3GB compressed (~30GB+ uncompressed) needs low concurrency
        fq_size=$(stat -c%s "${fq1}" 2>/dev/null || echo 0)
        fq_gb=$(echo "scale=1; ${fq_size}/1024/1024/1024" | bc)
        if (( $(echo "${fq_gb} > 3" | bc -l) )); then
            missing_large+=("${task_num}")
        else
            missing_other+=("${task_num}")
        fi
    else
        ((complete++))
    fi

    ((task_num++))

done < <(tail -n +2 "${MANIFEST}")

TOTAL=$(( task_num - 1 ))
MISSING=${#missing_tasks[@]}
N_LARGE=${#missing_large[@]}
N_OTHER=${#missing_other[@]}

echo "======================================="
echo "Pass 2 completion check — $(date)"
echo "======================================="
echo "Total samples     : ${TOTAL}"
echo "Complete          : ${complete}"
echo "Missing/incomplete: ${MISSING}"
echo "  Leucegene       : ${N_LARGE}"
echo "  Other cohorts   : ${N_OTHER}"
echo "  Empty BAMs removed: ${empty_bam}"
echo ""

if [[ "${MISSING}" -eq 0 ]]; then
    echo "All samples complete — ready for StringTie"
    exit 0
fi

# ── Leucegene submission ──────────────────────────────────────────────────────
if [[ "${N_LARGE}" -gt 0 ]]; then
    LARGE_SPEC=$(IFS=','; echo "${missing_large[*]}")
    echo "── Large FASTQs (${N_LARGE} samples >3GB compressed, ~100-200GB scratch) ────"
    echo "   Tasks : ${LARGE_SPEC}" | fold -w 72 -s | sed '2,$s/^/           /'
    echo "   Submit: sbatch --array=${LARGE_SPEC}%4 04_star_pass2.sh"
    echo "   (1TB scratch / 200GB = max 5 concurrent, using %4 for safety)"
    echo ""
fi

# ── Other cohorts submission ──────────────────────────────────────────────────
if [[ "${N_OTHER}" -gt 0 ]]; then
    OTHER_SPEC=$(IFS=','; echo "${missing_other[*]}")
    echo "── Other cohorts (${N_OTHER} samples, ~50GB scratch each) ──"
    echo "   Tasks : ${OTHER_SPEC}" | fold -w 72 -s | sed '2,$s/^/           /'
    echo "   Submit: sbatch --array=${OTHER_SPEC}%15 04_star_pass2.sh"
    echo "   (1TB scratch / 50GB = max 20 concurrent, using %15 for safety)"
    echo ""
fi

# ── Combined if preferred ─────────────────────────────────────────────────────
ALL_SPEC=$(IFS=','; echo "${missing_tasks[*]}")
echo "── Or submit all together at conservative concurrency ───"
echo "   sbatch --array=${ALL_SPEC}%4 04_star_pass2.sh"
echo "======================================="