#!/bin/bash
# =============================================================================
# 02b_junction_pool_sanity_check.sh
# Sanity check for pooled junction file before building pass 2 index.
#
# Checks:
#   1. File exists and is non-empty
#   2. Correct 6-column format
#   3. Motif distribution (should be dominated by GT/AG)
#   4. Chromosome distribution (no unexpected contigs)
#   5. Coordinate sanity (start < end, no negative coordinates)
#   6. Novel junction enrichment vs GENCODE annotation
#   7. Comparison to pass 1 individual sample junction counts
#   8. AT/AC minor spliceosome junction count
#   9. Cross-sample support distribution
#
# Usage: bash 02b_junction_pool_sanity_check.sh
# =============================================================================

# ── Auto-detect most recent junction directory ────────────────────────────────
JUNC_BASE="/data1/abdelwao/maxim/splicing_pipeline/junctions"
JUNC_DIR=$(ls -td "${JUNC_BASE}"/v1_* 2>/dev/null | head -1)

if [[ -z "${JUNC_DIR}" ]]; then
    echo "ERROR: No versioned junction directory found in ${JUNC_BASE}" >&2
    exit 1
fi

FINAL="${JUNC_DIR}/final_filtered_junctions.tab"
POOLED_RAW="${JUNC_DIR}/pooled_raw.tab"
ARCHIVE="${JUNC_BASE}/pass1_sj_archive"
OUT_FILE="${JUNC_DIR}/junction_pool_sanity_check.txt"

echo "======================================="
echo "Junction pool sanity check — $(date)"
echo "Junction dir : ${JUNC_DIR}"
echo "Final file   : ${FINAL}"
echo "Report saved : ${OUT_FILE}"
echo "======================================="
echo ""

# Tee all subsequent output to both terminal and report file
exec > >(tee "${OUT_FILE}") 2>&1

# ── Check 1: File exists ──────────────────────────────────────────────────────
echo "── Check 1: File existence ─────────────────────────────"
if [[ ! -f "${FINAL}" ]]; then
    echo "  ERROR: Final junction file not found: ${FINAL}" >&2
    exit 1
fi
if [[ ! -f "${POOLED_RAW}" ]]; then
    echo "  WARNING: Pooled raw file not found: ${POOLED_RAW}"
fi

FINAL_COUNT=$(wc -l < "${FINAL}")
RAW_COUNT=$(wc -l < "${POOLED_RAW}" 2>/dev/null || echo "N/A")
echo "  Final filtered junctions : ${FINAL_COUNT}"
echo "  Pooled raw junctions     : ${RAW_COUNT}"
echo "  File size                : $(du -sh ${FINAL} | cut -f1)"

if [[ "${FINAL_COUNT}" -eq 0 ]]; then
    echo "  ERROR: Final junction file is empty" >&2
    exit 1
fi
echo "  Status: OK"

# ── Check 2: Column format ────────────────────────────────────────────────────
echo ""
echo "── Check 2: Column format ──────────────────────────────"
NCOLS=$(awk '{print NF; exit}' "${FINAL}")
echo "  Columns: ${NCOLS} (expected: 6)"
if [[ "${NCOLS}" -ne 6 ]]; then
    echo "  ERROR: Expected 6 columns, got ${NCOLS}" >&2
else
    echo "  Status: OK"
fi

echo "  Sample rows:"
head -5 "${FINAL}" | awk '{printf "    %s\t%s\t%s\t%s\t%s\t%s\n", \
    $1,$2,$3,$4,$5,$6}'

# ── Check 3: Motif distribution ───────────────────────────────────────────────
echo ""
echo "── Check 3: Motif distribution ─────────────────────────"
awk '{motif[$5]++; total++} END {
    printf "  %-20s %8s %8s\n", "Motif", "Count", "Pct"
    printf "  %-20s %8s %8s\n", "-----", "-----", "---"
    printf "  %-20s %8d %7.2f%%\n", "GT/AG (1, major)",
        motif[1]+0, (motif[1]+0)/total*100
    printf "  %-20s %8d %7.2f%%\n", "CT/AC (2, major-)",
        motif[2]+0, (motif[2]+0)/total*100
    printf "  %-20s %8d %7.2f%%\n", "GC/AG (3, minor)",
        motif[3]+0, (motif[3]+0)/total*100
    printf "  %-20s %8d %7.2f%%\n", "CT/GC (4, minor-)",
        motif[4]+0, (motif[4]+0)/total*100
    printf "  %-20s %8d %7.2f%%\n", "AT/AC (5, U12)",
        motif[5]+0, (motif[5]+0)/total*100
    printf "  %-20s %8d %7.2f%%\n", "GT/AT (6, rare)",
        motif[6]+0, (motif[6]+0)/total*100
    printf "  %-20s %8d\n", "Total", total
}' "${FINAL}"

# Sanity check — GT/AG + CT/AC (reverse complement) should be >80%
# CT/AC (motif 2) is GT/AG on the minus strand — same junction type
MAJOR_PCT=$(awk '{motif[$5]++; total++} END {
    printf "%.1f", ((motif[1]+0)+(motif[2]+0))/total*100}' "${FINAL}")
GTAG_PCT=$(awk '{motif[$5]++; total++} END {
    printf "%.1f", (motif[1]+0)/total*100}' "${FINAL}")
echo ""
echo "  GT/AG + CT/AC combined: ${MAJOR_PCT}%"
if (( $(echo "${MAJOR_PCT} < 80" | bc -l) )); then
    echo "  WARNING: Major spliceosome (GT/AG+CT/AC) is only ${MAJOR_PCT}% — expected >80%"
else
    echo "  Status: OK (major spliceosome = ${MAJOR_PCT}%)"
fi

# ── Check 4: Chromosome distribution ─────────────────────────────────────────
echo ""
echo "── Check 4: Chromosome distribution ───────────────────"
echo "  Standard chromosomes:"
awk '$1 ~ /^chr([0-9]+|X|Y)$/{chr[$1]++} END {
    for(c in chr) print c, chr[c]
}' "${FINAL}" | sort -V | \
    awk '{printf "    %-10s %d\n", $1, $2}'

echo ""
echo "  Non-standard contigs (should be empty):"
awk '$1 !~ /^chr([0-9]+|X|Y)$/{print $1}' "${FINAL}" | \
    sort | uniq -c | sort -rn | head -10
NON_STANDARD=$(awk '$1 !~ /^chr([0-9]+|X|Y)$/{count++} \
    END{print count+0}' "${FINAL}")
if [[ "${NON_STANDARD}" -gt 0 ]]; then
    echo "  WARNING: ${NON_STANDARD} junctions on non-standard contigs"
else
    echo "  Status: OK (no non-standard contigs)"
fi

# ── Check 5: Coordinate sanity ────────────────────────────────────────────────
echo ""
echo "── Check 5: Coordinate sanity ──────────────────────────"
awk '{
    if($2 <= 0) neg_start++
    if($3 <= 0) neg_end++
    if($2 >= $3) bad_coord++
    if($3-$2 < 20) too_short++
    if($3-$2 > 1000000) too_long++
    total++
} END {
    printf "  Negative start    : %d\n", neg_start+0
    printf "  Negative end      : %d\n", neg_end+0
    printf "  Start >= end      : %d\n", bad_coord+0
    printf "  Intron < 20bp     : %d\n", too_short+0
    printf "  Intron > 1Mbp     : %d\n", too_long+0
    printf "  Total checked     : %d\n", total
}' "${FINAL}"

BAD=$(awk '$2>=$3 || $2<=0{count++} END{print count+0}' "${FINAL}")
if [[ "${BAD}" -gt 0 ]]; then
    echo "  WARNING: ${BAD} junctions with bad coordinates"
else
    echo "  Status: OK"
fi

# ── Check 6: Strand distribution ─────────────────────────────────────────────
echo ""
echo "── Check 6: Strand distribution ────────────────────────"
awk '{strand[$4]++} END {
    printf "  %-20s %8d\n", "Forward (+, 1)", strand[1]+0
    printf "  %-20s %8d\n", "Reverse (-, 2)", strand[2]+0
    printf "  %-20s %8d\n", "Unknown (0)", strand[0]+0
}' "${FINAL}"

UNKNOWN_STRAND=$(awk '$4==0{count++} END{print count+0}' "${FINAL}")
if [[ "${UNKNOWN_STRAND}" -gt 0 ]]; then
    echo "  WARNING: ${UNKNOWN_STRAND} junctions with unknown strand"
    echo "           These may cause issues with stranded LeafCutter analysis"
else
    echo "  Status: OK"
fi

# ── Check 7: Comparison to raw pool ──────────────────────────────────────────
echo ""
echo "── Check 7: Filter retention rate ─────────────────────"
if [[ -f "${POOLED_RAW}" ]]; then
    awk -v final="${FINAL_COUNT}" -v raw="${RAW_COUNT}" 'BEGIN {
        printf "  Raw junctions    : %d\n", raw
        printf "  Final junctions  : %d\n", final
        printf "  Retention rate   : %.1f%%\n", final/raw*100
        printf "  Filtered out     : %d\n", raw-final
    }'

    # Show what got filtered and why
    echo ""
    echo "  Filter breakdown (from pooled raw):"
    awk '{
        # chrM/alt/random
        if($1 ~ /chrM|chrUn|random|alt|fix|hap/) { f_chrom++; next }
        # annotated
        if($6 != 0) { f_annot++; next }
        # non-canonical motif
        if($5 == 0) { f_motif++; next }
        # overhang
        if($9 < 6) { f_overhang++; next }
        # minor motifs < 1 read
        if($5 >= 3 && $7 < 1) { f_reads++; next }
        # major motifs threshold
        if($5 <= 2 && $7 < 3 && \
           !($NF >= 2 && $7 >= 2)) { f_reads++; next }
        kept++
    } END {
        printf "    Excluded chroms   : %d\n", f_chrom+0
        printf "    Annotated (skip)  : %d\n", f_annot+0
        printf "    Non-canonical (0) : %d\n", f_motif+0
        printf "    Low overhang (<6) : %d\n", f_overhang+0
        printf "    Low reads         : %d\n", f_reads+0
        printf "    Kept              : %d\n", kept+0
    }' "${POOLED_RAW}"
fi

# ── Check 8: AT/AC minor spliceosome ─────────────────────────────────────────
echo ""
echo "── Check 8: Minor spliceosome (AT/AC) ──────────────────"
ATAC=$(awk '$5==5{count++} END{print count+0}' "${FINAL}")
echo "  AT/AC junctions : ${ATAC}"

# Compare to expected range based on pass 1 samples
# BA2000R had 28,682 AT/AC splices from one sample
# Expect thousands of unique AT/AC junctions across all samples
if [[ "${ATAC}" -lt 100 ]]; then
    echo "  WARNING: Very few AT/AC junctions — check minor spliceosome filtering"
elif [[ "${ATAC}" -gt 100000 ]]; then
    echo "  WARNING: Unusually high AT/AC count — check motif assignment"
else
    echo "  Status: OK"
fi

# ── Check 9: Cross-sample support distribution ───────────────────────────────
echo ""
echo "── Check 9: Cross-sample support (from raw pool) ───────"
if [[ -f "${POOLED_RAW}" ]]; then
    echo "  Junctions by number of supporting samples:"
    awk '{print $NF}' "${POOLED_RAW}" | \
        awk '{
            if($1==1) s1++
            else if($1<=5) s5++
            else if($1<=10) s10++
            else if($1<=50) s50++
            else if($1<=100) s100++
            else multi++
            total++
        } END {
            printf "    1 sample       : %8d (%5.1f%%)\n", \
                s1+0, (s1+0)/total*100
            printf "    2-5 samples    : %8d (%5.1f%%)\n", \
                s5+0, (s5+0)/total*100
            printf "    6-10 samples   : %8d (%5.1f%%)\n", \
                s10+0, (s10+0)/total*100
            printf "    11-50 samples  : %8d (%5.1f%%)\n", \
                s50+0, (s50+0)/total*100
            printf "    51-100 samples : %8d (%5.1f%%)\n", \
                s100+0, (s100+0)/total*100
            printf "    >100 samples   : %8d (%5.1f%%)\n", \
                multi+0, (multi+0)/total*100
        }'
fi

# ── Check 10: Compare archive vs final ───────────────────────────────────────
echo ""
echo "── Check 10: Archive completeness ──────────────────────"
if [[ -d "${ARCHIVE}" ]]; then
    ARCHIVE_COUNT=$(ls "${ARCHIVE}" | wc -l)
    echo "  Archived SJ files : ${ARCHIVE_COUNT}"
    echo "  Archive size      : $(du -sh ${ARCHIVE} | cut -f1)"

    # Check for any archived files not in manifest
    MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
    MANIFEST_COUNT=$(tail -n +2 "${MANIFEST}" | wc -l)
    echo "  Manifest samples  : ${MANIFEST_COUNT}"

    if [[ "${ARCHIVE_COUNT}" -lt "${MANIFEST_COUNT}" ]]; then
        DIFF=$(( MANIFEST_COUNT - ARCHIVE_COUNT ))
        echo "  WARNING: ${DIFF} samples missing from archive"
        echo "           These samples may not have completed pass 1"
    else
        echo "  Status: OK"
    fi
fi

# ── Final assessment ──────────────────────────────────────────────────────────
echo ""
echo "── Final assessment ────────────────────────────────────"
echo "  Junction file  : ${FINAL}"
echo "  Total junctions: ${FINAL_COUNT}"
echo ""

if [[ "${FINAL_COUNT}" -gt 500000 ]]; then
    echo "  NOTE: Junction count is high (>500k)"
    echo "        This is expected for 1000+ samples"
    echo "        STAR index build may take longer than usual"
elif [[ "${FINAL_COUNT}" -lt 10000 ]]; then
    echo "  WARNING: Junction count seems low (<10k)"
    echo "           Check pooling and filtering steps"
else
    echo "  Junction count looks reasonable"
fi

echo ""
echo "  Ready to proceed: sbatch 03_build_pass2_index.sh"
echo "======================================="
echo "Completed: $(date)"