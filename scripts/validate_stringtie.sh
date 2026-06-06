#!/bin/bash
# =============================================================================
# validate_stringtie.sh
# Validates per-sample and merged StringTie outputs before running rMATS.
#
# Usage: bash validate_stringtie.sh
# =============================================================================

STRINGTIE_DIR=/data1/abdelwao/maxim/splicing_pipeline/stringtie
MERGED_FINAL=${STRINGTIE_DIR}/merge/merged_final.gtf
MERGED_RAW=${STRINGTIE_DIR}/merge/merged.gtf
GFFCMP=${STRINGTIE_DIR}/merge/merged_compared.annotated.gtf
SUMMARY=${STRINGTIE_DIR}/merge/merge_summary.txt
MANIFEST=/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv

echo "================================================"
echo "1. PER-SAMPLE GTF COMPLETION"
echo "================================================"
total_samples=$(tail -n +2 "${MANIFEST}" | wc -l)
total_gtfs=$(find ${STRINGTIE_DIR}/gtf/ -name "*.gtf" -size +0c | wc -l)
empty_gtfs=$(find ${STRINGTIE_DIR}/gtf/ -name "*.gtf" -empty | wc -l)
echo "  Samples in manifest : ${total_samples}"
echo "  GTFs produced       : ${total_gtfs}"
echo "  Empty GTFs          : ${empty_gtfs}"

# List any missing samples
missing=0
while IFS=$'\t' read -r sample_id cohort rest; do
    if [[ ! -f "${STRINGTIE_DIR}/gtf/${sample_id}.gtf" ]] || \
       [[ ! -s "${STRINGTIE_DIR}/gtf/${sample_id}.gtf" ]]; then
        echo "  MISSING: ${sample_id} (${cohort})"
        missing=$((missing + 1))
    fi
done < <(tail -n +2 "${MANIFEST}")
[[ "$missing" -eq 0 ]] && echo "  All samples complete"
echo ""

echo "================================================"
echo "2. MERGED GTF SANITY"
echo "================================================"
echo "  Transcript counts:"
echo "    merged_raw   : $(awk '$3=="transcript"' ${MERGED_RAW} | wc -l)"
echo "    merged_final : $(awk '$3=="transcript"' ${MERGED_FINAL} | wc -l)"
echo ""
echo "  Chromosomes in merged_final (should all be chr*):"
cut -f1 ${MERGED_FINAL} | grep -v "^#" | sort -u
echo ""
echo "  Strand distribution (should have both + and -):"
cut -f7 ${MERGED_FINAL} | grep -v "^#" | sort | uniq -c
echo ""

echo "================================================"
echo "3. NOVEL TRANSCRIPT CONTENT"
echo "================================================"
n_enst=$(awk '$3=="transcript"' ${MERGED_FINAL} | grep -c 'transcript_id "ENST' || true)
n_mstrg=$(awk '$3=="transcript"' ${MERGED_FINAL} | grep -c 'transcript_id "MSTRG' || true)
echo "  GENCODE transcripts (ENST*) : ${n_enst}"
echo "  Novel transcripts (MSTRG*)  : ${n_mstrg}"
echo "  Expected novel range        : 20,000 - 60,000"
[[ "$n_mstrg" -lt 20000 ]] && echo "  WARNING: Novel count low — assembly filters may be too strict or samples missing from merge"
[[ "$n_mstrg" -gt 200000 ]] && echo "  WARNING: Novel count high — low-coverage sample noise may be inflating the GTF"
echo ""

echo "================================================"
echo "4. GFFCOMPARE CLASS CODE SUMMARY"
echo "================================================"
if [[ -f "${SUMMARY}" ]]; then
    cat "${SUMMARY}"
else
    echo "  WARNING: ${SUMMARY} not found — was step 3 (Python filtering) run?"
fi
echo ""

echo "================================================"
echo "5. KNOWN GENE SPOT CHECK"
echo "================================================"
for gene in ACTB SF3B1 SRSF2 U2AF1 ZRSR2 IL17RC SLC3A2 TFR2; do
    n=$(grep -c "gene_name \"${gene}\"" ${MERGED_FINAL} || true)
    status="OK"
    [[ "$n" -eq 0 ]] && status="MISSING — check GTF"
    echo "  ${gene}: ${n} features  [${status}]"
done
echo ""

echo "================================================"
echo "6. STRANDEDNESS CALL DISTRIBUTION ACROSS SAMPLES"
echo "================================================"
if ls ${STRINGTIE_DIR}/strand/*.strand.txt 1>/dev/null 2>&1; then
    cut -f3 ${STRINGTIE_DIR}/strand/*.strand.txt | sort | uniq -c | sed 's/^/  /'
    echo ""
    # Flag if unstranded is unexpectedly dominant
    n_unstranded=$(cut -f3 ${STRINGTIE_DIR}/strand/*.strand.txt | grep -c "unstranded" || true)
    n_total=$(cut -f3 ${STRINGTIE_DIR}/strand/*.strand.txt | wc -l)
    pct_unstranded=$(echo "scale=1; ${n_unstranded} * 100 / ${n_total}" | bc -l 2>/dev/null || echo "?")
    echo "  Unstranded fraction: ${pct_unstranded}%"
    [[ "$n_unstranded" -gt $(( n_total / 2 )) ]] && \
        echo "  WARNING: >50% samples called unstranded — check strandedness logic in extract_bam_info.sh"
else
    echo "  WARNING: No strand files found in ${STRINGTIE_DIR}/strand/"
fi
echo ""

echo "================================================"
echo "7. NOVEL MULTI-EXON TRANSCRIPTS (rMATS targets)"
echo "================================================"
if [[ -f "${GFFCMP}" ]]; then
    n_j=$(grep 'class_code "j"' ${GFFCMP} | awk '$3=="transcript"' | wc -l)
    n_i=$(grep 'class_code "i"' ${GFFCMP} | awk '$3=="transcript"' | wc -l)
    n_u=$(grep 'class_code "u"' ${GFFCMP} | awk '$3=="transcript"' | wc -l)
    echo "  class j (novel junction, known gene)     : ${n_j}  ← primary rMATS targets"
    echo "  class i (intronic / retained intron)     : ${n_i}  ← ZRSR2 targets"
    echo "  class u (intergenic — mostly excluded)   : ${n_u}"
    [[ "$n_j" -eq 0 ]] && echo "  WARNING: No class j transcripts — something is wrong with GFFCompare or the merge"
else
    echo "  WARNING: ${GFFCMP} not found — was GFFCompare (step 2) run?"
fi
echo ""

echo "================================================"
echo "8. GTF READY FOR RMATS"
echo "================================================"
if [[ -f "${MERGED_FINAL}" && -s "${MERGED_FINAL}" ]]; then
    echo "  merged_final.gtf exists and is non-empty"
    echo "  Path to use as --gtf in rMATS:"
    echo "    ${MERGED_FINAL}"
else
    echo "  ERROR: ${MERGED_FINAL} missing or empty — do not proceed with rMATS"
fi
echo ""
echo "Done — $(date)"