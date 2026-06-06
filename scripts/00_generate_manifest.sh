#!/bin/bash
# =============================================================================
# 00_generate_manifest.sh
# Generate Stage 1 sample manifest from FASTQs de novo.
#
# Scans all cohort FASTQ directories, detects paired/single-end from naming,
# assigns cohort/tier/disease from directory path.
#
# FASTQ naming conventions:
#   _R1.fastq.gz / _R2.fastq.gz  : BEAT AML, TCGA
#   _1.fastq.gz  / _2.fastq.gz   : all other cohorts
#
# Output columns:
#   sample_id | cohort | tier | disease | seq_type | fastq_R1 | fastq_R2
#
# Tier assignments:
#   1 = MDS/AML
#   2 = Normal hematopoietic
#   3 = Normal other tissue
#
# Disease assignments:
#   MDS_AML              : all AML/MDS cohorts
#   normal_hematopoietic : Maiga, Madan normal
#   normal_tissue        : HPA, Bodymap2
#
# Usage: bash 00_generate_manifest.sh
# =============================================================================

BASE="/data1/abdelwao/shared/splicing_analysis/hg38"
OUT_DIR="/data1/abdelwao/shared/splicing_analysis/splicing_pipeline/metadata"
mkdir -p "${OUT_DIR}"
MANIFEST="${OUT_DIR}/sample_manifest.tsv"

echo "Generating sample manifest..."
echo ""

# ── Header ───────────────────────────────────────────────────────────────────
echo -e "sample_id\tcohort\ttier\tdisease\tseq_type\tfastq_R1\tfastq_R2" \
    > "${MANIFEST}"

# ── Helper: write one sample row ──────────────────────────────────────────────
add_sample() {
    local sample_id="$1"
    local cohort="$2"
    local tier="$3"
    local disease="$4"
    local fq1="$5"
    local fq2="${6:-}"
    local seq_type="single"
    [[ -n "${fq2}" && -f "${fq2}" ]] && seq_type="paired"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sample_id}" "${cohort}" "${tier}" "${disease}" \
        "${seq_type}" "${fq1}" "${fq2}" \
        >> "${MANIFEST}"
}

# ── Helper: process directory with _1/_2 naming ──────────────────────────────
process_1_2() {
    local dir="$1" cohort="$2" tier="$3" disease="$4"
    local count=0
    while IFS= read -r -d '' fq1; do
        fq2="${fq1/_1.fastq.gz/_2.fastq.gz}"
        sample=$(basename "${fq1}" _1.fastq.gz)
        add_sample "${sample}" "${cohort}" "${tier}" "${disease}" "${fq1}" "${fq2}"
        ((count++))
    done < <(find "${dir}" -maxdepth 1 -name "*_1.fastq.gz" -type f -print0 | sort -z)
    echo "  ${cohort}: ${count} samples"
}

# ── Helper: process directory with _R1/_R2 naming ────────────────────────────
process_R1_R2() {
    local dir="$1" cohort="$2" tier="$3" disease="$4"
    local count=0
    while IFS= read -r -d '' fq1; do
        fq2="${fq1/_R1.fastq.gz/_R2.fastq.gz}"
        sample=$(basename "${fq1}" _R1.fastq.gz)
        add_sample "${sample}" "${cohort}" "${tier}" "${disease}" "${fq1}" "${fq2}"
        ((count++))
    done < <(find "${dir}" -maxdepth 1 -name "*_R1.fastq.gz" -type f -print0 | sort -z)
    echo "  ${cohort}: ${count} samples"
}

# =============================================================================
# TIER 1 — MDS/AML
# =============================================================================
echo "=== TIER 1: MDS/AML ==="

process_R1_R2 "${BASE}/beat"                    "beat"            "1" "MDS_AML"
process_R1_R2 "${BASE}/tcga"                    "tcga"            "1" "MDS_AML"
process_1_2   "${BASE}/Leucegene/SRSF2"         "leucegene_SRSF2" "1" "MDS_AML"
process_1_2   "${BASE}/Leucegene/SF3B1"         "leucegene_SF3B1" "1" "MDS_AML"
process_1_2   "${BASE}/Leucegene/U2AF1"         "leucegene_U2AF1" "1" "MDS_AML"
process_1_2   "${BASE}/Leucegene/ZRSR2"         "leucegene_ZRSR2" "1" "MDS_AML"
process_1_2   "${BASE}/Pellagatti_2018/SF3B1"   "pellagatti_SF3B1" "1" "MDS_AML"
process_1_2   "${BASE}/Pellagatti_2018/SRSF2"   "pellagatti_SRSF2" "1" "MDS_AML"
process_1_2   "${BASE}/Pellagatti_2018/U2AF1"   "pellagatti_U2AF1" "1" "MDS_AML"
process_1_2   "${BASE}/Pellagatti_2018/WT"      "pellagatti_WT"   "1" "MDS_AML"
process_1_2   "${BASE}/Madan_2015/ZRSR2"        "madan_2015_ZRSR2" "1" "MDS_AML"
process_1_2   "${BASE}/Madan_2015/WT"           "madan_2015_WT"   "1" "MDS_AML"

# =============================================================================
# TIER 2 — Normal hematopoietic
# =============================================================================
echo ""
echo "=== TIER 2: Normal hematopoietic ==="

process_1_2 "${BASE}/Maiga_et_al_normal" "maiga_normal" "2" "normal_hematopoietic"
process_1_2 "${BASE}/Madan_et_al_normal" "madan_normal" "2" "normal_hematopoietic"

# =============================================================================
# TIER 3 — Normal other tissue
# =============================================================================
echo ""
echo "=== TIER 3: Normal other tissue ==="

process_1_2 "${BASE}/hpa_2014"  "hpa_2014"  "3" "normal_tissue"
process_1_2 "${BASE}/bodymap2"  "bodymap2"  "3" "normal_tissue"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=== SUMMARY ==="
TOTAL=$(tail -n +2 "${MANIFEST}" | wc -l)
PAIRED=$(awk -F'\t' 'NR>1 && $5=="paired"' "${MANIFEST}" | wc -l)
SINGLE=$(awk -F'\t' 'NR>1 && $5=="single"' "${MANIFEST}" | wc -l)
T1=$(awk -F'\t' 'NR>1 && $3=="1"' "${MANIFEST}" | wc -l)
T2=$(awk -F'\t' 'NR>1 && $3=="2"' "${MANIFEST}" | wc -l)
T3=$(awk -F'\t' 'NR>1 && $3=="3"' "${MANIFEST}" | wc -l)

echo "Total samples : ${TOTAL}"
echo "  Paired-end  : ${PAIRED}"
echo "  Single-end  : ${SINGLE}"
echo "  Tier 1 (MDS/AML)          : ${T1}"
echo "  Tier 2 (Normal hematopoietic): ${T2}"
echo "  Tier 3 (Normal tissue)    : ${T3}"
echo ""

# ── Sanity check: missing R2 files ───────────────────────────────────────────
echo "Checking for missing R2 files..."
MISSING=0
while IFS=$'\t' read -r sample cohort tier disease seq_type fq1 fq2; do
    if [[ "${seq_type}" == "paired" && ( -z "${fq2}" || ! -f "${fq2}" ) ]]; then
        echo "  WARNING: Missing R2 for ${sample}: ${fq2}"
        ((MISSING++))
    fi
done < <(tail -n +2 "${MANIFEST}")

if [[ "${MISSING}" -eq 0 ]]; then
    echo "All R2 files present."
else
    echo "WARNING: ${MISSING} missing R2 files — check naming conventions."
fi

# ── Sanity check: duplicate sample IDs ───────────────────────────────────────
echo ""
echo "Checking for duplicate sample IDs..."
DUPS=$(awk -F'\t' 'NR>1{print $1}' "${MANIFEST}" | sort | uniq -d | wc -l)
if [[ "${DUPS}" -eq 0 ]]; then
    echo "No duplicate sample IDs found."
else
    echo "WARNING: ${DUPS} duplicate sample IDs:"
    awk -F'\t' 'NR>1{print $1}' "${MANIFEST}" | sort | uniq -d
fi

echo ""
echo "Manifest written: ${MANIFEST}"
echo ""
echo "Preview (first 5 samples):"
head -6 "${MANIFEST}" | column -t