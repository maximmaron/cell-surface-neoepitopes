#!/bin/bash
# =============================================================================
# 00_prep_bam_lists.sh
# Parse sf_final_metadata_with_del.tsv and write per-group BAM lists for
# LeafCutter differential splicing comparisons.
#
# Outputs one-per-line BAM lists and prints sbatch array commands.
# Run ONCE from the working directory before submitting 01_bam_to_junc.sh.
#
# Usage:
#   bash 00_prep_bam_lists.sh [--metadata FILE] [--outdir DIR] [--comparison NAME]
#
# Comparisons available (set via --comparison):
#   sf3b1_vs_normal_heme   SF3B1 hotspot vs normal hematopoietic (default)
#   sf3b1_vs_normal_all    SF3B1 hotspot vs normal hematopoietic + normal tissue
#   srsf2_vs_normal_heme   SRSF2 hotspot vs normal hematopoietic
#   u2af1_vs_normal_heme   U2AF1 hotspot vs normal hematopoietic
#   zrsr2_vs_normal_heme   ZRSR2 hotspot vs normal hematopoietic
#   all_sf_vs_normal       All SF hotspot (any gene) vs normal hematopoietic
#   custom                 Uses CASE_GROUP and CTRL_GROUP env vars (see below)
#
# For 'custom':
#   CASE_GROUP: comma-separated sample_group values for cases
#   CTRL_GROUP: comma-separated sample_group values for controls
#   Example:
#     CASE_GROUP="SF3B1_hotspot,SRSF2_hotspot" CTRL_GROUP="normal_hematopoietic" \
#       bash 00_prep_bam_lists.sh --comparison custom
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
METADATA="/data1/abdelwao/maxim/splicing_pipeline/gatk/sf_final_metadata_with_del.tsv"
#WORK_DIR="/data1/abdelwao/shared/splicing_analysis/leafcutter/ZRSR2_vs_normal"
LIST_DIR="${WORK_DIR}/bam_lists"
#COMPARISON="sf3b1_vs_normal_heme"

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --metadata) METADATA="$2"; shift 2 ;;
        --outdir)   LIST_DIR="$2";  shift 2 ;;
        --comparison) COMPARISON="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

mkdir -p "${LIST_DIR}"

# ── Column indices (0-based for awk, extracted from header) ───────────────────
# Rather than hardcoding positions, resolve dynamically from the header row.
get_col() {
    # get_col <header_line> <field_name>  → 1-based awk column index
    local header="$1" field="$2"
    echo "${header}" | tr '\t' '\n' | grep -n "^${field}$" | cut -d: -f1
}

HEADER=$(head -1 "${METADATA}")

COL_GROUP=$(get_col "${HEADER}" "sample_group")
COL_BAM=$(get_col "${HEADER}" "bam_path")
COL_SF3B1_CALL=$(get_col "${HEADER}" "SF3B1_call")
COL_SRSF2_CALL=$(get_col "${HEADER}" "SRSF2_call")
COL_U2AF1_CALL=$(get_col "${HEADER}" "U2AF1_call")
COL_ZRSR2_CALL=$(get_col "${HEADER}" "ZRSR2_call")
COL_DISEASE=$(get_col "${HEADER}" "disease")

echo "Metadata:   ${METADATA}"
echo "Output dir: ${LIST_DIR}"
echo "Comparison: ${COMPARISON}"
echo ""

# ── Define case/control sample_group values per comparison ───────────────────
case "${COMPARISON}" in
    sf3b1_vs_normal_heme)
        CASE_LABEL="sf3b1"
        CTRL_LABEL="normal"
        CASE_GROUPS="SF3B1_hotspot"
        CTRL_GROUPS="normal_hematopoietic"
        ;;
    sf3b1_vs_normal_all)
        CASE_LABEL="sf3b1"
        CTRL_LABEL="normal"
        CASE_GROUPS="SF3B1_hotspot"
        CTRL_GROUPS="normal_hematopoietic,normal_tissue"
        ;;
    srsf2_vs_normal_heme)
        CASE_LABEL="srsf2"
        CTRL_LABEL="normal"
        CASE_GROUPS="SRSF2_hotspot"
        CTRL_GROUPS="normal_hematopoietic"
        ;;
    u2af1_vs_normal_heme)
        CASE_LABEL="u2af1"
        CTRL_LABEL="normal"
        CASE_GROUPS="U2AF1_hotspot"
        CTRL_GROUPS="normal_hematopoietic"
        ;;
    zrsr2_vs_normal_heme)
        CASE_LABEL="zrsr2"
        CTRL_LABEL="normal"
        CASE_GROUPS="ZRSR2_hotspot"
        CTRL_GROUPS="normal_hematopoietic"
        ;;
    all_sf_vs_normal)
        CASE_LABEL="any_sf_mut"
        CTRL_LABEL="normal"
        CASE_GROUPS="SF3B1_hotspot,SRSF2_hotspot,U2AF1_hotspot,ZRSR2_hotspot"
        CTRL_GROUPS="normal_hematopoietic"
        ;;
    custom)
        CASE_LABEL="cases"
        CTRL_LABEL="controls"
        if [[ -z "${CASE_GROUP:-}" ]]; then
            echo "ERROR: --comparison custom requires CASE_GROUP env var to be set"
            echo "  Example: CASE_GROUP=\"SF3B1_hotspot\" CTRL_GROUP=\"normal_hematopoietic\" bash 00_prep_bam_lists.sh --comparison custom"
            exit 1
        fi
        if [[ -z "${CTRL_GROUP:-}" ]]; then
            echo "ERROR: --comparison custom requires CTRL_GROUP env var to be set"
            echo "  Example: CASE_GROUP=\"SF3B1_hotspot\" CTRL_GROUP=\"normal_hematopoietic\" bash 00_prep_bam_lists.sh --comparison custom"
            exit 1
        fi
        CASE_GROUPS="${CASE_GROUP}"
        CTRL_GROUPS="${CTRL_GROUP}"
        ;;
    *)
        echo "ERROR: Unknown comparison '${COMPARISON}'"
        echo "Valid options: sf3b1_vs_normal_heme, sf3b1_vs_normal_all, srsf2_vs_normal_heme,"
        echo "               u2af1_vs_normal_heme, zrsr2_vs_normal_heme, all_sf_vs_normal, custom"
        exit 1
        ;;
esac

# ── Helper: extract BAM paths matching a pipe-separated group list ────────────
# Converts comma-separated group string to awk regex: "A,B,C" → "^(A|B|C)$"
groups_to_regex() {
    echo "$1" | tr ',' '|' | sed 's/|/\\|/g' | awk '{print "^(" $0 ")$"}'
}

extract_bams() {
    local groups="$1"
    local outfile="$2"
    local regex
    regex=$(groups_to_regex "${groups}")

    awk -F'\t' -v col_group="${COL_GROUP}" \
               -v col_bam="${COL_BAM}" \
               -v regex="${regex}" \
    'NR > 1 {
        grp = $col_group
        bam = $col_bam
        if (grp ~ regex && bam != "" && bam != ".") {
            print bam
        }
    }' "${METADATA}" | sort -u > "${outfile}"
}

# ── Extract BAM lists ─────────────────────────────────────────────────────────
CASE_FILE="${LIST_DIR}/${CASE_LABEL}_bams.txt"
CTRL_FILE="${LIST_DIR}/${CTRL_LABEL}_bams.txt"

extract_bams "${CASE_GROUPS}" "${CASE_FILE}"
extract_bams "${CTRL_GROUPS}" "${CTRL_FILE}"

N_CASE=$(wc -l < "${CASE_FILE}")
N_CTRL=$(wc -l < "${CTRL_FILE}")

# ── Validate BAM paths exist ──────────────────────────────────────────────────
echo "Validating BAM paths..."
MISSING=0
while IFS= read -r bam; do
    if [[ ! -f "${bam}" ]]; then
        echo "  MISSING: ${bam}"
        MISSING=$((MISSING + 1))
    fi
done < <(cat "${CASE_FILE}" "${CTRL_FILE}")

if [[ ${MISSING} -gt 0 ]]; then
    echo "WARNING: ${MISSING} BAM file(s) not found on disk. Check paths above."
else
    echo "All BAM paths verified."
fi

echo ""
echo "Case  (${CASE_LABEL}): ${N_CASE} samples  → ${CASE_FILE}"
echo "  Groups: ${CASE_GROUPS}"
echo ""
echo "Control (${CTRL_LABEL}): ${N_CTRL} samples → ${CTRL_FILE}"
echo "  Groups: ${CTRL_GROUPS}"
echo ""

# ── Print sample group summary from metadata ─────────────────────────────────
echo "Sample group counts in metadata:"
awk -F'\t' -v col="${COL_GROUP}" 'NR > 1 {print $col}' "${METADATA}" \
    | sort | uniq -c | sort -rn | awk '{printf "  %6d  %s\n", $1, $2}'
echo ""

# ── Print sbatch commands ─────────────────────────────────────────────────────
echo "Submit jobs with:"
echo "  sbatch --export=COHORT=${CASE_LABEL}  --array=1-${N_CASE}  01_bam_to_junc.sh"
echo "  sbatch --export=COHORT=${CTRL_LABEL} --array=1-${N_CTRL} 01_bam_to_junc.sh"