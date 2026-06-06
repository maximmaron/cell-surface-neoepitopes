#!/bin/bash
# =============================================================================
# submit_rmats_persample.sh
#
# Builds rMATS prep manifest directly from the splicing factor metadata TSV
# (sf_final_metadata_with_del.tsv), then submits prep as a SLURM array job
# and post as a dependent job.
#
# b1 = controls : sample_group IN (normal_hematopoietic, normal_tissue)
# b2 = mutants  : sample_group == <mutant_group argument>
#
# dPSI = b2 - b1, so activated cryptic events in mutants have positive dPSI.
#
# Usage:
#   bash submit_rmats_persample.sh \
#       -m /path/to/sf_final_metadata_with_del.tsv \
#       -g SF3B1_hotspot \
#       -q /path/to/bam_summary.tsv \
#       -o /path/to/output_base \
#       [-n analysis_name]   # defaults to <mutant_group>_vs_normal
#
# Metadata TSV columns used:
#   sample_id (1), cohort (2), sample_group (5), bam_path (10)
#
# Requires: bam_summary.tsv from BAM QC pipeline for libtype/readlen/paired lookup
# =============================================================================

set -euo pipefail

# --- Fixed paths ---
rmats_container="/data1/abdelwao/shared/containers/rmats_latest.sif"
star_indices="/data1/abdelwao/maxim/annotations/Homo_sapiens/STAR/gencode_v49_pass2_rl100/"

STRINGTIE_GTF="/data1/abdelwao/maxim/splicing_pipeline/stringtie/merge/merged_final.gtf"
GENCODE_GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"

PREP_SCRIPT="/data1/abdelwao/maxim/splicing_pipeline/rmats/scripts/rmats_prep_persample.sh"
POST_SCRIPT="/data1/abdelwao/maxim/splicing_pipeline/rmats/scripts/rmats_post.sh"

# Control groups — edit here if sample_group labels change
CONTROL_GROUPS=("normal_hematopoietic" "normal_tissue")

# --- Argument parsing ---
metadata_tsv=""
mutant_group=""
qc_manifest=""
out_base=""
analysis_name=""

usage() {
    echo "Usage: $0 -m metadata.tsv -g SF3B1_hotspot -q bam_summary.tsv -o out_dir [-n analysis_name]"
    echo ""
    echo "  -m  sf_final_metadata_with_del.tsv (from mutation calling pipeline)"
    echo "  -g  mutant sample_group value (e.g. SF3B1_hotspot, SRSF2_hotspot)"
    echo "  -q  bam_summary.tsv from BAM QC pipeline"
    echo "  -o  output base directory"
    echo "  -n  analysis name (default: <mutant_group>_vs_normal)"
    echo ""
    echo "Available sample_group values can be checked with:"
    echo "  awk -F'\t' 'NR>1 {print \$5}' <metadata.tsv> | sort -u"
    exit 1
}

while getopts "m:g:q:o:n:h" opt; do
    case $opt in
        m) metadata_tsv="$OPTARG" ;;
        g) mutant_group="$OPTARG" ;;
        q) qc_manifest="$OPTARG" ;;
        o) out_base="$OPTARG" ;;
        n) analysis_name="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -z "$metadata_tsv" || ! -f "$metadata_tsv" ]] && { echo "ERROR: metadata TSV not found: '$metadata_tsv'"; usage; }
[[ -z "$mutant_group" ]]                          && { echo "ERROR: mutant group required (-g)"; usage; }
[[ -z "$qc_manifest"  || ! -f "$qc_manifest"  ]] && { echo "ERROR: QC manifest not found: '$qc_manifest'"; usage; }
[[ -z "$out_base" ]]                              && { echo "ERROR: output directory required (-o)"; usage; }

[[ -z "$analysis_name" ]] && analysis_name="${mutant_group}_vs_normal"

analysis_dir="${out_base}/${analysis_name}"
mkdir -p "${analysis_dir}"/{prep,post,merged_tmp,logs,bamlists}

# --- GTF selection ---
if [[ -f "$STRINGTIE_GTF" ]]; then
    gtf="$STRINGTIE_GTF"
    echo "GTF: StringTie merged -> $gtf"
else
    echo "WARNING: StringTie GTF not found at $STRINGTIE_GTF"
    echo "Falling back to GENCODE v49. Novel assembled junctions will not be in event universe."
    gtf="$GENCODE_GTF"
fi

echo "=============================================="
echo "rMATS per-sample prep setup"
echo "Analysis      : $analysis_name"
echo "Mutant group  : $mutant_group (b2)"
echo "Control groups: ${CONTROL_GROUPS[*]} (b1)"
echo "Metadata TSV  : $metadata_tsv"
echo "QC manifest   : $qc_manifest"
echo "GTF           : $gtf"
echo "Output        : $analysis_dir"
echo "=============================================="
echo ""

# =============================================================================
# Extract BAM paths from metadata TSV
# Columns: sample_id(1) cohort(2) tier(3) disease(4) sample_group(5)
#          mutated_gene(6) protein_change(7) vaf(8) co_mutation(9) bam_path(10)
#
# Build an awk condition string from CONTROL_GROUPS array so this works
# regardless of how many control groups are defined.
# =============================================================================

# Build awk OR-condition from CONTROL_GROUPS array: ($5=="g1" || $5=="g2" || ...)
ctrl_awk_cond='('
for g in "${CONTROL_GROUPS[@]}"; do
    ctrl_awk_cond+='$5=="'"$g"'" || '
done
ctrl_awk_cond="${ctrl_awk_cond% || })"   # strip trailing ' || ' and close paren

# Extract controls
mapfile -t CTRL_BAMS < <(awk -F'\t' \
    "NR>1 && \$10 != \"\" && \$10 != \".\" && ${ctrl_awk_cond} {print \$10}" \
    "$metadata_tsv" | sort -u)

# Extract mutants
mapfile -t MUT_BAMS < <(awk -F'\t' -v grp="$mutant_group" \
    'NR>1 && $5==grp && $10 != "" && $10 != "." {print $10}' \
    "$metadata_tsv" | sort -u)

echo "Controls found in metadata : ${#CTRL_BAMS[@]}"
echo "Mutants found in metadata  : ${#MUT_BAMS[@]}"

if [[ "${#CTRL_BAMS[@]}" -eq 0 ]]; then
    echo "ERROR: No control samples found."
    echo "  Expected sample_group values: ${CONTROL_GROUPS[*]}"
    echo "  Observed sample_group values:"
    awk -F'\t' 'NR>1 {print $5}' "$metadata_tsv" | sort -u | sed 's/^/    /'
    exit 1
fi

if [[ "${#MUT_BAMS[@]}" -eq 0 ]]; then
    echo "ERROR: No mutant samples found for group: $mutant_group"
    echo "  Available sample_group values:"
    awk -F'\t' 'NR>1 {print $5}' "$metadata_tsv" | sort -u | sed 's/^/    /'
    exit 1
fi

# =============================================================================
# QC manifest lookup — get libtype, readlen, paired for each BAM
# bam_summary.tsv columns:
#   1=sample  2=bam_path  3=read_length_max  4=read_length_median
#   5=paired_end  6=total_reads  7=mapped_reads  8=mapping_rate
#   9=forward_pct  10=reverse_pct  11=undetermined_pct
#   12=stranded_call  13=rmats_libtype  14=rmats_lib_type_flag
# =============================================================================

MANIFEST="${analysis_dir}/persample_manifest.txt"
MISSING_LOG="${analysis_dir}/logs/missing_from_qc_manifest.txt"
> "$MANIFEST"
> "$MISSING_LOG"

lookup_and_write() {
    local bam="$1"
    local condition="$2"

    local hit
    hit=$(awk -F'\t' -v b="$bam" 'NR>1 && $2==b {print; exit}' "$qc_manifest")

    if [[ -z "$hit" ]]; then
        hit=$(awk -F'\t' -v b="${bam,,}" 'NR>1 && tolower($2)==b {print; exit}' "$qc_manifest")
    fi

    if [[ -z "$hit" ]]; then
        echo "  WARNING: not in QC manifest: $bam"
        echo "$bam" >> "$MISSING_LOG"
        return 1
    fi

    local sample readlen paired libtype
    sample=$(echo "$hit"  | cut -f1)
    readlen=$(echo "$hit" | cut -f3)
    paired=$(echo "$hit"  | cut -f5)
    libtype=$(echo "$hit" | cut -f13)

    [[ -z "$libtype" || "$libtype" == "NA" ]] && libtype="fr-unstranded"
    [[ -z "$readlen" || "$readlen" == "NA" ]] && readlen="100"
    [[ -z "$paired"  || "$paired"  == "NA" ]] && paired="yes"

    echo "${sample}|${bam}|${libtype}|${readlen}|${paired}|${condition}" >> "$MANIFEST"
}

echo ""
echo "Looking up QC info from bam_summary.tsv..."

n_ctrl_found=0; n_ctrl_miss=0
for bam in "${CTRL_BAMS[@]}"; do
    [[ -z "$bam" ]] && continue
    if lookup_and_write "$bam" "control"; then
        n_ctrl_found=$((n_ctrl_found + 1))
    else
        n_ctrl_miss=$((n_ctrl_miss + 1))
    fi
done
echo "  Controls : ${n_ctrl_found} found, ${n_ctrl_miss} missing"

n_mut_found=0; n_mut_miss=0
for bam in "${MUT_BAMS[@]}"; do
    [[ -z "$bam" ]] && continue
    if lookup_and_write "$bam" "mutant"; then
        n_mut_found=$((n_mut_found + 1))
    else
        n_mut_miss=$((n_mut_miss + 1))
    fi
done
echo "  Mutants  : ${n_mut_found} found, ${n_mut_miss} missing"

N_SAMPLES=$(wc -l < "$MANIFEST")

if [[ "$N_SAMPLES" -eq 0 ]]; then
    echo ""
    echo "ERROR: No samples in manifest."
    echo "  Check that BAM paths in metadata TSV match paths in bam_summary.tsv"
    echo "  Example metadata path : ${CTRL_BAMS[0]:-none}"
    echo "  Example QC manifest   : $(awk -F'\t' 'NR==2{print $2}' "$qc_manifest")"
    exit 1
fi

n_missing_total=$((n_ctrl_miss + n_mut_miss))
if [[ "$n_missing_total" -gt 0 ]]; then
    echo ""
    echo "WARNING: $n_missing_total BAMs not in QC manifest -> $MISSING_LOG"
    echo "  These are EXCLUDED. Re-run BAM QC on missing samples or check paths."
fi

echo ""
echo "Library type breakdown:"
awk -F'|' '{print $3}' "$MANIFEST" | sort | uniq -c | sed 's/^/  /'
echo ""
echo "Read length breakdown:"
awk -F'|' '{print $4}' "$MANIFEST" | sort | uniq -c | sed 's/^/  /'
echo ""
echo "Condition breakdown:"
awk -F'|' '{print $6}' "$MANIFEST" | sort | uniq -c | sed 's/^/  /'
echo ""
echo "Total samples in manifest: $N_SAMPLES"
echo ""

# =============================================================================
# Save comma-separated b1/b2 BAM lists for post step
# =============================================================================
b1_post="${analysis_dir}/bamlists/b1_controls_post.txt"
b2_post="${analysis_dir}/bamlists/b2_mutants_post.txt"

awk -F'|' '$6=="control" {print $2}' "$MANIFEST" | paste -sd',' > "$b1_post"
awk -F'|' '$6=="mutant"  {print $2}' "$MANIFEST" | paste -sd',' > "$b2_post"

n_b1=$(tr ',' '\n' < "$b1_post" | grep -c .)
n_b2=$(tr ',' '\n' < "$b2_post" | grep -c .)
echo "Post BAM lists:"
echo "  b1 controls : $n_b1 -> $b1_post"
echo "  b2 mutants  : $n_b2 -> $b2_post"
echo ""

# =============================================================================
# Submit prep array — pass explicit log paths so logs land in analysis_dir/logs/
# =============================================================================
MAX_INDEX=$(( N_SAMPLES - 1 ))

PREP_JOB_ID=$(sbatch \
    --array=0-${MAX_INDEX} \
    --output="${analysis_dir}/logs/rmats_prep_%A_%a.out" \
    --error="${analysis_dir}/logs/rmats_prep_%A_%a.err" \
    --export=ALL,MANIFEST="${MANIFEST}",analysis_dir="${analysis_dir}",rmats_container="${rmats_container}",gtf="${gtf}",star_indices="${star_indices}" \
    "$PREP_SCRIPT" \
    | awk '{print $NF}')

echo "Prep array submitted : $PREP_JOB_ID ($N_SAMPLES tasks)"

# =============================================================================
# Submit post dependent on all prep tasks completing
# =============================================================================
POST_JOB_ID=$(sbatch \
    --dependency=afterok:${PREP_JOB_ID} \
    --output="${analysis_dir}/logs/rmats_post_%j.out" \
    --error="${analysis_dir}/logs/rmats_post_%j.err" \
    --export=ALL,analysis_dir="${analysis_dir}",rmats_container="${rmats_container}",gtf="${gtf}",star_indices="${star_indices}",MANIFEST="${MANIFEST}",b1_post="${b1_post}",b2_post="${b2_post}" \
    "$POST_SCRIPT" \
    | awk '{print $NF}')

echo "Post job submitted   : $POST_JOB_ID (runs after all prep tasks complete)"
echo ""
echo "Monitor: squeue -u $USER"
echo ""
echo "To resubmit only failed prep tasks:"
echo "  sbatch --array=<failed_indices> \\"
echo "    --output=${analysis_dir}/logs/rmats_prep_%A_%a.out \\"
echo "    --error=${analysis_dir}/logs/rmats_prep_%A_%a.err \\"
echo "    --export=ALL,MANIFEST=${MANIFEST},analysis_dir=${analysis_dir},rmats_container=${rmats_container},gtf=${gtf},star_indices=${star_indices} \\"
echo "    $PREP_SCRIPT"