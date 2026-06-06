#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=rmats_post
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=24:00:00
# NOTE: --output and --error are NOT set here — they are passed explicitly
# by submit_rmats_persample.sh so logs land in analysis_dir/logs/.

# =============================================================================
# rmats_post.sh
# 1. Merges per-sample tmp directories into one merged_tmp
# 2. Runs rMATS post using the original comma-separated b1/b2 BAM lists
#    (only including BAMs that were successfully prepped)
#
# Passed via --export:
#   analysis_dir, rmats_container, gtf, star_indices, MANIFEST,
#   b1_post, b2_post
# =============================================================================

echo "Host    : $(hostname)"
echo "Started : $(date)"
echo "GTF     : $gtf"
echo ""

# =============================================================================
# 1. Merge per-sample tmp directories
# =============================================================================

MERGED_TMP=${analysis_dir}/merged_tmp
mkdir -p "$MERGED_TMP"

echo "Merging per-sample tmp directories..."

failed_samples=""
n_copied=0

while IFS='|' read -r sample bam libtype readlen paired condition; do

    sample_tmp="${analysis_dir}/prep/${sample}/tmp"

    if [[ ! -d "$sample_tmp" ]]; then
        echo "  WARNING: tmp dir not found for $sample"
        failed_samples="${failed_samples} ${sample}"
        continue
    fi

    zero_rmat=("${sample_tmp}"/*_0.rmats)
    if [[ ${#zero_rmat[@]} -eq 0 || ! -f "${zero_rmat[0]}" ]]; then
        echo "  WARNING: No _0.rmats file found for $sample — prep likely failed"
        failed_samples="${failed_samples} ${sample}"
        continue
    fi

    dest="${MERGED_TMP}/${sample}.rmats"

    if [[ -f "$dest" ]]; then
        echo "  WARNING: ${sample}.rmats already exists in merged tmp - skipping duplicate"
    else
        cp "${zero_rmat[0]}" "$dest"
        n_copied=$((n_copied + 1))
    fi

done < "$MANIFEST"

N_TOTAL=$(wc -l < "$MANIFEST")
echo ""
echo "Samples in manifest   : $N_TOTAL"
echo ".rmats files merged   : $n_copied"

if [[ -n "$failed_samples" ]]; then
    echo ""
    echo "WARNING: Failed samples excluded from post:"
    echo "$failed_samples" | tr ' ' '\n' | grep -v '^$' | sed 's/^/  /'
    echo ""
    echo "Re-run prep for these samples then resubmit post manually"
fi

# =============================================================================
# 2. If any samples failed, rebuild b1/b2 excluding them
# =============================================================================
if [[ -n "$failed_samples" ]]; then
    echo "Rebuilding b1/b2 lists excluding failed samples..."

    b1_final="${analysis_dir}/bamlists/b1_controls_final.txt"
    b2_final="${analysis_dir}/bamlists/b2_mutants_final.txt"

    > "${b1_final}.tmp"
    > "${b2_final}.tmp"

    while IFS='|' read -r sample bam libtype readlen paired condition; do
        echo "$failed_samples" | grep -qw "$sample" && continue
        if [[ "$condition" == "control" ]]; then
            echo "$bam" >> "${b1_final}.tmp"
        else
            echo "$bam" >> "${b2_final}.tmp"
        fi
    done < "$MANIFEST"

    paste -sd',' "${b1_final}.tmp" > "$b1_final"
    paste -sd',' "${b2_final}.tmp" > "$b2_final"
    rm "${b1_final}.tmp" "${b2_final}.tmp"

    b1_post="$b1_final"
    b2_post="$b2_final"
fi

n_b1=$(tr ',' '\n' < "$b1_post" | grep -c .)
n_b2=$(tr ',' '\n' < "$b2_post" | grep -c .)
echo ""
echo "Controls (b1): $n_b1 samples"
echo "Mutants  (b2): $n_b2 samples"

# =============================================================================
# 3. Derive run parameters from manifest
# =============================================================================

# --readLength is used only to compute IncFormLen/SkipFormLen normalization lengths.
# --variable-read-length handles actual per-read counting regardless of this value.
# Set to 150 (longest read length in cohort) as a safe upper bound.
readlen_post=100
echo "Read length for normalization: ${readlen_post} (--variable-read-length handles actual variation)"
echo "Read length distribution in manifest:"
awk -F'|' '{print $4}' "$MANIFEST" | sort | uniq -c | sort -rn | sed 's/^/  /'

# -t flag: if ANY sample is paired-end, use paired (rMATS post requires
# consistent library type; if you have a true mix, split into separate runs)
has_single=$(awk -F'|' 'NR>1 && $5=="no" {found=1; exit} END {print found+0}' "$MANIFEST")
has_paired=$(awk -F'|' 'NR>1 && $5=="yes" {found=1; exit} END {print found+0}' "$MANIFEST")

if [[ "$has_single" -eq 1 && "$has_paired" -eq 1 ]]; then
    echo "WARNING: manifest contains both paired and single-end samples."
    echo "rMATS post requires a consistent -t flag. Defaulting to -t paired."
    echo "If single-end samples are present intentionally, split into separate runs."
    lib_type_flag="-t paired"
elif [[ "$has_single" -eq 1 ]]; then
    lib_type_flag="-t single"
else
    lib_type_flag="-t paired"
fi
echo "-t flag: $lib_type_flag"
echo ""

# =============================================================================
# 4. Run rMATS post
# =============================================================================
post_out="${analysis_dir}/post/out"
mkdir -p "$post_out"

echo "Running rMATS post..."
echo "Output: $post_out"
echo ""

# --novelSS must match prep — prep used --novelSS so post must too.
# --libType: formally required but ignored by post; fr-unstranded is a safe placeholder.
singularity exec \
    --bind /data1/abdelwao/:/data1/abdelwao/ \
    "$rmats_container" \
    python /rmats/rmats.py \
    --b1 "$b1_post" \
    --b2 "$b2_post" \
    --od "$post_out" \
    --tmp "$MERGED_TMP" \
    --gtf "$gtf" \
    --readLength "$readlen_post" \
    --variable-read-length \
    --allow-clipping \
    --novelSS \
    --libType fr-unstranded \
    $lib_type_flag \
    --nthread "$SLURM_CPUS_PER_TASK" \
    --individual-counts \
    --task post

exit_code=$?

echo ""
echo "Finished : $(date)"
echo "Exit code: $exit_code"

[[ $exit_code -ne 0 ]] && { echo "ERROR: rMATS post failed"; exit $exit_code; }

# =============================================================================
# 5. Results summary
# =============================================================================
# Column 20 in MATS.JC.txt is FDR (q-value) for rMATS v4.x output
echo ""
echo "=============================="
echo "RESULTS SUMMARY"
echo "=============================="
for event in SE A5SS A3SS MXE RI; do
    f="${post_out}/${event}.MATS.JC.txt"
    if [[ -f "$f" ]]; then
        total=$(tail -n +2 "$f" | wc -l)
        sig=$(tail -n +2 "$f" | awk -F'\t' '$20 != "NA" && $20+0 < 0.05' | wc -l)
        echo "  $event : $total total, $sig significant (FDR<0.05)"
    else
        echo "  $event : output not found"
    fi
done

echo ""
echo "Results: $post_out"
echo "Done."