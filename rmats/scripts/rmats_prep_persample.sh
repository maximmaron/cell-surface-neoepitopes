#!/bin/bash
#SBATCH --job-name=rmats_prep
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
# NOTE: --output and --error are NOT set here — they are passed explicitly
# by submit_rmats_persample.sh so logs land in analysis_dir/logs/.
# Setting them here would cause them to resolve relative to wherever sbatch
# is called from, which is wrong.

# =============================================================================
# rmats_prep_persample.sh
# One array task per sample, each prepped with its own correct libType
# and readLength looked up from the QC manifest.
#
# Manifest format: sample|bam|libtype|readlen|paired|condition
#
# Passed via --export:
#   MANIFEST, analysis_dir, rmats_container, gtf, star_indices
# =============================================================================

echo "Array task : $SLURM_ARRAY_TASK_ID"
echo "Host       : $(hostname)"
echo "Started    : $(date)"
echo ""

# --- Get this sample ---
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")

if [[ -z "$LINE" ]]; then
    echo "ERROR: No entry for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

sample=$(echo "$LINE"    | cut -d'|' -f1)
bam=$(echo "$LINE"       | cut -d'|' -f2)
libtype=$(echo "$LINE"   | cut -d'|' -f3)
readlen=$(echo "$LINE"   | cut -d'|' -f4)
paired=$(echo "$LINE"    | cut -d'|' -f5)
condition=$(echo "$LINE" | cut -d'|' -f6)

[[ "$paired" == "yes" ]] && lib_type_flag="-t paired" || lib_type_flag="-t single"

echo "Sample    : $sample"
echo "BAM       : $bam"
echo "GTF       : $gtf"
echo "libType   : $libtype"
echo "readLength: $readlen"
echo "Paired    : $paired ($lib_type_flag)"
echo "Condition : $condition"
echo ""

# --- Validate ---
if [[ ! -f "$bam" ]]; then
    echo "ERROR: BAM not found: $bam"
    exit 1
fi

if [[ ! -f "$gtf" ]]; then
    echo "ERROR: GTF not found: $gtf"
    echo "Check that StringTie merge has completed or GENCODE GTF path is correct"
    exit 1
fi

# --- Skip if already completed ---
sample_dir=${analysis_dir}/prep/${sample}
if [[ -d "${sample_dir}/tmp" ]] && \
   [[ $(ls "${sample_dir}/tmp/"*.rmats 2>/dev/null | wc -l) -gt 0 ]]; then
    echo "Already complete, skipping: $sample"
    exit 0
fi

mkdir -p "${sample_dir}"/{out,tmp,bamlists}

# For per-sample prep, same BAM goes into both b1 and b2.
# rMATS requires both even for single-sample prep step.
echo "$bam" > "${sample_dir}/bamlists/b1.txt"
echo "$bam" > "${sample_dir}/bamlists/b2.txt"

# --- Run rMATS prep ---
# --novelSS: uses XS tag (set by --outSAMstrandField intronMotif in STAR) to
#   assign strand to novel splice sites. Correct and required combination.
# --allow-clipping: handles soft-clipped reads at novel junctions in
#   pooled two-pass BAMs.
# --libType and -t: set per-sample from QC manifest.
echo "Running rMATS prep..."

singularity exec \
    --bind /data1/abdelwao/:/data1/abdelwao/ \
    "$rmats_container" \
    python /rmats/rmats.py \
    --b1 "${sample_dir}/bamlists/b1.txt" \
    --b2 "${sample_dir}/bamlists/b2.txt" \
    --od "${sample_dir}/out" \
    --tmp "${sample_dir}/tmp" \
    --gtf "$gtf" \
    --readLength "$readlen" \
    --variable-read-length \
    --allow-clipping \
    --novelSS \
    --libType "$libtype" \
    $lib_type_flag \
    --nthread "$SLURM_CPUS_PER_TASK" \
    --individual-counts \
    --task prep

exit_code=$?

echo ""
echo "Finished : $(date)"
echo "Exit code: $exit_code"

if [[ $exit_code -ne 0 ]]; then
    echo "ERROR: rMATS prep failed for $sample"
    exit $exit_code
fi

n_rmat=$(ls "${sample_dir}/tmp/"*.rmats 2>/dev/null | wc -l)
echo "rmats files: $n_rmat"
[[ "$n_rmat" -eq 0 ]] && { echo "ERROR: No .rmats files generated"; exit 1; }

echo "Prep complete: $sample"