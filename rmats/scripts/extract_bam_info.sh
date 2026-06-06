#!/bin/bash
#SBATCH --job-name=bam_qc
#SBATCH --array=0-%NBAMS%        # placeholder — overridden by submit_bam_qc.sh
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
# NOTE: --output and --error are NOT set here — they are passed explicitly
# by submit_bam_qc.sh so logs land in OUT_DIR/logs/.

# =============================================================================
# extract_bam_info.sh
# SLURM array job to extract QC metrics from STAR-aligned BAMs for rMATS prep
#
# Do not submit this script directly. Use the submit_bam_qc.sh wrapper which:
#   1. Builds the BAM manifest from sf_final_metadata_with_del.tsv
#   2. Sets the correct --array range
#   3. Passes --output/--error pointing to OUT_DIR/logs/
#   4. Exports BED_FILE and OUT_DIR via --export
#   5. Submits the aggregation job with --dependency
#
# Variables received via --export from submit_bam_qc.sh:
#   BED_FILE  — GENCODE v49 BED12 for infer_experiment.py
#   OUT_DIR   — output root (bam_qc/)
# =============================================================================

N_READS=200000      # reads to sample for infer_experiment
THREADS=$SLURM_CPUS_PER_TASK
MANIFEST=${OUT_DIR}/bam_manifest.txt

# --- Load tools ---
source /admin/software/anaconda/1.11.1/bin/activate
conda activate gatk

# --- Get this job's BAM from the manifest ---
bam=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")

if [[ -z "$bam" ]]; then
    echo "ERROR: No BAM found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

sample=$(basename "$bam" "_Aligned.sortedByCoord.out.bam")
echo "Array index : $SLURM_ARRAY_TASK_ID"
echo "Sample      : $sample"
echo "BAM         : $bam"
echo "BED_FILE    : $BED_FILE"
echo "OUT_DIR     : $OUT_DIR"
echo ""

# --- Output paths ---
mkdir -p "${OUT_DIR}"/{strand,flagstat,per_sample}
SAMPLE_OUT="${OUT_DIR}/per_sample/${sample}.tsv"

# Skip if already done
if [[ -f "$SAMPLE_OUT" ]]; then
    echo "Output already exists, skipping: $SAMPLE_OUT"
    exit 0
fi

# =============================================================================
# 1. Read length
# =============================================================================
echo "Extracting read lengths..."
read_lengths=$(samtools view -@ "$THREADS" -F 4 -s 0.1 "$bam" 2>/dev/null \
    | head -n 100000 \
    | awk '{print length($10)}' \
    | sort -n)

read_length_max=$(echo "$read_lengths" | tail -1)
read_length_median=$(echo "$read_lengths" | awk '{
    lines[NR] = $1
} END {
    mid = int(NR/2)
    print lines[mid]
}')

read_length_max=${read_length_max:-"NA"}
read_length_median=${read_length_median:-"NA"}

# =============================================================================
# 2. Paired vs single end
# =============================================================================
echo "Checking paired/single end..."
paired_count=$(samtools view -@ "$THREADS" -f 1 -c "$bam" 2>/dev/null || echo 0)
if [[ "$paired_count" -gt 0 ]]; then
    paired_end="yes"
    lib_type_flag="-t paired"
else
    paired_end="no"
    lib_type_flag="-t single"
fi

# =============================================================================
# 3. Flagstat
# =============================================================================
echo "Running flagstat..."
flagstat_file="${OUT_DIR}/flagstat/${sample}_flagstat.txt"
samtools flagstat -@ "$THREADS" "$bam" > "$flagstat_file" 2>/dev/null

total_reads=$(grep "in total"    "$flagstat_file" | awk '{print $1}')
mapped_reads=$(grep "mapped ("   "$flagstat_file" | head -1 | awk '{print $1}')
mapping_rate=$(grep "mapped ("   "$flagstat_file" | head -1 | grep -oP '\(\K[0-9.]+(?=%)')

total_reads=${total_reads:-"NA"}
mapped_reads=${mapped_reads:-"NA"}
mapping_rate=${mapping_rate:-"NA"}

# =============================================================================
# 4. Strandedness via RSeQC
# =============================================================================
# BAMs aligned with --outSAMstrandField intronMotif carry the XS tag derived
# from splice site motif (GT-AG, AT-AC etc.), making infer_experiment.py more
# reliable — especially for reads spanning novel junctions.
echo "Inferring strandedness..."
strand_file="${OUT_DIR}/strand/${sample}_strand.txt"
infer_experiment.py \
    -r "$BED_FILE" \
    -i "$bam" \
    -s "$N_READS" \
    > "$strand_file" 2>/dev/null

forward_pct=$(grep "1++,1--,2+-,2-+" "$strand_file" \
    | grep -oP '(?<=: )[0-9.]+')
reverse_pct=$(grep "1+-,1-+,2++,2--" "$strand_file" \
    | grep -oP '(?<=: )[0-9.]+')
undetermined_pct=$(grep "Fraction of reads failed" "$strand_file" \
    | grep -oP '(?<=: )[0-9.]+')

# Apply NA explicitly for empty results (grep -oP exits 0 even on no match
# in some versions, so || echo "NA" is not a reliable fallback)
[[ -z "$forward_pct" ]]     && forward_pct="NA"
[[ -z "$reverse_pct" ]]     && reverse_pct="NA"
[[ -z "$undetermined_pct" ]] && undetermined_pct="NA"

# Interpret strandedness — guard requires non-empty numeric values
if [[ "$forward_pct" != "NA" && "$reverse_pct" != "NA" ]]; then
    fwd=$(echo "$forward_pct > 0.75" | bc -l 2>/dev/null || echo 0)
    rev=$(echo "$reverse_pct > 0.75" | bc -l 2>/dev/null || echo 0)

    if   [[ "$fwd" -eq 1 ]]; then
        stranded_call="forward"
        rmats_libtype="fr-secondstrand"
    elif [[ "$rev" -eq 1 ]]; then
        stranded_call="reverse"
        rmats_libtype="fr-firststrand"
    else
        stranded_call="unstranded"
        rmats_libtype="fr-unstranded"
    fi
else
    stranded_call="NA"
    rmats_libtype="fr-unstranded"
fi

# =============================================================================
# 5. Write per-sample result
# =============================================================================
echo -e "sample\tbam_path\tread_length_max\tread_length_median\tpaired_end\ttotal_reads\tmapped_reads\tmapping_rate\tforward_pct\treverse_pct\tundetermined_pct\tstranded_call\trmats_libtype\trmats_lib_type_flag" \
    > "$SAMPLE_OUT"

echo -e "${sample}\t${bam}\t${read_length_max}\t${read_length_median}\t${paired_end}\t${total_reads}\t${mapped_reads}\t${mapping_rate}\t${forward_pct}\t${reverse_pct}\t${undetermined_pct}\t${stranded_call}\t${rmats_libtype}\t${lib_type_flag}" \
    >> "$SAMPLE_OUT"

echo ""
echo "Done: $sample"
echo "  Read length (max/median) : ${read_length_max} / ${read_length_median}"
echo "  Paired-end               : ${paired_end}"
echo "  Mapping rate             : ${mapping_rate}%"
echo "  Strandedness             : ${stranded_call} (fwd=${forward_pct}, rev=${reverse_pct})"
echo "  rMATS libType            : ${rmats_libtype}"

conda deactivate