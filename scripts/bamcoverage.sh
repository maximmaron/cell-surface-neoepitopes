#!/bin/bash
# ---------------------------------------------------------------------------
# bamcoverage.sh — SLURM array job for bamCoverage.
# Reads bam_path directly from sf_final_metadata_with_del.tsv.
# BigWigs written alongside each BAM.
#
# Usage:
#   bash bamcoverage.sh \
#       --metadata   sf_final_metadata_with_del.tsv \
#       [--bin-size  1] \
#       [--normalize CPM] \
#       [--cpus      16] \
#       [--mem       32G] \
#       [--time      2:00:00] \
#       [--partition cpu]
# ---------------------------------------------------------------------------
set -euo pipefail

METADATA=""
BIN_SIZE=1
NORMALIZE="CPM"
CPUS=16
MEM="32G"
TIME="2:00:00"
PARTITION="cpu"

while [[ $# -gt 0 ]]; do
    case $1 in
        --metadata)   METADATA="$2";  shift 2 ;;
        --bin-size)   BIN_SIZE="$2";  shift 2 ;;
        --normalize)  NORMALIZE="$2"; shift 2 ;;
        --cpus)       CPUS="$2";      shift 2 ;;
        --mem)        MEM="$2";       shift 2 ;;
        --time)       TIME="$2";      shift 2 ;;
        --partition)  PARTITION="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

[[ -z "$METADATA" ]] && { echo "ERROR: --metadata required"; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: not found: $METADATA"; exit 1; }

METADATA=$(realpath "$METADATA")
WORKDIR="$(dirname "$METADATA")/bamcoverage"
LOGDIR="${WORKDIR}/slurm_logs"
mkdir -p "$LOGDIR"

# Pull sample_id (col1) and bam_path (col9) — skip header and blank/dot paths
BAM_LIST="${WORKDIR}/bam_list.tsv"
tail -n +2 "$METADATA" | awk -F'\t' '$9 != "." && $9 != "" {print $1"\t"$9}' > "$BAM_LIST"

N_BAMS=$(wc -l < "$BAM_LIST")
echo "Found $N_BAMS BAMs"
[[ $N_BAMS -eq 0 ]] && { echo "ERROR: No BAMs found"; exit 1; }
LAST_IDX=$(( N_BAMS - 1 ))

SLURM_SCRIPT="${WORKDIR}/slurm_bamcoverage.sh"
cat > "$SLURM_SCRIPT" << SLURM
#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --array=0-${LAST_IDX}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --partition=${PARTITION}
#SBATCH --output=${LOGDIR}/bamcoverage_%A_%a.out
#SBATCH --error=${LOGDIR}/bamcoverage_%A_%a.err

source /admin/software/anaconda/1.11.1/bin/activate
conda activate deeptools

LINE=\$(sed -n "\$((SLURM_ARRAY_TASK_ID + 1))p" "${BAM_LIST}")
SAMPLE_ID=\$(echo "\$LINE" | cut -f1)
BAM=\$(echo "\$LINE"       | cut -f2)

OUT_BW="\$(dirname \$BAM)/\${SAMPLE_ID}.bw"

echo "Task \$SLURM_ARRAY_TASK_ID | \$SAMPLE_ID"
echo "  BAM: \$BAM"
echo "  BW:  \$OUT_BW"

bamCoverage \\
    --bam                "\$BAM" \\
    --outFileName        "\$OUT_BW" \\
    --outFileFormat      bigwig \\
    --binSize            ${BIN_SIZE} \\
    --normalizeUsing     ${NORMALIZE} \\
    --minMappingQuality  20 \\
    --numberOfProcessors ${CPUS}
SLURM

chmod +x "$SLURM_SCRIPT"

echo "Submitting array job for $N_BAMS BAMs..."
JOB_ID=$(sbatch --parsable "$SLURM_SCRIPT")
echo "Submitted job: $JOB_ID"
echo "BigWigs will be written alongside each BAM"