MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
missing_tasks=()
task_num=1

while IFS=$'\t' read -r sample_id cohort tier disease seq_type fq1 fq2; do
    out_dir=$(dirname "${fq1}")
    sj="${out_dir}/${sample_id}_SJ.out.tab"
    log="${out_dir}/${sample_id}_Log.final.out"
    if [[ ! -f "${sj}" || ! -f "${log}" ]]; then
        missing_tasks+=("${task_num}")
    fi
    ((task_num++))
done < <(tail -n +2 "${MANIFEST}")

# Print count
echo "Missing: ${#missing_tasks[@]}"

# Format as comma-separated for sbatch --array
ARRAY_SPEC=$(IFS=','; echo "${missing_tasks[*]}")
echo ""
echo "Submit with:"
echo "  sbatch --array=${ARRAY_SPEC}%10 01_star_pass1.sh"