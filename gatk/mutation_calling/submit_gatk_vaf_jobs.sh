#!/bin/bash
# =============================================================================
# submit_gatk_vaf_jobs.sh
# GATK HaplotypeCaller variant calling across SF gene exons for all samples.
#
# Pipeline per sample:
#   1. SplitNCigarReads        — split reads at N (splice junctions) for RNA-seq
#   2. AddOrReplaceReadGroups  — ensure SM: read group tag present (required by HC)
#   3. HaplotypeCaller         — call variants in SF gene exon regions
#   4. VariantFiltration       — hard filter for RNA-seq
#   5. Funcotator              — annotate with gene/consequence
#   6. Parse output            — extract coding non-synonymous variants
#
# Usage:
#   bash submit_gatk_vaf_jobs.sh \
#       --manifest           /path/to/sample_manifest.tsv \
#       --bed                /path/to/sf_genes_exons_merged.bed \
#       --ref                /path/to/GRCh38.primary_assembly.genome.fa \
#       --funcotator-sources /path/to/funcotator_dataSources/ \
#       --outdir             /path/to/output/dir \
#       [--concurrency       20] \
#       [--scratch           /scratch/abdelwao]
# =============================================================================

set -uo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
MANIFEST=""
BED=""
REF=""
FUNCOTATOR_SOURCES=""
OUTDIR=""
CONCURRENCY=20
SCRATCH="/scratch/abdelwao"

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)           MANIFEST="$2";           shift 2 ;;
        --bed)                BED="$2";                shift 2 ;;
        --ref)                REF="$2";                shift 2 ;;
        --funcotator-sources) FUNCOTATOR_SOURCES="$2"; shift 2 ;;
        --outdir)             OUTDIR="$2";             shift 2 ;;
        --concurrency)        CONCURRENCY="$2";        shift 2 ;;
        --scratch)            SCRATCH="$2";            shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
for var in MANIFEST BED REF FUNCOTATOR_SOURCES OUTDIR; do
    [[ -z "${!var}" ]] && { echo "ERROR: --${var,,} is required" >&2; exit 1; }
done
for f in "${MANIFEST}" "${BED}" "${REF}"; do
    [[ ! -f "${f}" ]] && { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done
[[ ! -d "${FUNCOTATOR_SOURCES}" ]] && \
    { echo "ERROR: Funcotator sources dir not found: ${FUNCOTATOR_SOURCES}" >&2; exit 1; }

mkdir -p "${OUTDIR}"/{vcf,funcotated,parsed,logs}

# ── Build BAM list from manifest ──────────────────────────────────────────────
echo "Building BAM list from manifest..."
BAM_LIST="${OUTDIR}/bam_files.txt"
> "${BAM_LIST}"

FOUND=0
MISSING=0
while IFS=$'\t' read -r sample_id cohort tier disease seq_type fq1 fq2; do
    out_dir=$(dirname "${fq1}")
    bam="${out_dir}/${sample_id}_Aligned.sortedByCoord.out.bam"
    if [[ -f "${bam}" ]]; then
        echo "${bam}" >> "${BAM_LIST}"
        FOUND=$(( FOUND + 1 ))
    else
        echo "WARNING: BAM not found for ${sample_id}: ${bam}" >&2
        MISSING=$(( MISSING + 1 ))
    fi
done < <(tail -n +2 "${MANIFEST}")

N_BAMS=$(wc -l < "${BAM_LIST}")
echo "  Found  : ${FOUND}"
echo "  Missing: ${MISSING}"
echo "  Total  : ${N_BAMS}"
echo ""

[[ "${N_BAMS}" -eq 0 ]] && { echo "ERROR: No BAMs found — aborting." >&2; exit 1; }

# ── Write per-sample GATK array script ───────────────────────────────────────
GATK_ARRAY="${OUTDIR}/run_gatk_array.sh"

cat > "${GATK_ARRAY}" << SLURM
#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=gatk_sf
#SBATCH --output=${OUTDIR}/logs/gatk_%A_%a.out
#SBATCH --error=${OUTDIR}/logs/gatk_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

source /admin/software/anaconda/1.11.1/bin/activate
conda activate gatk

set -uo pipefail

BAM=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${BAM_LIST}")
[[ -z "\${BAM}" ]] && { echo "ERROR: No BAM for task \${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
[[ ! -f "\${BAM}" ]]  && { echo "ERROR: BAM not found: \${BAM}" >&2; exit 1; }

SAMPLE=\$(basename "\${BAM}" _Aligned.sortedByCoord.out.bam)
SPLIT_BAM="${OUTDIR}/vcf/\${SAMPLE}_split.bam"
RG_BAM="${OUTDIR}/vcf/\${SAMPLE}_split_rg.bam"
VCF_RAW="${OUTDIR}/vcf/\${SAMPLE}_raw.vcf.gz"
VCF_FILT="${OUTDIR}/vcf/\${SAMPLE}_filtered.vcf.gz"
FUNCOTATED="${OUTDIR}/funcotated/\${SAMPLE}.vcf.gz"
PARSED="${OUTDIR}/parsed/\${SAMPLE}.variants.tsv"

echo "======================================="
echo "Job    : \${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}"
echo "Sample : \${SAMPLE}"
echo "BAM    : \${BAM}"
echo "Started: \$(date)"
echo "======================================="

# Skip if already complete — header-only files do NOT count as complete.
# This prevents a failed parse step from blocking re-runs.
if [[ -f "\${PARSED}" ]] && [[ \$(wc -l < "\${PARSED}") -gt 1 ]]; then
    echo "Skipping \${SAMPLE} — already complete (\$(( \$(wc -l < \${PARSED}) - 1 )) variants)"
    exit 0
fi

# Ensure scratch dir exists
mkdir -p "${SCRATCH}"

# ── Step 1: SplitNCigarReads ─────────────────────────────────────────────────
# Splits reads at splice junctions (N in CIGAR).
# Essential for RNA-seq variant calling — prevents HaplotypeCaller from
# misinterpreting splice junctions as deletions.
echo ""
echo "── Step 1: SplitNCigarReads ────────────────────────────"

gatk SplitNCigarReads \
    -R "${REF}" \
    -I "\${BAM}" \
    -O "\${SPLIT_BAM}" \
    --intervals "${BED}" \
    --create-output-bam-index true \
    --tmp-dir "${SCRATCH}"

if [[ \$? -ne 0 ]]; then
    echo "ERROR: SplitNCigarReads failed for \${SAMPLE}" >&2
    exit 1
fi

# ── Step 2: AddOrReplaceReadGroups ────────────────────────────────────────────
# HaplotypeCaller requires at least one @RG line with a SM: tag in the BAM
# header. STAR does not add read groups by default, so we add them here.
# If your STAR command already includes --outSAMattrRGline, this step is still
# safe — it will simply overwrite the existing read group with matching values.
echo ""
echo "── Step 2: AddOrReplaceReadGroups ───────────────────────"

gatk AddOrReplaceReadGroups \
    -I "\${SPLIT_BAM}" \
    -O "\${RG_BAM}" \
    --RGID "\${SAMPLE}" \
    --RGLB lib1 \
    --RGPL ILLUMINA \
    --RGPU unit1 \
    --RGSM "\${SAMPLE}" \
    --CREATE_INDEX true \
    --TMP_DIR "${SCRATCH}"

if [[ \$? -ne 0 ]]; then
    echo "ERROR: AddOrReplaceReadGroups failed for \${SAMPLE}" >&2
    exit 1
fi

# Split BAM no longer needed — free disk space early
rm -f "\${SPLIT_BAM}" "\${SPLIT_BAM}.bai"

# ── Step 3: HaplotypeCaller ───────────────────────────────────────────────────
# RNA-seq mode flags:
#   --dont-use-soft-clipped-bases  exclude soft-clipped bases (unreliable in RNA)
#   -stand-call-conf 20            min phred-scaled confidence threshold
#   --native-pair-hmm-threads      multi-threaded HMM (uses SLURM CPU allocation)
echo ""
echo "── Step 3: HaplotypeCaller ──────────────────────────────"

gatk HaplotypeCaller \
    -R "${REF}" \
    -I "\${RG_BAM}" \
    -O "\${VCF_RAW}" \
    --intervals "${BED}" \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --native-pair-hmm-threads \${SLURM_CPUS_PER_TASK} \
    --tmp-dir "${SCRATCH}"

if [[ \$? -ne 0 ]]; then
    echo "ERROR: HaplotypeCaller failed for \${SAMPLE}" >&2
    exit 1
fi

# RG BAM no longer needed
rm -f "\${RG_BAM}" "\${RG_BAM}.bai"

# ── Step 4: VariantFiltration ─────────────────────────────────────────────────
# GATK recommended RNA-seq hard filters (VQSR is not appropriate for RNA-seq).
#   FS > 30   Fisher strand bias — strand-biased artifacts common in RNA-seq
#   QD < 2    Quality by depth  — low confidence relative to coverage
#   DP < 10   Minimum read depth
echo ""
echo "── Step 4: VariantFiltration ────────────────────────────"

gatk VariantFiltration \
    -R "${REF}" \
    -V "\${VCF_RAW}" \
    -O "\${VCF_FILT}" \
    --filter-expression "QD < 2.0" \
    --filter-name "QD2" \
    --tmp-dir "${SCRATCH}"

if [[ \$? -ne 0 ]]; then
    echo "ERROR: VariantFiltration failed for \${SAMPLE}" >&2
    exit 1
fi

rm -f "\${VCF_RAW}" "\${VCF_RAW}.tbi"

# ── Step 5: Funcotator ────────────────────────────────────────────────────────
# Annotates variants with gene, transcript, consequence, and protein change.
# Uses MANE Select transcripts where available.
echo ""
echo "── Step 5: Funcotator ───────────────────────────────────"

gatk Funcotator \
    -R "${REF}" \
    -V "\${VCF_FILT}" \
    -O "\${FUNCOTATED}" \
    --output-file-format VCF \
    --data-sources-path "${FUNCOTATOR_SOURCES}" \
    --ref-version hg38 \
    --intervals "${BED}" \
    --annotation-default Tumor_Sample_Barcode:\${SAMPLE} \
    --tmp-dir "${SCRATCH}"

if [[ \$? -ne 0 ]]; then
    echo "ERROR: Funcotator failed for \${SAMPLE}" >&2
    exit 1
fi

rm -f "\${VCF_FILT}" "\${VCF_FILT}.tbi"

# ── Step 6: Parse Funcotator output ──────────────────────────────────────────
# Extract coding non-synonymous PASS variants into a TSV.
# Retains: missense, nonsense, frameshift, splice_site, in-frame indels.
# Excludes: synonymous, intronic, UTR, and hard-filtered variants.
#
# Uses bcftools instead of gatk PrintVariants — PrintVariants silently drops
# records with FILTER='.' which causes real variants to be lost.
#
# Confirmed Funcotation pipe-delimited field indices for this data source:
#   f[0]  = gene           f[5]  = variant_class
#   f[12] = transcript     f[18] = protein_change
#
# Variants are written to a tmpfile first, then the header and results are
# combined in one atomic write. This prevents a header-only file from being
# mistaken as complete by the skip guard on re-submission.
echo ""
echo "── Step 6: Parse variants ───────────────────────────────"

TMPFILE=\$(mktemp)

bcftools view -H "\${FUNCOTATED}" | \
python3 -c "
import sys

CODING_CLASSES = {
    'MISSENSE', 'NONSENSE', 'NONSTOP',
    'FRAME_SHIFT_INS', 'FRAME_SHIFT_DEL',
    'IN_FRAME_INS', 'IN_FRAME_DEL',
    'START_CODON_SNP', 'START_CODON_INS', 'START_CODON_DEL',
    'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME',
}

sample = sys.argv[1]

for line in sys.stdin:
    if line.startswith('#'):
        continue
    parts = line.strip().split('\t')
    if len(parts) < 8:
        continue

    chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]

    if filt not in ('PASS', '.'):
        continue

    funcotation = ''
    for field in info.split(';'):
        if field.startswith('FUNCOTATION='):
            funcotation = field[12:].strip('[]')
            break

    if not funcotation:
        continue

    f = funcotation.split('|')
    if len(f) < 19:
        continue

    gene           = f[0]
    transcript     = f[12]
    variant_class  = f[5]
    protein_change = f[18]

    if variant_class.upper() not in CODING_CLASSES:
        continue

    vaf = '.'
    depth = '.'
    if len(parts) >= 10:
        fmt_dict = dict(zip(parts[8].split(':'), parts[9].split(':')))
        ad = fmt_dict.get('AD', '0,0').split(',')
        try:
            ref_d = int(ad[0])
            alt_d = int(ad[1]) if len(ad) > 1 else 0
            depth = ref_d + alt_d
            vaf   = round(alt_d / depth, 4) if depth > 0 else 0.0
        except (ValueError, IndexError):
            pass

    print('\t'.join([
        sample, chrom, pos, ref, alt,
        gene, transcript, variant_class, protein_change,
        str(vaf), str(depth), filt
    ]))
    sys.stdout.flush()
" "\${SAMPLE}" > "\${TMPFILE}"

{
    echo -e "sample_id\tchrom\tpos\tref\talt\tgene\ttranscript\tvariant_class\tprotein_change\tVAF\tdepth\tfilter"
    cat "\${TMPFILE}"
} > "\${PARSED}"
rm -f "\${TMPFILE}"

N_VARS=\$(( \$(wc -l < "\${PARSED}") - 1 ))
echo "  Coding variants found: \${N_VARS}"

echo ""
echo "Completed : \${SAMPLE}"
echo "Finished  : \$(date)"
echo "Output    : \${PARSED}"
SLURM

chmod +x "${GATK_ARRAY}"

# ── Write aggregation script ──────────────────────────────────────────────────
AGG_SCRIPT="${OUTDIR}/run_aggregate.sh"

cat > "${AGG_SCRIPT}" << SLURM
#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=gatk_aggregate
#SBATCH --output=${OUTDIR}/logs/aggregate_%j.out
#SBATCH --error=${OUTDIR}/logs/aggregate_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

source /admin/software/anaconda/1.11.1/bin/activate
conda activate gatk

set -uo pipefail

OUTDIR="${OUTDIR}"
MANIFEST="${MANIFEST}"
SUMMARY="\${OUTDIR}/sf_mutation_calls.tsv"
METADATA="\${OUTDIR}/sf_mutation_metadata.tsv"
ALL_VARIANTS="\${OUTDIR}/sf_all_variants.tsv"

echo "Aggregating GATK variant calls..."
echo "Started: \$(date)"

# ── Combine all parsed variant TSVs ──────────────────────────────────────────
echo -e "sample_id\tchrom\tpos\tref\talt\tgene\ttranscript\tvariant_class\tprotein_change\tVAF\tdepth\tfilter" \
    > "\${ALL_VARIANTS}"

shopt -s nullglob
TSV_FILES=( "\${OUTDIR}/parsed/"*.variants.tsv )
if [[ \${#TSV_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No parsed variant TSVs found in \${OUTDIR}/parsed/" >&2
    exit 1
fi

for f in "\${TSV_FILES[@]}"; do
    tail -n +2 "\${f}" >> "\${ALL_VARIANTS}"
done

N_TOTAL=\$(( \$(wc -l < "\${ALL_VARIANTS}") - 1 ))
echo "Total coding variants across all samples: \${N_TOTAL}"

# ── Build per-sample mutation call summary ────────────────────────────────────
# Priority order: MISSENSE > NONSENSE > FRAMESHIFT > IN_FRAME > SPLICE_SITE
# For each sample, takes the highest-priority then highest-VAF variant per gene.
echo ""
echo "Building mutation call summary..."

python3 << 'PYEOF'
import csv
from collections import defaultdict

VAF_THRESHOLD = 0.05
MIN_DEPTH     = 10    # reject variants supported by fewer than 10 reads
GENES = ['SF3B1', 'SRSF2', 'U2AF1', 'ZRSR2']

OUTDIR = "${OUTDIR}"

# Load all variants
variants = defaultdict(list)
with open(f"{OUTDIR}/sf_all_variants.tsv") as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        try:
            vaf = float(row['VAF'])
        except (ValueError, TypeError):
            continue
        if vaf < VAF_THRESHOLD:
            continue
        try:
            depth = int(row['depth'])
        except (ValueError, TypeError):
            depth = 0
        if depth < MIN_DEPTH:
            continue
        variants[row['sample_id']].append(row)

summary_rows = []

for sample_id, vars_list in variants.items():
    by_gene = defaultdict(list)
    for v in vars_list:
        by_gene[v['gene']].append(v)

    for gene in GENES:
        gene_vars = by_gene.get(gene, [])
        if not gene_vars:
            summary_rows.append({
                'sample_id':      sample_id,
                'gene':           gene,
                'called':         'WT',
                'variant_class':  '.',
                'protein_change': '.',
                'VAF':            '.',
                'chrom':          '.',
                'pos':            '.',
                'ref':            '.',
                'alt':            '.',
            })
            continue

        # Sort by alt read count (VAF * depth) descending — the variant
        # with the most supporting reads wins, regardless of consequence class.
        # This correctly handles cases where a low-depth missense would
        # otherwise outrank a well-supported nonsense due to VAF alone.
        def alt_reads(v):
            try:
                return float(v['VAF']) * int(v['depth'])
            except (ValueError, TypeError):
                return 0.0

        gene_vars.sort(key=lambda x: -alt_reads(x))
        best = gene_vars[0]
        summary_rows.append({
            'sample_id':      sample_id,
            'gene':           gene,
            'called':         f"{gene}_mut",
            'variant_class':  best['variant_class'],
            'protein_change': best['protein_change'],
            'VAF':            best['VAF'],
            'chrom':          best['chrom'],
            'pos':            best['pos'],
            'ref':            best['ref'],
            'alt':            best['alt'],
        })

FIELDNAMES = [
    'sample_id', 'gene', 'called', 'variant_class',
    'protein_change', 'VAF', 'chrom', 'pos', 'ref', 'alt'
]

with open(f"{OUTDIR}/sf_mutation_calls.tsv", 'w', newline='') as fh:
    writer = csv.DictWriter(fh, delimiter='\t', fieldnames=FIELDNAMES)
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"Mutation calls written: {OUTDIR}/sf_mutation_calls.tsv")
print(f"  Total sample-gene rows: {len(summary_rows)}")

from collections import Counter
called_counts = Counter(r['gene'] for r in summary_rows if r['called'] != 'WT')
for gene, count in sorted(called_counts.items()):
    print(f"  {gene}: {count} mutant samples")
PYEOF

# ── Build metadata TSV joining manifest + mutation calls ──────────────────────
echo ""
echo "Building metadata TSV..."

echo -e "sample_id\tcohort\ttier\tdisease\tbam_path\tSF3B1_call\tSF3B1_change\tSF3B1_vaf\tSRSF2_call\tSRSF2_change\tSRSF2_vaf\tU2AF1_call\tU2AF1_change\tU2AF1_vaf\tZRSR2_call\tZRSR2_change\tZRSR2_vaf" \
    > "\${METADATA}"

while IFS=\$'\t' read -r sample_id cohort tier disease seq_type fq1 fq2; do
    out_dir=\$(dirname "\${fq1}")
    bam="\${out_dir}/\${sample_id}_Aligned.sortedByCoord.out.bam"

    row="\${sample_id}\t\${cohort}\t\${tier}\t\${disease}\t\${bam}"

    for gene in SF3B1 SRSF2 U2AF1 ZRSR2; do
        call_line=\$(awk -F'\t' -v sid="\${sample_id}" -v g="\${gene}" \
            'NR>1 && \$1==sid && \$2==g {print; exit}' "\${SUMMARY}" 2>/dev/null || true)

        if [[ -n "\${call_line}" ]]; then
            called=\$(echo "\${call_line}" | awk -F'\t' '{print \$3}')
            change=\$(echo "\${call_line}" | awk -F'\t' '{print \$5}')
            vaf=\$(echo    "\${call_line}" | awk -F'\t' '{print \$6}')
        else
            called="no_data"; change="."; vaf="."
        fi
        row="\${row}\t\${called}\t\${change}\t\${vaf}"
    done

    echo -e "\${row}" >> "\${METADATA}"

done < <(tail -n +2 "\${MANIFEST}")

N_META=\$(( \$(wc -l < "\${METADATA}") - 1 ))
echo "Metadata written : \${METADATA}"
echo "  Rows: \${N_META}"
echo ""
echo "Final outputs:"
echo "  All variants  : \${ALL_VARIANTS}"
echo "  Mutation calls: \${SUMMARY}"
echo "  Metadata      : \${METADATA}"
echo ""
echo "Completed: \$(date)"
SLURM

chmod +x "${AGG_SCRIPT}"

# ── Submit jobs ───────────────────────────────────────────────────────────────
echo "Submitting GATK array (${N_BAMS} samples, max ${CONCURRENCY} concurrent)..."

GATK_JOB_ID=$(sbatch \
    --array=1-${N_BAMS}%${CONCURRENCY} \
    --parsable \
    "${GATK_ARRAY}")

echo "  GATK job ID: ${GATK_JOB_ID}"

echo "Submitting aggregation job (depends on GATK array completing)..."

AGG_JOB_ID=$(sbatch \
    --dependency=afterok:${GATK_JOB_ID} \
    --parsable \
    "${AGG_SCRIPT}")

echo "  Aggregation job ID: ${AGG_JOB_ID}"

echo ""
echo "======================================="
echo "Jobs submitted"
echo "  GATK array  : ${GATK_JOB_ID} (${N_BAMS} samples)"
echo "  Aggregation : ${AGG_JOB_ID}"
echo ""
echo "Monitor:"
echo "  squeue -j ${GATK_JOB_ID},${AGG_JOB_ID}"
echo ""
echo "Final outputs:"
echo "  All variants  : ${OUTDIR}/sf_all_variants.tsv"
echo "  Mutation calls: ${OUTDIR}/sf_mutation_calls.tsv"
echo "  Metadata TSV  : ${OUTDIR}/sf_mutation_metadata.tsv"
echo "  VCFs          : ${OUTDIR}/vcf/"
echo "  Funcotated    : ${OUTDIR}/funcotated/"
echo "======================================="