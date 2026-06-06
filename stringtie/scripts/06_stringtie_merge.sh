#!/bin/bash
# =============================================================================
# 06_stringtie_merge.sh
# Merge per-sample StringTie GTFs into a single unified GTF for rMATS.
#
# Three steps, all run fresh:
#   1. StringTie --merge: collapse redundant transcripts across all samples
#   1b. Contig filter: remove unplaced/haplotype scaffolds (GL*, KI*)
#   2. GFFCompare: classify merged transcripts against GENCODE v49
#   3. Python filter: exclude noise classes, write merged_final.gtf
#
# Exclusion logic (tightened from prior version):
#   Excluded regardless of exon count:
#     k — fully contained within reference intron (background transcription)
#     x — antisense overlap (not relevant for surface neoepitopes)
#     s — possible pre-mRNA fragment artifact
#   Excluded if single-exon only:
#     u — intergenic/unknown
#     o — generic overlap
#
#   Kept (including novel):
#     = — exact GENCODE match
#     c — contained in reference exon
#     j — novel junction, known gene  ← primary rMATS targets
#     i — intronic / retained intron  ← ZRSR2 targets
#     m/n — retained intron variants
#     e — single exon partial overlap (kept if multi-exon)
#     p — polymerase run-on (potential read-through neoepitopes)
#     y — contains reference within intron
#
# Usage: sbatch 06_stringtie_merge.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=stringtie_merge
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

source /admin/software/anaconda/1.11.1/bin/activate
conda activate stringtie

# ── Paths ─────────────────────────────────────────────────────────────────────
GTF_DIR="/data1/abdelwao/maxim/splicing_pipeline/stringtie/gtf"
OUTDIR="/data1/abdelwao/maxim/splicing_pipeline/stringtie/merge"
GENCODE_GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"
MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"

mkdir -p "${OUTDIR}"
mkdir -p /data1/abdelwao/maxim/splicing_pipeline/stringtie/logs

GTF_LIST="${OUTDIR}/gtf_list.txt"
MERGED_RAW="${OUTDIR}/merged.gtf"
MERGED_CANONICAL="${OUTDIR}/merged_canonical.gtf"
GFFCMP_PREFIX="${OUTDIR}/merged_compared"
MERGED_FINAL="${OUTDIR}/merged_final.gtf"
SUMMARY="${OUTDIR}/merge_summary.txt"

# Tee all output to a log file
exec > >(tee "${OUTDIR}/merge_log.txt") 2>&1

echo "======================================="
echo "StringTie merge — $(date)"
echo "GTF dir  : ${GTF_DIR}"
echo "Output   : ${OUTDIR}"
echo "GENCODE  : ${GENCODE_GTF}"
echo "======================================="
echo ""

# ── Build GTF list from manifest ──────────────────────────────────────────────
echo "── Building GTF list ───────────────────────────────────"

> "${GTF_LIST}"
FOUND=0
MISSING=0

while IFS=$'\t' read -r sample_id cohort tier disease seq_type fq1 fq2; do
    gtf="${GTF_DIR}/${sample_id}.gtf"
    if [[ -f "${gtf}" && -s "${gtf}" ]]; then
        echo "${gtf}" >> "${GTF_LIST}"
        FOUND=$(( FOUND + 1 ))
    else
        echo "  WARNING: Missing GTF for ${sample_id} (${cohort})"
        MISSING=$(( MISSING + 1 ))
    fi
done < <(tail -n +2 "${MANIFEST}")

echo "  GTFs found   : ${FOUND}"
echo "  GTFs missing : ${MISSING}"
echo ""

if [[ "${FOUND}" -eq 0 ]]; then
    echo "ERROR: No GTF files found in ${GTF_DIR}" >&2
    exit 1
fi

if [[ "${MISSING}" -gt 0 ]]; then
    echo "  WARNING: ${MISSING} samples missing GTFs — proceeding with ${FOUND}"
    echo "  Resubmit 05_stringtie_per_sample.sh for missing samples"
    echo ""
fi

# ── Step 1: StringTie --merge ─────────────────────────────────────────────────
echo "── Step 1: StringTie --merge ───────────────────────────"
echo "  Input GTFs : ${FOUND}"
echo "  Filters:"
echo "    -m 200   minimum transcript length"
echo "    -c 1.5   minimum coverage"
echo "    -F 1.0   minimum FPKM"
echo "    -T 1.0   minimum TPM"
echo "    -f 0.01  minimum isoform fraction"
echo "    -i       keep retained intron transcripts"
echo ""

stringtie --merge \
    -G "${GENCODE_GTF}" \
    -o "${MERGED_RAW}" \
    -m 200 \
    -c 1.5 \
    -F 1.0 \
    -T 1.0 \
    -f 0.01 \
    -i \
    -l MSTRG \
    "${GTF_LIST}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: StringTie merge failed" >&2
    exit 1
fi

N_MERGED=$(awk '$3=="transcript"' "${MERGED_RAW}" | wc -l)
echo "  Raw merged transcripts: ${N_MERGED}"
echo ""

# ── Step 1b: Remove unplaced/haplotype contigs ────────────────────────────────
# Keep only canonical chromosomes: chr1-22, chrX, chrY, chrM
# Removes GL*, KI* unplaced scaffolds which add noise and are not interpretable
echo "── Step 1b: Filtering unplaced contigs ─────────────────"

awk '$1 ~ /^chr([0-9]+|X|Y|M)$/ || $0 ~ /^#/' "${MERGED_RAW}" \
    > "${MERGED_CANONICAL}"

N_CANONICAL=$(awk '$3=="transcript"' "${MERGED_CANONICAL}" | wc -l)
N_REMOVED=$(( N_MERGED - N_CANONICAL ))
echo "  Transcripts removed (unplaced contigs) : ${N_REMOVED}"
echo "  Transcripts on canonical chromosomes   : ${N_CANONICAL}"
echo ""

# ── Step 2: GFFCompare ────────────────────────────────────────────────────────
echo "── Step 2: GFFCompare annotation ──────────────────────"

gffcompare \
    -r "${GENCODE_GTF}" \
    -G \
    -o "${GFFCMP_PREFIX}" \
    "${MERGED_CANONICAL}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: GFFCompare failed" >&2
    exit 1
fi

echo "  GFFCompare output: ${GFFCMP_PREFIX}.*"
echo ""

# ── Step 3: Python filtering ──────────────────────────────────────────────────
echo "── Step 3: Filtering to final GTF ──────────────────────"

python3 - << PYEOF
import re
from collections import defaultdict

annotated_gtf = "${GFFCMP_PREFIX}.annotated.gtf"
merged_input  = "${MERGED_CANONICAL}"
merged_final  = "${MERGED_FINAL}"
summary_file  = "${SUMMARY}"

# ── Parse GFFCompare annotated GTF ───────────────────────────────────────────
class_codes = {}
exon_counts  = defaultdict(int)

with open(annotated_gtf) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        feature = parts[2]
        attrs   = parts[8]

        tid_m = re.search(r'transcript_id "([^"]+)"', attrs)
        if not tid_m:
            continue
        tid = tid_m.group(1)

        if feature == 'transcript':
            cc_m = re.search(r'class_code "([^"]+)"', attrs)
            class_codes[tid] = cc_m.group(1) if cc_m else 'u'
        elif feature == 'exon':
            exon_counts[tid] += 1

# ── Exclusion logic ───────────────────────────────────────────────────────────
# Excluded regardless of exon count:
#   k — fully intronic (background transcription, not real isoforms)
#   x — antisense overlap (irrelevant for surface neoepitopes)
#   s — pre-mRNA fragment artifact
EXCLUDE_ALL = {'k', 'x', 's'}

# Excluded only if single-exon:
#   u — intergenic/unknown (multi-exon u may be real novel genes, keep)
#   o — generic overlap (multi-exon o may be real alternative isoforms, keep)
EXCLUDE_SINGLE_EXON = {'u', 'o'}

kept     = set()
excluded = set()
class_summary = defaultdict(int)

for tid, cc in class_codes.items():
    n_exons = exon_counts.get(tid, 0)
    class_summary[cc] += 1

    if cc in EXCLUDE_ALL:
        excluded.add(tid)
    elif cc in EXCLUDE_SINGLE_EXON and n_exons <= 1:
        excluded.add(tid)
    else:
        kept.add(tid)

# ── Write final GTF ───────────────────────────────────────────────────────────
with open(merged_input) as fin, open(merged_final, 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        tid_m = re.search(r'transcript_id "([^"]+)"', parts[8])
        if not tid_m:
            fout.write(line)
            continue
        if tid_m.group(1) not in excluded:
            fout.write(line)

# ── Write summary ─────────────────────────────────────────────────────────────
code_desc = {
    '=': 'Exact GENCODE match',
    'c': 'Contained in reference exon',
    'j': 'Novel junction, known gene  ← primary rMATS targets',
    'e': 'Single exon partial overlap',
    'i': 'Intronic (potential retained intron)  ← ZRSR2 targets',
    'k': 'Fully intronic (background transcription) [EXCLUDED]',
    'm': 'Retained intron — all sites match reference',
    'n': 'Retained intron — not all sites match reference',
    'o': 'Generic overlap',
    'p': 'Polymerase run-on (potential read-through neoepitopes)',
    's': 'Pre-mRNA fragment artifact [EXCLUDED]',
    'u': 'Intergenic/unknown',
    'x': 'Antisense overlap [EXCLUDED]',
    'y': 'Contains a reference within its intron',
}

with open(summary_file, 'w') as fh:
    fh.write("StringTie merge summary\n")
    fh.write("=" * 50 + "\n\n")
    fh.write(f"Total transcripts (canonical) : {len(class_codes)}\n")
    fh.write(f"Kept in final GTF             : {len(kept)}\n")
    fh.write(f"Excluded                      : {len(excluded)}\n\n")
    fh.write("GFFCompare class code distribution:\n")
    fh.write("-" * 50 + "\n")
    for cc in sorted(class_summary.keys()):
        fh.write(f"  {cc}  {class_summary[cc]:>8}  {code_desc.get(cc, 'other')}\n")

print(f"Final GTF written : {merged_final}")
print(f"  Kept            : {len(kept)}")
print(f"  Excluded        : {len(excluded)}")
print(f"\nClass code distribution:")
for cc in sorted(class_summary.keys()):
    print(f"  {cc}  {class_summary[cc]:>8}  {code_desc.get(cc, 'other')}")
PYEOF

if [[ $? -ne 0 ]]; then
    echo "ERROR: Python filtering failed" >&2
    exit 1
fi

# ── Final summary ─────────────────────────────────────────────────────────────
echo ""
echo "── Final summary ───────────────────────────────────────"
N_FINAL=$(awk '$3=="transcript"' "${MERGED_FINAL}" | wc -l)
N_FINAL_ENST=$(awk '$3=="transcript"' "${MERGED_FINAL}" | grep -c 'transcript_id "ENST' || true)
N_FINAL_NOVEL=$(awk '$3=="transcript"' "${MERGED_FINAL}" | grep -c 'transcript_id "MSTRG' || true)

echo "  merged.gtf (raw)       : ${N_MERGED} transcripts"
echo "  merged_canonical.gtf   : ${N_CANONICAL} transcripts (after contig filter)"
echo "  merged_final.gtf       : ${N_FINAL} transcripts (after class filter)"
echo "    GENCODE (ENST*)      : ${N_FINAL_ENST}"
echo "    Novel   (MSTRG*)     : ${N_FINAL_NOVEL}"
echo ""
cat "${SUMMARY}"
echo ""
echo "Use for rMATS --gtf:"
echo "  ${MERGED_FINAL}"
echo ""
echo "Completed: $(date)"
echo "======================================="

conda deactivate