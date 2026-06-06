#!/bin/bash
# =============================================================================
# 02_pool_junctions.sh
# Single SLURM job: archive, pool and filter SJ.out.tab from all pass 1 samples
#
# Steps:
#   1. Archive all pass 1 SJ.out.tab files to permanent location
#   2. Pool all samples with uniform read thresholds
#   3. Filter and produce final junction file for pass 2 index
#
# Pooling strategy: uniform across all tiers
#   - This paper is helpful reference: PMID: 25123659
#   - All samples pooled — maximum junction discovery
#   - Thresholds set by genomic dinucleotide frequency:
#     GT/AG, CT/AC (motifs 1,2) — common   : ≥3 reads OR ≥2 samples+≥2 reads
#     GC/AG, CT/GC (motifs 3,4) — moderate : ≥2 reads
#     AT/AC, GT/AT (motifs 5,6) — rare     : ≥1 read
#   - Minimum overhang: 6bp (matches --alignSJoverhangMin 6 in pass 1)
#   - Minimum intron length: 20bp (matches --alignIntronMin 20 in pass 1)
#
# SJ.out.tab columns:
#   1: chr
#   2: intron start (1-based)
#   3: intron end (1-based)
#   4: strand (0=undefined, 1=+, 2=-)
#   5: intron motif (0=non-canon, 1=GT/AG, 2=CT/AC, 3=GC/AG,
#                    4=CT/GC, 5=AT/AC, 6=GT/AT)
#   6: annotated junction (0=novel, 1=annotated)
#   7: number of uniquely mapping reads
#   8: number of multimapping reads
#   9: maximum spliced alignment overhang
#
# Outputs:
#   junctions/v1_{DATE}/pooled_raw.tab
#   junctions/v1_{DATE}/final_filtered_junctions.tab  ← for pass 2 index
#   junctions/pass1_sj_archive/{SAMPLE}_SJ.out.tab    ← permanent archive
#
# Usage: sbatch 02_pool_junctions.sh
# =============================================================================
#SBATCH -p cpu
#SBATCH --job-name=pool_junctions
#SBATCH --output=logs/pool_junctions_%j.out
#SBATCH --error=logs/pool_junctions_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

MANIFEST="/data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv"
JUNC_BASE="/data1/abdelwao/maxim/splicing_pipeline/junctions"
ARCHIVE_DIR="${JUNC_BASE}/pass1_sj_archive"

# ── Versioned output directory ────────────────────────────────────────────────
DATE=$(date +%Y%m%d)
JUNC_DIR="${JUNC_BASE}/v1_${DATE}"
mkdir -p "${JUNC_DIR}" "${ARCHIVE_DIR}" logs

echo "======================================="
echo "Junction pooling — $(date)"
echo "Output  : ${JUNC_DIR}"
echo "Archive : ${ARCHIVE_DIR}"
echo "Threads : ${SLURM_CPUS_PER_TASK}"
echo "======================================="
echo ""

# ── Step 1: Collect SJ.out.tab paths + archive ───────────────────────────────
echo "── Step 1: Collecting + archiving SJ.out.tab files ────"

SJ_LIST="${JUNC_DIR}/all_sj_paths.txt"
> "${SJ_LIST}"

FOUND=0
MISSING=0
ARCHIVED=0

while IFS=$'\t' read -r sample_id cohort tier disease seq_type fq1 fq2; do
    out_dir=$(dirname "${fq1}")
    sj="${out_dir}/${sample_id}_SJ.out.tab"

    if [[ ! -f "${sj}" ]]; then
        echo "  WARNING: Missing SJ.out.tab for ${sample_id}"
        ((MISSING++))
        continue
    fi

    # Add to pool list
    echo "${sj}" >> "${SJ_LIST}"
    ((FOUND++))

    # Archive — copy to permanent location if not already there
    archive_dest="${ARCHIVE_DIR}/${sample_id}_SJ.out.tab"
    if [[ ! -f "${archive_dest}" ]]; then
        cp "${sj}" "${archive_dest}"
        ((ARCHIVED++))
    fi

done < <(tail -n +2 "${MANIFEST}")

echo "  Found   : ${FOUND}"
echo "  Missing : ${MISSING}"
echo "  Archived: ${ARCHIVED} new files"
echo "  Archive : ${ARCHIVE_DIR} ($(ls ${ARCHIVE_DIR} | wc -l) total files)"
echo "  Size    : $(du -sh ${ARCHIVE_DIR} | cut -f1)"

if [[ "${MISSING}" -gt 0 ]]; then
    echo ""
    echo "WARNING: ${MISSING} SJ.out.tab files missing."
    echo "Set ABORT_ON_MISSING=1 to stop here."
    ABORT_ON_MISSING=0
    if [[ "${ABORT_ON_MISSING}" -eq 1 ]]; then
        echo "Aborting — rerun pass 1 for missing samples first." >&2
        exit 1
    fi
fi

if [[ "${FOUND}" -eq 0 ]]; then
    echo "ERROR: No SJ.out.tab files found — cannot pool." >&2
    exit 1
fi

# ── Step 2: Pool all samples ──────────────────────────────────────────────────
echo ""
echo "── Step 2: Pooling ${FOUND} samples ────────────────────"

POOLED_RAW="${JUNC_DIR}/pooled_raw.tab"

# Use xargs cat to handle large file lists safely
# avoids argument length limits from cat $(cat file_list)
# xargs splits the list into safe-sized batches automatically
echo "  Concatenating and aggregating SJ files..."
echo "  (using xargs to handle ${FOUND} files safely)"

xargs cat < "${SJ_LIST}" | \
    awk 'OFS="\t" {
        key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6
        reads[key]  += $7
        samples[key]++
        if($9 > overhang[key]) overhang[key] = $9
    }
    END {
        for(k in reads)
            print k, reads[k], 0, overhang[k], samples[k]
    }' | \
    sort --parallel="${SLURM_CPUS_PER_TASK}" -k1,1 -k2,2n \
    > "${POOLED_RAW}"

echo "  Unique junctions (before filter): $(wc -l < "${POOLED_RAW}")"

# ── Step 3: Filter ────────────────────────────────────────────────────────────
echo ""
echo "── Step 3: Filtering ───────────────────────────────────"
echo "  Thresholds (by genomic dinucleotide frequency):"
echo "    GT/AG, CT/AC (motifs 1,2, common)  : ≥3 reads OR ≥2 samples+≥2 reads"
echo "    GC/AG, CT/GC (motifs 3,4, moderate): ≥2 reads"
echo "    AT/AC, GT/AT (motifs 5,6, rare)    : ≥1 read"
echo "    Minimum overhang: 6bp"
echo "    Minimum intron : 20bp"
echo "    Exclude: non-standard chroms, non-canonical (motif 0)"
echo "    Novel junctions only (annotated already in GENCODE index)"

FINAL="${JUNC_DIR}/final_filtered_junctions.tab"

awk 'OFS="\t" {
    chr=$1; start=$2; end=$3; strand=$4
    motif=$5; annot=$6; reads=$7; overhang=$9

    # Exclude non-standard chromosomes — keep only chr1-22, chrX, chrY
    if(chr !~ /^chr([0-9]+|X|Y)$/) next

    # Novel junctions only — annotated ones already in GENCODE index
    if(annot != 0) next

    # Exclude unrecognized motifs (motif 0)
    if(motif == 0) next

    # Minimum intron length — matches --alignIntronMin 20 in pass 1
    if(end - start < 20) next

    # Minimum overhang — matches --alignSJoverhangMin 6 in pass 1
    if(overhang < 6) next

    # ── Motif-specific thresholds ──────────────────────────────────────────
    #
    # Rationale: threshold is determined by genomic frequency of the
    # splice site dinucleotide. Rare dinucleotides have low false positive
    # risk so single-read evidence is meaningful. Common dinucleotides
    # require more evidence to distinguish real junctions from artifacts.
    #
    # AT/AC (5) — U12 minor spliceosome, canonical ZRSR2 target
    #   Dinucleotide AT...AC is very rare in genome → ≥1 read sufficient
    #   ~700 known U12 introns in human genome; single-read events in
    #   ZRSR2 mutants are biologically meaningful
    #
    # GT/AT (6) — rare alternative splice site
    #   Dinucleotide GT...AT is rare in genome → ≥1 read sufficient
    #   Similar false positive profile to AT/AC based on read distribution
    #   (15,670 single-read events vs 445,745 for GC/AG — 30x less)
    #
    # GC/AG (3) and CT/GC (4) — U2-type variant splice sites
    #   Dinucleotide GC...AG is moderately common → ≥2 reads required
    #   445,745 single-read GC/AG events (66% of total) — high noise risk
    #   CT/GC is reverse complement of GC/AG, same reasoning applies
    #
    # GT/AG (1) and CT/AC (2) — major spliceosome
    #   Most common splice site — highest false positive risk
    #   Require ≥3 pooled reads OR cross-sample corroboration
    #   CT/AC is reverse complement of GT/AG on minus strand
    # ────────────────────────────────────────────────────────────────────

    # AT/AC (5) and GT/AT (6): ≥1 read
    if((motif == 5 || motif == 6) && reads >= 1) {
        print chr, start, end, strand, motif, annot
        next
    }

    # GC/AG (3) and CT/GC (4): ≥2 reads
    if((motif == 3 || motif == 4) && reads >= 2) {
        print chr, start, end, strand, motif, annot
        next
    }

    # GT/AG (1) and CT/AC (2): ≥3 reads OR ≥2 samples with ≥2 reads
    if((motif == 1 || motif == 2) &&        (reads >= 3 || (samples >= 2 && reads >= 2))) {
        print chr, start, end, strand, motif, annot
    }

}' "${POOLED_RAW}" > "${FINAL}"

echo "  Final junctions: $(wc -l < "${FINAL}")"

# ── Step 4: Summary ───────────────────────────────────────────────────────────
echo ""
echo "── Summary ─────────────────────────────────────────────"
echo "  Samples pooled    : ${FOUND}"
echo "  Raw junctions     : $(wc -l < "${POOLED_RAW}")"
echo "  Final junctions   : $(wc -l < "${FINAL}")"
echo ""

echo "  Motif breakdown (final):"
awk '{motif[$5]++} END {
    print "    GT/AG  (1) :", motif[1]+0, "(major spliceosome)"
    print "    CT/AC  (2) :", motif[2]+0, "(major spliceosome, minus strand)"
    print "    GC/AG  (3) :", motif[3]+0, "(minor variant)"
    print "    CT/GC  (4) :", motif[4]+0, "(minor variant, minus strand)"
    print "    AT/AC  (5) :", motif[5]+0, "(minor spliceosome, U12)"
    print "    GT/AT  (6) :", motif[6]+0, "(rare alternative)"
}' "${FINAL}"

echo ""
echo "  Sample support distribution (pooled raw):"
awk '{print $NF}' "${POOLED_RAW}" | \
    awk '{
        if($1==1) s1++
        else if($1<=5) s5++
        else if($1<=10) s10++
        else if($1<=50) s50++
        else multi++
    } END {
        print "    1 sample       :", s1+0
        print "    2-5 samples    :", s5+0
        print "    6-10 samples   :", s10+0
        print "    11-50 samples  :", s50+0
        print "    >50 samples    :", multi+0
    }'

echo ""
echo "  Read support distribution (pooled raw):"
awk '{print $7}' "${POOLED_RAW}" | \
    awk '{
        if($1<3) r1++
        else if($1<10) r10++
        else if($1<100) r100++
        else if($1<1000) r1000++
        else rhigh++
    } END {
        print "    <3 reads    :", r1+0, "(filtered out for major motifs)"
        print "    3-9 reads   :", r10+0
        print "    10-99 reads :", r100+0
        print "    100-999     :", r1000+0
        print "    ≥1000 reads :", rhigh+0
    }'

echo ""
echo "Archive : ${ARCHIVE_DIR}"
echo "  Files : $(ls ${ARCHIVE_DIR} | wc -l)"
echo "  Size  : $(du -sh ${ARCHIVE_DIR} | cut -f1)"
echo ""
echo "Final junction file: ${FINAL}"
echo ""
echo "Next: sbatch 03_build_pass2_index.sh"
echo "======================================="
echo "Completed: $(date)"