#!/bin/bash
# =============================================================================
# make_sf_gene_bed.sh
# Generate exon BED file for SF genes from GENCODE v49 using MANE Select
# transcripts only.
#
# Genes: SF3B1, SRSF2, U2AF1, ZRSR2
#
# Output:
#   sf_genes_exons.bed       — exon intervals (0-based, for variant calling)
#   sf_genes_exons_padded.bed — exons padded 10bp each side (catches splice sites)
#
# Usage: bash make_sf_gene_bed.sh
# =============================================================================

GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"
OUTDIR="/data1/abdelwao/maxim/splicing_pipeline/gatk"
mkdir -p "${OUTDIR}"

BED="${OUTDIR}/sf_genes_exons.bed"
PADDED="${OUTDIR}/sf_genes_exons_padded.bed"

echo "Extracting SF gene exons from GENCODE v49 (MANE Select)..."
echo "GTF: ${GTF}"
echo ""

# MANE Select transcript IDs for each SF gene (GRCh38 / GENCODE v49)
# These are the canonical clinical transcripts
# SF3B1 : ENST00000335508.10
# SRSF2 : ENST00000329501.8
# U2AF1 : ENST00000291552.8
# ZRSR2 : ENST00000370225.8

awk '
$3 == "exon" {
    # Extract transcript_id
    match($0, /transcript_id "([^"]+)"/, arr)
    tid = arr[1]

    # Extract gene_name
    match($0, /gene_name "([^"]+)"/, arr)
    gene = arr[1]

    # Check if MANE Select transcript for our genes of interest
    is_target = 0
    if (gene == "SF3B1" && tid ~ /ENST00000335508/) is_target = 1
    if (gene == "SRSF2" && tid ~ /ENST00000359995/) is_target = 1
    if (gene == "U2AF1" && tid ~ /ENST00000291552/) is_target = 1
    if (gene == "ZRSR2" && tid ~ /ENST00000307771/) is_target = 1

    if (is_target) {
        # Output 0-based BED format
        chrom = $1
        start = $4 - 1   # GTF is 1-based, BED is 0-based
        end   = $5
        strand = $7

        # Extract exon_number
        match($0, /exon_number ([0-9]+)/, arr)
        exon_num = arr[1]

        print chrom"\t"start"\t"end"\t"gene"_exon"exon_num"\t0\t"strand
    }
}' "${GTF}" | sort -k1,1 -k2,2n > "${BED}"

echo "Exon BED written: ${BED}"
echo "  Total exons: $(wc -l < ${BED})"
echo ""
echo "Per gene:"
awk '{split($4,a,"_"); print a[1]}' "${BED}" | sort | uniq -c | sort -rn
echo ""

# Padded BED — add 10bp each side to capture splice site variants
awk 'OFS="\t" {
    start = $2 - 10
    end   = $3 + 10
    if(start < 0) start = 0
    print $1, start, end, $4, $5, $6
}' "${BED}" | sort -k1,1 -k2,2n > "${PADDED}"

echo "Padded BED (±10bp): ${PADDED}"
echo "  Total intervals: $(wc -l < ${PADDED})"
echo ""

# Merge overlapping intervals for variant calling
# (some exons may overlap after padding)
echo "Merging overlapping intervals..."
MERGED="${OUTDIR}/sf_genes_exons_merged.bed"

sort -k1,1 -k2,2n "${PADDED}" | \
    awk 'OFS="\t" {
        if(NR==1) {
            chrom=$1; start=$2; end=$3; name=$4
        } else if($1==chrom && $2<=end) {
            if($3>end) end=$3
            name=name","$4
        } else {
            print chrom, start, end, name
            chrom=$1; start=$2; end=$3; name=$4
        }
    } END {
        print chrom, start, end, name
    }' > "${MERGED}"

echo "Merged BED: ${MERGED}"
echo "  Total intervals: $(wc -l < ${MERGED})"
echo ""

echo "Preview:"
cat "${BED}" | column -t
echo ""
echo "Done: $(date)"