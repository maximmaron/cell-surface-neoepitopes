#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --output=logs/01_kallisto_index_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# ==============================================================================
# 01_build_kallisto_index.sh
#
# Steps:
#   1a. Pass 1 — filter merged_final.gtf by ref_gene_id (version-stripped)
#       → catches reference-anchored StringTie transcripts of surface genes
#   1b. Pass 2 — coordinate overlap for novel MSTRG loci lacking ref_gene_id
#       → builds surface gene BED from GENCODE v49, intersects with MSTRG loci,
#         rewrites gene_id/gene_name attributes with the matched surface gene
#       → catches novel isoforms of surface genes that StringTie assembled
#         de novo without anchoring to the reference transcript
#   1c. Merge pass 1 + pass 2 into surface_merged.gtf
#   2.  Extract transcript FASTA with gffread
#   3.  Build kallisto index
#   4.  Build t2g map (uses ref_gene_id where available, overlapping gene otherwise)
#
# Dependencies: gffread, kallisto, pybedtools (conda install pybedtools)
# ==============================================================================

set -euo pipefail

PIPELINE_DIR="/data1/abdelwao/maxim/splicing_pipeline/kallisto"
GTF="/data1/abdelwao/maxim/splicing_pipeline/stringtie/merge/merged_final.gtf"
GENCODE_GTF="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.gtf"
GENOME="/data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/GRCh38.primary_assembly.genome.fa"
REF_DIR="${PIPELINE_DIR}/reference"

mkdir -p "${REF_DIR}" logs

# ── Step 1a: Pass 1 — ref_gene_id matching (version-stripped) ────────────────
echo "[$(date)] Pass 1: filtering by ref_gene_id..."

python3 - <<PYEOF
import re

surface_ensg = set(
    line.strip() for line in open("${REF_DIR}/surface_ensg.txt") if line.strip()
)
print(f"Surface ENSG IDs loaded: {len(surface_ensg)}", flush=True)

n_written = 0
pass1_mstrg = set()  # MSTRG gene_ids already captured via ref_gene_id

with open("${GTF}") as fh, open("${REF_DIR}/pass1.gtf", "w") as out:
    for line in fh:
        if line.startswith("#"):
            out.write(line)
            continue
        ref_m = re.search(r'ref_gene_id "([^"]+)"', line)
        gid_m = re.search(r'gene_id "([^"]+)"',     line)
        if ref_m:
            ref_gene_id = ref_m.group(1).split(".")[0]  # strip version suffix
            if ref_gene_id in surface_ensg:
                out.write(line)
                n_written += 1
                if gid_m:
                    pass1_mstrg.add(gid_m.group(1))

print(f"Pass 1 lines written:        {n_written}", flush=True)
print(f"Pass 1 MSTRG gene IDs found: {len(pass1_mstrg)}", flush=True)

# Save pass1 MSTRG IDs so pass 2 can exclude them
with open("${REF_DIR}/pass1_mstrg_ids.txt", "w") as f:
    for gid in pass1_mstrg:
        f.write(gid + "\n")
PYEOF

echo "[$(date)] Pass 1 complete: ${REF_DIR}/pass1.gtf"
wc -l "${REF_DIR}/pass1.gtf"

# ── Step 1b: Pass 2 — coordinate overlap for novel MSTRG loci ────────────────
echo "[$(date)] Pass 2: coordinate-based recovery of novel MSTRG loci..."

python3 - <<PYEOF
import re
import pybedtools

# ── Build surface gene BED from GENCODE v49 (gene-level intervals) ────────────
surface_ensg = set(
    line.strip() for line in open("${REF_DIR}/surface_ensg.txt") if line.strip()
)
pass1_mstrg = set(
    line.strip() for line in open("${REF_DIR}/pass1_mstrg_ids.txt") if line.strip()
)

print("Building surface gene BED from GENCODE v49...", flush=True)
gene_bed_lines = []
gene_info = {}  # ensg_bare -> (gene_name, ensg_versioned)

with open("${GENCODE_GTF}") as fh:
    for line in fh:
        if line.startswith("#") or "\tgene\t" not in line:
            continue
        fields = line.strip().split("\t")
        chrom, start, end, strand = fields[0], fields[3], fields[4], fields[6]
        gid_m  = re.search(r'gene_id "([^"]+)"',   line)
        gsym_m = re.search(r'gene_name "([^"]+)"', line)
        if not gid_m:
            continue
        gid_versioned = gid_m.group(1)
        gid_bare      = gid_versioned.split(".")[0]
        gsym          = gsym_m.group(1) if gsym_m else gid_bare
        if gid_bare in surface_ensg:
            # BED: 0-based half-open
            gene_bed_lines.append(
                f"{chrom}\t{int(start)-1}\t{end}\t{gid_bare}\t0\t{strand}\t{gsym}"
            )
            gene_info[gid_bare] = (gsym, gid_versioned)

print(f"  Surface gene intervals: {len(gene_bed_lines)}", flush=True)
surface_bed = pybedtools.BedTool("\n".join(gene_bed_lines) + "\n", from_string=True)

# ── Build BED of MSTRG loci that have NO ref_gene_id (truly novel) ────────────
print("Extracting novel MSTRG transcript intervals...", flush=True)
novel_lines = []
novel_transcript_lines = {}  # transcript_id -> list of GTF lines
novel_tid_to_mstrg     = {}  # transcript_id -> MSTRG gene_id

with open("${GTF}") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        gid_m = re.search(r'gene_id "([^"]+)"', line)
        if not gid_m:
            continue
        gid = gid_m.group(1)
        # Skip if not MSTRG, already in pass1, or has a ref_gene_id
        if not gid.startswith("MSTRG"):
            continue
        if gid in pass1_mstrg:
            continue
        if 'ref_gene_id' in line:
            continue

        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        chrom, feat, start, end, strand = (
            fields[0], fields[2], fields[3], fields[4], fields[6]
        )

        # Use transcript-level records to define locus intervals for overlap
        # BED6 only — store gid separately in a dict to avoid column offset issues
        if feat == "transcript":
            tid_m = re.search(r'transcript_id "([^"]+)"', line)
            tid   = tid_m.group(1) if tid_m else None
            if tid:
                novel_lines.append(
                    f"{chrom}\t{int(start)-1}\t{end}\t{tid}\t0\t{strand}"
                )
                novel_transcript_lines[tid] = []
                novel_tid_to_mstrg[tid] = gid

        # Collect all lines per transcript for later writing
        tid_m = re.search(r'transcript_id "([^"]+)"', line)
        if tid_m and tid_m.group(1) in novel_transcript_lines:
            novel_transcript_lines[tid_m.group(1)].append(line)

print(f"  Novel MSTRG transcripts: {len(novel_lines)}", flush=True)

if not novel_lines:
    print("  No novel MSTRG transcripts to process — pass 2 output will be empty.")
    open("${REF_DIR}/pass2.gtf", "w").close()
    open("${REF_DIR}/pass2_novel_gene_map.tsv", "w").write(
        "transcript_id\tmstrg_gene_id\tmatched_ensg\tmatched_gene_symbol\n"
    )
else:
    novel_bed = pybedtools.BedTool(
        "\n".join(novel_lines) + "\n", from_string=True
    )

    # ── Intersect novel MSTRG transcripts with surface gene intervals ─────────
    # -s: require same strand
    # -f 0.5: novel transcript must overlap >= 50% of its length with gene body
    # (catches transcripts that extend slightly beyond gene boundaries)
    hits = novel_bed.intersect(surface_bed, s=True, f=0.5, wo=True)

    # Map each transcript to its best-overlapping surface gene
    # (longest overlap wins if multiple surface genes overlap)
    tx_to_gene = {}  # transcript_id -> (ensg_bare, gene_symbol, overlap_bp)
    for hit in hits:
        fields    = str(hit).strip().split("\t")
        # novel BED6:   fields[0-5]  = chrom,start,end,tid,0,strand
        # surface BED7: fields[6-12] = chrom,start,end,ensg,0,strand,gsym
        # overlap bp:   fields[13]
        tid       = fields[3]
        ensg_bare = fields[9]    # name field from surface_bed (col 4 of surface = index 9)
        gsym      = fields[12]   # gsym field from surface_bed (col 7 of surface = index 12)
        overlap   = int(fields[13])
        if tid not in tx_to_gene or overlap > tx_to_gene[tid][2]:
            tx_to_gene[tid] = (ensg_bare, gsym, overlap)

    print(f"  Novel MSTRG transcripts matched to surface genes: {len(tx_to_gene)}", flush=True)

    # ── Write pass2.gtf: rewrite gene_id and gene_name attributes ────────────
    # Keep MSTRG transcript_id intact (needed for kallisto quant matching)
    # but replace gene_id with ENSG and add ref_gene_id + gene_name
    n_written = 0
    with open("${REF_DIR}/pass2.gtf", "w") as out:
        for tid, (ensg, gsym, _) in tx_to_gene.items():
            for gtf_line in novel_transcript_lines.get(tid, []):
                # Replace gene_id "MSTRG.X" with the matched ENSG
                new_line = re.sub(
                    r'gene_id "MSTRG\.[^"]+"',
                    f'gene_id "{ensg}"',
                    gtf_line
                )
                # Add gene_name if not present
                if 'gene_name' not in new_line:
                    new_line = new_line.rstrip("\n") + f' gene_name "{gsym}";\n'
                # Add ref_gene_id to mark as coordinate-recovered
                if 'ref_gene_id' not in new_line:
                    new_line = new_line.rstrip("\n") + \
                        f' ref_gene_id "{ensg}"; novel_recovery "coordinate_overlap";\n'
                out.write(new_line)
                n_written += 1

    print(f"  Pass 2 lines written: {n_written}", flush=True)

    # Save mapping for provenance
    with open("${REF_DIR}/pass2_novel_gene_map.tsv", "w") as f:
        f.write("transcript_id\tmstrg_gene_id\tmatched_ensg\tmatched_gene_symbol\n")
        for tid, (ensg, gsym, _) in tx_to_gene.items():
            mstrg = novel_tid_to_mstrg.get(tid, "unknown")
            f.write(f"{tid}\t{mstrg}\t{ensg}\t{gsym}\n")

pybedtools.cleanup()
print("Pass 2 complete.", flush=True)
PYEOF

echo "[$(date)] Pass 2 complete: ${REF_DIR}/pass2.gtf"
wc -l "${REF_DIR}/pass2.gtf"
echo "[$(date)] Novel gene map: ${REF_DIR}/pass2_novel_gene_map.tsv"
wc -l "${REF_DIR}/pass2_novel_gene_map.tsv"

# ── Step 1c: Merge pass1 + pass2 ─────────────────────────────────────────────
echo "[$(date)] Merging pass1 + pass2 into surface_merged.gtf..."

# Carry over header from pass1, then append pass2 (no header duplication)
grep "^#" "${REF_DIR}/pass1.gtf"      >  "${REF_DIR}/surface_merged.gtf"
grep -v "^#" "${REF_DIR}/pass1.gtf"   >> "${REF_DIR}/surface_merged.gtf"
grep -v "^#" "${REF_DIR}/pass2.gtf"   >> "${REF_DIR}/surface_merged.gtf"

echo "[$(date)] surface_merged.gtf total lines:"
wc -l "${REF_DIR}/surface_merged.gtf"

PASS1_TX=$(grep -c $'\ttranscript\t' "${REF_DIR}/pass1.gtf" || echo 0)
PASS2_TX=$(grep -c $'\ttranscript\t' "${REF_DIR}/pass2.gtf" || echo 0)
echo "[$(date)] Transcripts — pass1 (ref_gene_id): ${PASS1_TX}  pass2 (novel overlap): ${PASS2_TX}"

# ── Step 2: Extract transcript FASTA with gffread ────────────────────────────
echo "[$(date)] Extracting transcript sequences with gffread..."

gffread "${REF_DIR}/surface_merged.gtf" \
    -g "${GENOME}" \
    -w "${REF_DIR}/surface_transcripts.fa" \
    --no-pseudo

N_TX=$(grep -c "^>" "${REF_DIR}/surface_transcripts.fa" || echo 0)
echo "[$(date)] Transcripts extracted: ${N_TX}"

if [[ "${N_TX}" -eq 0 ]]; then
    echo "ERROR: surface_transcripts.fa is empty. Check chromosome name" \
         "compatibility between GTF and genome FASTA." >&2
    exit 1
fi

# ── Step 3: Build kallisto index ─────────────────────────────────────────────
echo "[$(date)] Building kallisto index..."

kallisto index \
    -i "${REF_DIR}/surface_kallisto.idx" \
    -k 31 \
    "${REF_DIR}/surface_transcripts.fa"

echo "[$(date)] Kallisto index: ${REF_DIR}/surface_kallisto.idx"

# ── Step 4: Build t2g map ─────────────────────────────────────────────────────
# For pass1 transcripts: use ref_gene_id (ENSG, version-stripped) + gene_name
# For pass2 transcripts: gene_id has already been rewritten to ENSG
# For truly novel loci without a match: fall back to MSTRG gene_id
echo "[$(date)] Building t2g map..."

python3 - <<PYEOF
import re

gtf_path = "${REF_DIR}/surface_merged.gtf"
out_path = "${REF_DIR}/t2g.tsv"

t2g = {}
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith("#") or "\ttranscript\t" not in line:
            continue
        tid_m  = re.search(r'transcript_id "([^"]+)"', line)
        gid_m  = re.search(r'gene_id "([^"]+)"',       line)
        gsym_m = re.search(r'gene_name "([^"]+)"',     line)
        rgid_m = re.search(r'ref_gene_id "([^"]+)"',   line)
        if not (tid_m and gid_m):
            continue
        tid  = tid_m.group(1)
        gid  = gid_m.group(1)
        gsym = gsym_m.group(1) if gsym_m else gid

        # Prefer ref_gene_id (version-stripped) over MSTRG gene_id
        if rgid_m:
            gid = rgid_m.group(1).split(".")[0]

        t2g[tid] = (gid, gsym)

with open(out_path, "w") as out:
    out.write("transcript_id\tgene_id\tgene_symbol\n")
    for tid, (gid, gsym) in t2g.items():
        out.write(f"{tid}\t{gid}\t{gsym}\n")

print(f"t2g entries written: {len(t2g)}")
PYEOF

echo "[$(date)] t2g map: ${REF_DIR}/t2g.tsv"
echo "[$(date)] Done."