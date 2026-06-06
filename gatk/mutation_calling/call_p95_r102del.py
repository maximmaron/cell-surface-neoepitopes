#!/usr/bin/env python3
"""
call_p95_r102del.py

Calls SRSF2 P95_R102del VAF from RNA-seq BAMs using the indel detection
logic from calculate_vaf.py. No SNV/codon calling is performed.

For each BAM listed in sf_final_metadata.tsv, scans the P95_R102 window
(chr17:76736805-76736929) for in-frame 24 bp deletions. STAR frequently
encodes this deletion as a novel CIGAR N rather than D; passing --gtf
enables detection of those reads.

Output: one .indel.tsv per sample in --outdir, then patch with
patch_p95_r102del.py.

Usage:
    python call_p95_r102del.py \\
        --metadata  sf_final_metadata.tsv \\
        --ref       /path/to/hg38.fa \\
        --outdir    ./del_calls \\
        [--gtf      /path/to/gencode.v44.annotation.gtf] \\
        [--min-mq   20] \\
        [--min-depth 10] \\
        [--min-indel-reads 3] \\
        [--min-del-size 3] \\
        [--max-del-size 500]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam is not installed.  Run:  pip install pysam")

# ---------------------------------------------------------------------------
# Locus constants (hg38, 0-based half-open)
# ---------------------------------------------------------------------------
CHROM        = "chr17"
REGION_START = 76736805   # window start (deletion start - 50 bp)
REGION_END   = 76736929   # window end   (deletion end   + 50 bp)
REGION_NAME  = "SRSF2_P95_R102del"
STRAND       = "-"

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate_sequence(seq):
    seq = seq.upper()
    return "".join(
        CODON_TABLE.get(seq[i:i+3], "?")
        for i in range(0, len(seq) - 2, 3)
    )


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--metadata",        required=True,
                   help="sf_final_metadata.tsv (must have sample_id and bam_path columns)")
    p.add_argument("--ref",             required=True,
                   help="Reference FASTA (must be .fai indexed)")
    p.add_argument("--outdir",          required=True,
                   help="Output directory for per-sample .indel.tsv files")
    p.add_argument("--gtf",             default=None,
                   help="GTF annotation. Required to detect deletions encoded "
                        "as novel CIGAR N by STAR (strongly recommended)")
    p.add_argument("--min-mq",          type=int, default=20)
    p.add_argument("--min-depth",       type=int, default=10)
    p.add_argument("--min-indel-reads", type=int, default=3)
    p.add_argument("--min-del-size",    type=int, default=3)
    p.add_argument("--max-del-size",    type=int, default=500)
    return p.parse_args()


# ---------------------------------------------------------------------------
# GTF junction loader (copied verbatim from calculate_vaf.py)
# ---------------------------------------------------------------------------
def load_known_junctions(gtf_path):
    import re as _re

    def _get_attr(attr_str):
        attrs = {}
        for part in attr_str.strip().split(";"):
            part = part.strip()
            if not part:
                continue
            m = (_re.match(r'(\S+)\s+"([^"]+)"', part)
                 or _re.match(r'(\S+)\s+(\S+)', part))
            if m:
                attrs[m.group(1)] = m.group(2)
        return attrs

    transcript_exons = defaultdict(list)
    print(f"Loading splice junctions from GTF: {gtf_path}", flush=True)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attrs = _get_attr(parts[8])
            tid   = attrs.get("transcript_id", "").split(".")[0]
            chrom = parts[0]
            start = int(parts[3]) - 1
            end   = int(parts[4])
            transcript_exons[(chrom, tid)].append((start, end))

    junctions = set()
    for (chrom, tid), exons in transcript_exons.items():
        for i, (_, e_end) in enumerate(sorted(exons)):
            if i + 1 < len(sorted(exons)):
                donor    = sorted(exons)[i][1]
                acceptor = sorted(exons)[i + 1][0]
                junctions.add((chrom, donor, acceptor))

    print(f"  Loaded {len(junctions):,} known splice junctions.", flush=True)
    return junctions


# ---------------------------------------------------------------------------
# Depth helper
# ---------------------------------------------------------------------------
def pileup_depth(bam, chrom, pos_0, min_mq):
    depth = 0
    for col in bam.pileup(chrom, pos_0, pos_0 + 1,
                          truncate=True, stepper="all",
                          min_mapping_quality=min_mq):
        if col.pos == pos_0:
            for pr in col.pileups:
                aln = pr.alignment
                if not aln.is_secondary and not aln.is_supplementary:
                    depth += 1
            break
    return depth


# ---------------------------------------------------------------------------
# Indel detection (copied verbatim from calculate_vaf.py, minus GTF load)
# ---------------------------------------------------------------------------
def detect_indels(bam, ref, chrom, start_0, end_0, region_name, strand,
                  min_mq, min_indel_reads, min_depth, known_junctions,
                  min_del_size, max_del_size):

    event_reads  = defaultdict(set)
    insert_seqs  = {}

    for aln in bam.fetch(chrom, start_0, end_0):
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        if aln.mapping_quality < min_mq:
            continue
        if aln.cigartuples is None:
            continue

        ref_pos   = aln.reference_start
        query_pos = 0

        for op, length in aln.cigartuples:
            if op in (0, 7, 8):
                ref_pos   += length
                query_pos += length
            elif op == 2:   # D
                pos = ref_pos
                if pos < end_0 and (pos + length) > start_0:
                    event_reads[("deletion", pos, length)].add(aln.query_name)
                ref_pos += length
            elif op == 3:   # N
                pos = ref_pos
                if known_junctions is not None:
                    acceptor = pos + length
                    if (chrom, pos, acceptor) not in known_junctions:
                        if pos < end_0 and (pos + length) > start_0:
                            event_reads[("deletion", pos, length)].add(aln.query_name)
                ref_pos += length
            elif op == 1:   # I
                pos = ref_pos
                if start_0 <= pos < end_0:
                    key = ("insertion", pos, length)
                    event_reads[key].add(aln.query_name)
                    if key not in insert_seqs and aln.query_sequence:
                        insert_seqs[key] = aln.query_sequence[query_pos:query_pos + length].upper()
                query_pos += length
            elif op == 4:
                query_pos += length

    rows = []
    mid  = (start_0 + end_0) // 2

    if not event_reads:
        depth = pileup_depth(bam, chrom, mid, min_mq)
        return [_empty_row(chrom, start_0, end_0, region_name, strand, depth)]

    for (event_type, pos, size), names in sorted(event_reads.items(),
                                                  key=lambda x: -len(x[1])):
        support = len(names)
        if support < min_indel_reads:
            continue
        if size < min_del_size or size > max_del_size:
            continue

        local_depth = pileup_depth(bam, chrom, pos, min_mq)
        if local_depth < min_depth:
            continue

        inframe = "yes" if size % 3 == 0 else "no"
        vaf     = round(support / local_depth, 6) if local_depth > 0 else 0.0

        if event_type == "deletion":
            try:
                seq_g = ref.fetch(chrom, pos, pos + size).upper()
            except Exception:
                seq_g = "N" * size
            seq_c   = reverse_complement(seq_g) if strand == "-" else seq_g
            alt_aas = translate_sequence(seq_c) if inframe == "yes" else "frameshift"
            row_start, row_end = pos + 1, pos + size
        else:
            seq_g   = insert_seqs.get(("insertion", pos, size), "N" * size)
            seq_c   = reverse_complement(seq_g) if strand == "-" else seq_g
            alt_aas = translate_sequence(seq_c) if inframe == "yes" else "frameshift"
            row_start = row_end = pos + 1

        rows.append({
            "chrom":        chrom,
            "region_start": start_0 + 1,
            "region_end":   end_0,
            "region_name":  region_name,
            "strand":       strand,
            "event_type":   event_type,
            "event_start":  row_start,
            "event_end":    row_end,
            "event_size":   size,
            "inframe":      inframe,
            "alt_seq":      seq_c,
            "alt_aas":      alt_aas,
            "total_depth":  local_depth,
            "alt_depth":    support,
            "vaf":          vaf,
        })

    if not rows:
        depth = pileup_depth(bam, chrom, mid, min_mq)
        return [_empty_row(chrom, start_0, end_0, region_name, strand, depth)]

    return rows


def _empty_row(chrom, start_0, end_0, region_name, strand, depth):
    return {
        "chrom": chrom, "region_start": start_0 + 1, "region_end": end_0,
        "region_name": region_name, "strand": strand,
        "event_type": ".", "event_start": ".", "event_end": ".",
        "event_size": ".", "inframe": ".", "alt_seq": ".", "alt_aas": ".",
        "total_depth": depth, "alt_depth": 0, "vaf": 0.0,
    }


FIELDNAMES = [
    "chrom", "region_start", "region_end", "region_name", "strand",
    "event_type", "event_start", "event_end", "event_size", "inframe",
    "alt_seq", "alt_aas", "total_depth", "alt_depth", "vaf",
]

# ---------------------------------------------------------------------------
# BAM index helper
# ---------------------------------------------------------------------------
def ensure_bam_index(bam_path):
    bam = Path(bam_path)
    indices = [bam.parent / (bam.name + ".bai"), bam.with_suffix(".bai"),
               bam.parent / (bam.name + ".csi")]
    if any(i.is_file() for i in indices):
        return
    print(f"  No index found — indexing {bam_path} ...", flush=True)
    try:
        pysam.index(bam_path)
    except pysam.SamtoolsError:
        sorted_tmp = str(bam.parent / f".{bam.name}.sorted.tmp")
        pysam.sort("-o", sorted_tmp, bam_path)
        os.replace(sorted_tmp, bam_path)
        pysam.index(bam_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args   = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not os.path.isfile(args.ref):
        sys.exit(f"ERROR: FASTA not found: {args.ref}")
    if not os.path.isfile(args.ref + ".fai"):
        sys.exit(f"ERROR: FASTA index not found: {args.ref}.fai\n"
                 f"       Run:  samtools faidx {args.ref}")

    # Load known junctions once for all samples
    known_junctions = None
    if args.gtf:
        if not os.path.isfile(args.gtf):
            sys.exit(f"ERROR: GTF not found: {args.gtf}")
        known_junctions = load_known_junctions(args.gtf)
    else:
        print("WARNING: --gtf not provided. Novel CIGAR N operations (STAR-encoded "
              "deletions) will NOT be detected. VAF may be underestimated.", flush=True)

    # Load sample list
    samples = []
    with open(args.metadata) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            bam = row.get("bam_path", ".")
            if bam and bam != ".":
                samples.append({"sample_id": row["sample_id"], "bam_path": bam})

    if not samples:
        sys.exit("ERROR: No samples with bam_path found in metadata.")
    print(f"Processing {len(samples)} samples ...", flush=True)

    ref = pysam.FastaFile(args.ref)

    for i, s in enumerate(samples, 1):
        sid      = s["sample_id"]
        bam_path = s["bam_path"]

        # Output filename matches calculate_vaf.py convention
        sample_name = Path(bam_path).stem.replace(".sorted", "")
        out_tsv     = outdir / f"{sample_name}.indel.tsv"

        if not Path(bam_path).is_file():
            print(f"  [{i}/{len(samples)}] SKIP (BAM not found): {bam_path}", flush=True)
            with open(out_tsv, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
                writer.writeheader()
                writer.writerow(_empty_row(CHROM, REGION_START, REGION_END,
                                           REGION_NAME, STRAND, "."))
            continue

        ensure_bam_index(bam_path)
        bam = pysam.AlignmentFile(bam_path, "rb")

        rows = detect_indels(
            bam, ref,
            chrom         = CHROM,
            start_0       = REGION_START,
            end_0         = REGION_END,
            region_name   = REGION_NAME,
            strand        = STRAND,
            min_mq        = args.min_mq,
            min_indel_reads = args.min_indel_reads,
            min_depth     = args.min_depth,
            known_junctions = known_junctions,
            min_del_size  = args.min_del_size,
            max_del_size  = args.max_del_size,
        )
        bam.close()

        with open(out_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        if i % 50 == 0 or i == len(samples):
            print(f"  {i}/{len(samples)} done", flush=True)

    ref.close()
    print(f"\nAll done. Results in: {outdir}", flush=True)


if __name__ == "__main__":
    main()