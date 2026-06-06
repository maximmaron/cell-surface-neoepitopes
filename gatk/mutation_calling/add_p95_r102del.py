#!/usr/bin/env python3
"""
patch_p95_r102del.py

Integrates SRSF2 P95_R102del calls from call_p95_r102del.py's .indel.tsv
outputs into sf_final_metadata.tsv.

Usage:
    python patch_p95_r102del.py \
        --metadata   sf_final_metadata.tsv \
        --indel-dir  ./del_calls \
        --outfile    sf_final_metadata_patched.tsv \
        [--min-vaf   0.02] \
        [--min-reads 3]
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

EXPECTED_SIZE   = 24
EXPECTED_REGION = "SRSF2_P95_R102del"

LOCKED_GROUPS = {
    "SRSF2_hotspot", "SRSF2_other",
    "SF3B1_hotspot", "SF3B1_other",
    "U2AF1_hotspot", "U2AF1_other",
    "ZRSR2_mut",
}

# Suffixes that call_p95_r102del.py may append between sample_id and .indel.tsv
# (derived from STAR BAM naming conventions)
STRIP_SUFFIXES = [
    "_Aligned.sortedByCoord.out",
    "_Aligned.sortedByCoord",
    "_AlignedByCoord.out",
    "_AlignedByCoord",
    ".sortedByCoord.out",
    ".sortedByCoord",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata",   required=True)
    p.add_argument("--indel-dir",  required=True)
    p.add_argument("--outfile",    required=True)
    p.add_argument("--min-vaf",    type=float, default=0.02)
    p.add_argument("--min-reads",  type=int,   default=3)
    return p.parse_args()


def find_del_call(indel_tsv: Path, min_vaf: float, min_reads: int):
    best = None
    try:
        with open(indel_tsv) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("region_name") != EXPECTED_REGION:
                    continue
                if row.get("event_type")  != "deletion":
                    continue
                if row.get("inframe")     != "yes":
                    continue
                try:
                    size      = int(row["event_size"])
                    alt_depth = int(row["alt_depth"])
                    vaf       = float(row["vaf"])
                except (ValueError, KeyError):
                    continue
                if size      != EXPECTED_SIZE:
                    continue
                if alt_depth <  min_reads:
                    continue
                if vaf       <  min_vaf:
                    continue
                if best is None or vaf > float(best["vaf"]):
                    best = row
    except FileNotFoundError:
        pass
    return best


def main():
    args      = parse_args()
    indel_dir = Path(args.indel_dir)

    # ── Load metadata ─────────────────────────────────────────────────────
    rows = []
    with open(args.metadata) as fh:
        reader         = csv.DictReader(fh, delimiter="\t")
        original_fields = list(reader.fieldnames)
        for row in reader:
            rows.append(row)
    print(f"Loaded {len(rows)} samples from {args.metadata}")

    # ── Index indel TSVs ──────────────────────────────────────────────────
    # Build two lookups:
    #   full_stem  → path   (e.g. "BA2006R_AlignedByCoord.out" → path)
    #   short_stem → path   (e.g. "BA2006R" → path, after stripping known suffixes)
    all_tsvs = list(indel_dir.glob("*.indel.tsv"))
    print(f"Found {len(all_tsvs)} .indel.tsv files in {indel_dir}")

    by_full  = {}   # full stem (everything before .indel.tsv)
    by_short = {}   # stem after stripping STAR suffixes

    for f in all_tsvs:
        full = f.name.replace(".indel.tsv", "")
        by_full[full] = f
        short = full
        for suf in STRIP_SUFFIXES:
            if short.endswith(suf):
                short = short[: -len(suf)]
                break
        by_short[short] = f

    new_field  = "SRSF2_P95_R102del_vaf"
    out_fields = original_fields + [new_field]

    reclassified   = 0
    skipped_locked = 0
    not_found      = 0

    for row in rows:
        sid        = row["sample_id"]
        row[new_field] = "."

        # Match: try exact full stem first, then short stem
        indel_tsv = by_full.get(sid) or by_short.get(sid)
        if indel_tsv is None:
            not_found += 1
            continue

        del_row = find_del_call(indel_tsv, args.min_vaf, args.min_reads)
        if del_row is None:
            continue

        row[new_field] = del_row["vaf"]

        current_group = row.get("sample_group", "")
        if current_group in LOCKED_GROUPS:
            skipped_locked += 1
            print(f"  SKIP (locked as {current_group}): {sid}")
            continue

        vaf = del_row["vaf"]
        row["sample_group"]   = "SRSF2_hotspot"
        row["mutated_gene"]   = "SRSF2"
        row["protein_change"] = "p.P95_R102del"
        row["vaf"]            = vaf
        row["SRSF2_call"]     = "SRSF2_mut"
        row["SRSF2_change"]   = "p.P95_R102del"
        row["SRSF2_vaf"]      = vaf
        reclassified += 1
        print(f"  RECLASSIFIED → SRSF2_hotspot: {sid}  "
              f"(vaf={vaf}, alt_depth={del_row['alt_depth']}, "
              f"total_depth={del_row['total_depth']})")

    # ── Write output ──────────────────────────────────────────────────────
    with open(args.outfile, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(r["sample_group"] for r in rows)
    print(f"\nSummary:")
    print(f"  Total samples            : {len(rows)}")
    print(f"  Reclassified             : {reclassified}")
    print(f"  Skipped (locked group)   : {skipped_locked}")
    print(f"  No indel TSV found       : {not_found}")
    print(f"  Output written to        : {args.outfile}")
    print(f"\nFinal sample_group counts:")
    for grp, cnt in sorted(counts.items()):
        print(f"  {grp:<25} {cnt}")


if __name__ == "__main__":
    main()