#!/usr/bin/env python3
"""
06_prioritize_candidates.py

Integrates DESeq2 transcript DE results with StringTie GTF annotation
to produce a tiered, ranked neoepitope candidate list.

Prioritization tiers (applied in order):
  Tier 1 — Novel MSTRG transcripts upregulated in mutant
            (peptide sequence absent from normal proteome)
  Tier 2 — Annotated transcripts with novel junctions (class code j/o)
            upregulated in mutant (novel peptide at junction)
  Tier 3 — Isoform switching genes (multiple upregulated transcripts
            from same gene — suggests broad isoform remodeling)
  Tier 4 — High-confidence DE surface transcripts (baseMean > 10,
            LFC > 2, padj < 0.01)

Additional filters applied:
  - Upregulated in mutant (LFC > 0)
  - baseMean >= 5 (detectable expression)
  - Normal tissue specificity filter (if vs_tissue results available:
    flag transcripts NOT upregulated in normal tissue comparison)

Outputs:
  results/prioritized/
    candidates_all.tsv              — all significant upregulated transcripts
    candidates_tier1_novel.tsv      — MSTRG novel transcripts
    candidates_tier2_novel_junc.tsv — novel junction transcripts
    candidates_tier3_isoswitch.tsv  — isoform switching genes
    candidates_top100.tsv           — top 100 ranked across all tiers
    igv_regions.txt                 — IGV-ready locus strings for top candidates
    igv_batch.txt                   — IGV batch script to screenshot all top hits
    summary.txt                     — human-readable summary
"""

import pandas as pd
import numpy as np
import re
import os
import glob

# ── Paths ─────────────────────────────────────────────────────────────────────
PIPELINE_DIR  = "/data1/abdelwao/maxim/splicing_pipeline/kallisto"
GTF_MERGED    = "/data1/abdelwao/maxim/splicing_pipeline/stringtie/merge/merged_final.gtf"
SURFACE_GENES = "/data1/abdelwao/maxim/splicing_pipeline/kallisto/reference/surface_genes.tsv"
OUT_DIR       = f"{PIPELINE_DIR}/results/prioritized"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
PADJ_THRESH    = 0.05
LFC_THRESH     = 1.0       # minimum log2FC for upregulation
LFC_HC_THRESH  = 2.0       # high-confidence threshold
BASEMEAN_MIN   = 10        # minimum expression
IGV_PADDING    = 2000      # bp padding around transcript for IGV view
TOP_N_IGV      = 100       # number of top candidates for IGV batch script
NOVEL_CODES    = {"j", "o", "u", "i", "x", "p", "e", "s", "c"}

print("=" * 70)
print("NEOEPITOPE CANDIDATE PRIORITIZATION")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load DESeq2 results
# ─────────────────────────────────────────────────────────────────────────────
deseq2_files = sorted(glob.glob(f"{PIPELINE_DIR}/results/deseq2/*_transcript_DE.tsv"))
if not deseq2_files:
    raise FileNotFoundError(f"No DESeq2 result files found in {PIPELINE_DIR}/results/deseq2/")

print(f"\nLoading {len(deseq2_files)} DESeq2 result files:")
dfs = []
for f in deseq2_files:
    df = pd.read_csv(f, sep="\t")
    print(f"  {os.path.basename(f)}: {len(df)} transcripts, "
          f"{(df['padj'] < PADJ_THRESH).sum()} sig")
    dfs.append(df)
deseq2 = pd.concat(dfs, ignore_index=True)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Filter to upregulated significant hits
# ─────────────────────────────────────────────────────────────────────────────
sig_up = deseq2[
    (deseq2["padj"] < PADJ_THRESH) &
    (deseq2["log2FoldChange"] > LFC_THRESH) &
    (deseq2["baseMean"] >= BASEMEAN_MIN)
].copy()

print(f"\nSignificant upregulated (padj<{PADJ_THRESH}, LFC>{LFC_THRESH}, "
      f"baseMean>={BASEMEAN_MIN}): {len(sig_up)}")
print(f"  Across {sig_up['comparison'].nunique()} comparisons:")
for comp, n in sig_up['comparison'].value_counts().items():
    print(f"    {comp}: {n}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Load DRIMSeq DTU results if available
# ─────────────────────────────────────────────────────────────────────────────
dtu_files = sorted(glob.glob(f"{PIPELINE_DIR}/results/drimseq_dtu/*_sig.tsv"))
dtu_sig_txs = set()
dtu_df = None
if dtu_files:
    print(f"\nLoading {len(dtu_files)} DRIMSeq significant result files:")
    dtu_dfs = []
    for f in dtu_files:
        df = pd.read_csv(f, sep="\t")
        if len(df) > 0:
            dtu_dfs.append(df)
            dtu_sig_txs.update(df["transcript_id"].tolist())
            print(f"  {os.path.basename(f)}: {len(df)} DTU transcripts")
    if dtu_dfs:
        dtu_df = pd.concat(dtu_dfs, ignore_index=True)
    print(f"  Total unique DTU transcripts: {len(dtu_sig_txs)}")
else:
    print("\nNo DRIMSeq results found — proceeding with DESeq2 only")

sig_up["has_dtu"] = sig_up["transcript_id"].isin(dtu_sig_txs)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Parse StringTie GTF for class codes and genomic coordinates
# ─────────────────────────────────────────────────────────────────────────────
print("\nParsing StringTie GTF for class codes and coordinates...")

# Only parse transcripts we care about — much faster
target_txs = set(sig_up["transcript_id"].unique())

tx_info = {}
with open(GTF_MERGED) as fh:
    for line in fh:
        if line.startswith("#") or "\ttranscript\t" not in line:
            continue
        tid_m = re.search(r'transcript_id "([^"]+)"', line)
        if not tid_m:
            continue
        tid = tid_m.group(1)
        if tid not in target_txs:
            continue

        fields = line.strip().split("\t")
        def attr(key):
            m = re.search(rf'{key} "([^"]+)"', line)
            return m.group(1) if m else None

        tx_info[tid] = {
            "chrom":          fields[0],
            "start":          int(fields[3]),
            "end":            int(fields[4]),
            "strand":         fields[6],
            "class_code":     attr("class_code"),
            "ref_transcript": attr("cmp_ref"),
            "ref_gene_id":    attr("ref_gene_id"),
        }

print(f"  Coordinates found for {len(tx_info)} / {len(target_txs)} transcripts")

# Add coordinate and class code info to sig_up
coord_df = pd.DataFrame.from_dict(tx_info, orient="index").reset_index()
coord_df.rename(columns={"index": "transcript_id"}, inplace=True)
sig_up = sig_up.merge(coord_df, on="transcript_id", how="left")

# Strip version suffix from ref_gene_id for matching
sig_up["ref_gene_id_bare"] = sig_up["ref_gene_id"].str.split(".").str[0]

# Novel flags
sig_up["is_novel_mstrg"]  = sig_up["transcript_id"].str.startswith("MSTRG")
sig_up["is_novel_junc"]   = sig_up["class_code"].isin(NOVEL_CODES) & \
                             ~sig_up["is_novel_mstrg"]
sig_up["is_novel_any"]    = sig_up["is_novel_mstrg"] | sig_up["is_novel_junc"]

# ─────────────────────────────────────────────────────────────────────────────
# 5. Load surface gene list
# ─────────────────────────────────────────────────────────────────────────────
surface = pd.read_csv(SURFACE_GENES, sep="\t")
surface_symbols = set(surface["gene_symbol"].dropna())
surface_ensg    = set(surface["ensg_id"].dropna())

sig_up["is_surface_confirmed"] = (
    sig_up["gene_symbol"].isin(surface_symbols) |
    sig_up["gene_id"].isin(surface_ensg) |
    sig_up["ref_gene_id_bare"].isin(surface_ensg)
)

# ─────────────────────────────────────────────────────────────────────────────
# 6. Isoform switching: genes with multiple upregulated transcripts
# ─────────────────────────────────────────────────────────────────────────────
# Count upregulated transcripts per gene per comparison
tx_per_gene = (
    sig_up.groupby(["gene_id", "comparison"])["transcript_id"]
    .count()
    .reset_index()
    .rename(columns={"transcript_id": "n_up_tx_in_comparison"})
)
sig_up = sig_up.merge(tx_per_gene, on=["gene_id", "comparison"], how="left")

# Across all comparisons: max upregulated transcripts per gene
max_tx_per_gene = (
    sig_up.groupby("gene_id")["n_up_tx_in_comparison"]
    .max()
    .reset_index()
    .rename(columns={"n_up_tx_in_comparison": "max_up_tx_per_gene"})
)
sig_up = sig_up.merge(max_tx_per_gene, on="gene_id", how="left")
sig_up["is_isoform_switch"] = sig_up["max_up_tx_per_gene"] >= 3

# ─────────────────────────────────────────────────────────────────────────────
# 7. Specificity: consistent across comparisons
# ─────────────────────────────────────────────────────────────────────────────
# Transcripts significant in multiple comparisons are more robust
n_sig_per_tx = (
    sig_up.groupby("transcript_id")["comparison"]
    .count()
    .reset_index()
    .rename(columns={"comparison": "n_sig_comparisons"})
)
sig_up = sig_up.merge(n_sig_per_tx, on="transcript_id", how="left")

# ─────────────────────────────────────────────────────────────────────────────
# 8. Composite priority score
# ─────────────────────────────────────────────────────────────────────────────
# Higher score = more attractive neoepitope candidate
# Weights reflect: novelty > DTU > consistency > magnitude

sig_up["priority_score"] = (
    sig_up["is_novel_mstrg"].astype(float)  * 5.0 +   # Novel MSTRG: highest
    sig_up["is_novel_junc"].astype(float)   * 3.0 +   # Novel junction
    sig_up["has_dtu"].astype(float)         * 2.0 +   # Isoform switching confirmed
    sig_up["is_isoform_switch"].astype(float) * 1.0 + # Gene-level switching
    np.log1p(sig_up["n_sig_comparisons"])   * 1.5 +   # Consistent across comparisons
    np.log1p(sig_up["log2FoldChange"])      * 1.0 +   # Magnitude
    np.log1p(sig_up["baseMean"])            * 0.5 -   # Expression level
    np.log1p(-np.log10(sig_up["padj"].clip(1e-300))) * 0.0  # (padj already filtered)
)

# ─────────────────────────────────────────────────────────────────────────────
# 9. Tier assignment
# ─────────────────────────────────────────────────────────────────────────────
def assign_tier(row):
    if row["is_novel_mstrg"]:
        return 1
    elif row["is_novel_junc"] and row["log2FoldChange"] >= LFC_HC_THRESH:
        return 2
    elif row["is_isoform_switch"] and row["log2FoldChange"] >= LFC_HC_THRESH:
        return 3
    else:
        return 4

sig_up["tier"] = sig_up.apply(assign_tier, axis=1)

tier_labels = {
    1: "Novel MSTRG transcript",
    2: "Novel junction (annotated gene)",
    3: "Isoform switching gene",
    4: "High-confidence DE surface transcript"
}
sig_up["tier_label"] = sig_up["tier"].map(tier_labels)

# ─────────────────────────────────────────────────────────────────────────────
# 10. IGV coordinates
# ─────────────────────────────────────────────────────────────────────────────
sig_up["igv_locus"] = sig_up.apply(
    lambda r: f"{r['chrom']}:{max(1, r['start'] - IGV_PADDING)}-{r['end'] + IGV_PADDING}"
    if pd.notna(r["chrom"]) else "NA",
    axis=1
)

# ─────────────────────────────────────────────────────────────────────────────
# 11. Save outputs
# ─────────────────────────────────────────────────────────────────────────────
# Column order for output tables
OUT_COLS = [
    "tier", "tier_label", "priority_score",
    "transcript_id", "gene_id", "gene_symbol",
    "chrom", "start", "end", "strand",
    "class_code", "ref_transcript", "ref_gene_id",
    "is_novel_mstrg", "is_novel_junc", "is_novel_any",
    "is_surface_confirmed", "is_isoform_switch", "has_dtu",
    "log2FoldChange", "baseMean", "padj",
    "n_sig_comparisons", "n_up_tx_in_comparison", "max_up_tx_per_gene",
    "comparison", "igv_locus"
]
OUT_COLS = [c for c in OUT_COLS if c in sig_up.columns]

# Sort by tier then priority score
sig_up_sorted = sig_up.sort_values(["tier", "priority_score"], ascending=[True, False])

# All candidates
sig_up_sorted[OUT_COLS].to_csv(f"{OUT_DIR}/candidates_all.tsv", sep="\t", index=False)

# Per-tier files
for tier, label in tier_labels.items():
    subset = sig_up_sorted[sig_up_sorted["tier"] == tier]
    fname = f"candidates_tier{tier}_{label.split()[0].lower()}.tsv"
    subset[OUT_COLS].to_csv(f"{OUT_DIR}/{fname}", sep="\t", index=False)
    print(f"\nTier {tier} ({label}): {len(subset)} transcripts from "
          f"{subset['gene_id'].nunique()} genes")
    if len(subset) > 0:
        print(f"  Top 5:")
        print(subset[["transcript_id","gene_symbol","log2FoldChange","baseMean","padj","comparison"]]
              .head(5).to_string(index=False))

# Top 100 overall
top100 = sig_up_sorted.drop_duplicates("transcript_id").head(TOP_N_IGV)
top100[OUT_COLS].to_csv(f"{OUT_DIR}/candidates_top100.tsv", sep="\t", index=False)

# IGV regions file (can paste directly into IGV search box or use Go to Locus)
igv_regions = top100[top100["igv_locus"] != "NA"]["igv_locus"].unique()
with open(f"{OUT_DIR}/igv_regions.txt", "w") as f:
    for locus in igv_regions:
        f.write(locus + "\n")

# IGV batch script — automates screenshots for all top hits
with open(f"{OUT_DIR}/igv_batch.txt", "w") as f:
    f.write("new\n")
    f.write("genome hg38\n")
    f.write("snapshotDirectory igv_snapshots\n\n")
    for _, row in top100[top100["igv_locus"] != "NA"].iterrows():
        gsym = row["gene_symbol"] if pd.notna(row["gene_symbol"]) else row["gene_id"]
        tid  = row["transcript_id"]
        f.write(f"goto {row['igv_locus']}\n")
        f.write(f"snapshot {gsym}_{tid}.png\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 12. Summary report
# ─────────────────────────────────────────────────────────────────────────────
summary_lines = [
    "NEOEPITOPE CANDIDATE PRIORITIZATION SUMMARY",
    "=" * 60,
    f"Input: {len(deseq2_files)} DESeq2 comparisons",
    f"Total transcripts tested: {len(deseq2):,}",
    f"Significant upregulated (padj<{PADJ_THRESH}, LFC>{LFC_THRESH}, "
    f"baseMean>={BASEMEAN_MIN}): {len(sig_up_sorted.drop_duplicates('transcript_id')):,}",
    "",
    "TIER BREAKDOWN (unique transcripts):",
]
for tier, label in tier_labels.items():
    n = len(sig_up_sorted[sig_up_sorted["tier"] == tier].drop_duplicates("transcript_id"))
    n_genes = sig_up_sorted[sig_up_sorted["tier"] == tier]["gene_id"].nunique()
    summary_lines.append(f"  Tier {tier} — {label}:")
    summary_lines.append(f"    {n:>6} transcripts from {n_genes} genes")

summary_lines += [
    "",
    "NOVELTY:",
    f"  MSTRG novel transcripts: "
    f"{sig_up_sorted['is_novel_mstrg'].sum()}",
    f"  Novel junction transcripts: "
    f"{sig_up_sorted['is_novel_junc'].sum()}",
    f"  DTU-confirmed (DRIMSeq): "
    f"{sig_up_sorted['has_dtu'].sum()}",
    "",
    "TOP 20 CANDIDATES BY PRIORITY SCORE:",
    sig_up_sorted.drop_duplicates("transcript_id").head(20)[
        ["tier","gene_symbol","transcript_id","log2FoldChange",
         "baseMean","padj","comparison","igv_locus"]
    ].to_string(index=False),
    "",
    "OUTPUT FILES:",
    f"  {OUT_DIR}/candidates_all.tsv",
    f"  {OUT_DIR}/candidates_tier1_novel.tsv",
    f"  {OUT_DIR}/candidates_tier2_novel.tsv",
    f"  {OUT_DIR}/candidates_tier3_isoform.tsv",
    f"  {OUT_DIR}/candidates_tier4_high.tsv",
    f"  {OUT_DIR}/candidates_top100.tsv",
    f"  {OUT_DIR}/igv_regions.txt",
    f"  {OUT_DIR}/igv_batch.txt",
]

summary = "\n".join(summary_lines)
print("\n" + summary)

with open(f"{OUT_DIR}/summary.txt", "w") as f:
    f.write(summary + "\n")

print(f"\nDone. All outputs in: {OUT_DIR}")