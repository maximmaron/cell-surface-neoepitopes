#!/usr/bin/env python3
"""
02_build_manifest.py

Generates sample_manifest.tsv for kallisto array by joining:
  - sf_final_metadata_with_del.tsv  (sample groups, mutation calls)
  - bam_summary.tsv                 (bam_path, strandedness, paired_end)

FASTQ files are in the same directory as the BAM files.
Naming conventions handled:
  {sample_id}_R1.fastq.gz / _R2
  {sample_id}_1.fastq.gz  / _2
  {sample_id}_R1_001.fastq.gz / _R2_001
  {sample_id}_R1.fq.gz    / _R2
  {sample_id}_1.fq.gz     / _2

Output columns:
  sample_id | fastq_r1 | fastq_r2 | strandedness | sf_group | cohort | condition
"""

import pandas as pd
import os
import glob

METADATA    = "/data1/abdelwao/maxim/splicing_pipeline/gatk/sf_final_metadata_with_del.tsv"
BAM_SUMMARY = "/data1/abdelwao/maxim/splicing_pipeline/rmats/bam_qc/bam_summary.tsv"  # adjust path if needed
OUT_DIR     = "/data1/abdelwao/maxim/splicing_pipeline/kallisto/reference"

# ── Load and join ─────────────────────────────────────────────────────────────
meta = pd.read_csv(METADATA, sep="\t")
bam  = pd.read_csv(BAM_SUMMARY, sep="\t")

print(f"Metadata rows:    {len(meta)}")
print(f"BAM summary rows: {len(bam)}")

# Rename bam_summary 'sample' to 'sample_id' for join if needed
if "sample" in bam.columns and "sample_id" not in bam.columns:
    bam = bam.rename(columns={"sample": "sample_id"})

# Join on sample_id — left join keeps all metadata samples
merged = meta.merge(
    bam[["sample_id", "bam_path", "stranded_call", "paired_end", "read_length_median"]],
    on="sample_id",
    how="left",
    suffixes=("_meta", "_bam")
)

# Prefer bam_summary bam_path over metadata bam_path (more reliable)
if "bam_path_bam" in merged.columns:
    merged["bam_path"] = merged["bam_path_bam"].fillna(merged.get("bam_path_meta", ""))
    merged = merged.drop(columns=["bam_path_bam", "bam_path_meta"], errors="ignore")

# Samples in metadata but missing from bam_summary
missing_bam = merged[merged["stranded_call"].isna()]
if len(missing_bam) > 0:
    print(f"\nWARNING: {len(missing_bam)} samples not in bam_summary — will default to unstranded:")
    print(missing_bam["sample_id"].tolist())
merged["stranded_call"] = merged["stranded_call"].fillna("unstranded")
merged["paired_end"]    = merged["paired_end"].fillna("yes")

# ── Map strandedness to kallisto flags ────────────────────────────────────────
# stranded_call values from RSeQC: 'reverse', 'forward', 'unstranded'
STRAND_MAP = {
    "reverse":    "fr-firststrand",   # dUTP / TruSeq stranded → --rf-stranded
    "forward":    "fr-secondstrand",  # ligation → --fr-stranded
    "unstranded": "unstranded",
}
merged["strandedness"] = merged["stranded_call"].str.lower().map(STRAND_MAP).fillna("unstranded")

unmapped = merged[~merged["stranded_call"].str.lower().isin(STRAND_MAP)]
if len(unmapped) > 0:
    print(f"\nWARNING: unrecognised stranded_call values (defaulting to unstranded):")
    print(unmapped["stranded_call"].unique())

# ── Derive FASTQ paths from BAM directory ─────────────────────────────────────
def get_fastq_paths(row):
    bam_path  = row["bam_path"]
    sample_id = row["sample_id"]

    if pd.isna(bam_path) or not isinstance(bam_path, str):
        return None, "NA"

    fastq_dir = os.path.dirname(bam_path)

    r1_patterns = [
        f"{sample_id}_R1.fastq.gz",
        f"{sample_id}_1.fastq.gz",
        f"{sample_id}_R1_001.fastq.gz",
        f"{sample_id}_R1.fq.gz",
        f"{sample_id}_1.fq.gz",
    ]
    r2_patterns = [
        f"{sample_id}_R2.fastq.gz",
        f"{sample_id}_2.fastq.gz",
        f"{sample_id}_R2_001.fastq.gz",
        f"{sample_id}_R2.fq.gz",
        f"{sample_id}_2.fq.gz",
    ]

    r1 = next(
        (os.path.join(fastq_dir, p) for p in r1_patterns
         if os.path.exists(os.path.join(fastq_dir, p))),
        None
    )
    r2 = next(
        (os.path.join(fastq_dir, p) for p in r2_patterns
         if os.path.exists(os.path.join(fastq_dir, p))),
        None
    )

    # Glob fallback
    if r1 is None:
        hits = sorted(glob.glob(os.path.join(fastq_dir, f"{sample_id}*_R1*")) +
                      glob.glob(os.path.join(fastq_dir, f"{sample_id}*_1.*")))
        hits = [h for h in hits if h.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))]
        r1 = hits[0] if hits else None

    if r2 is None:
        hits = sorted(glob.glob(os.path.join(fastq_dir, f"{sample_id}*_R2*")) +
                      glob.glob(os.path.join(fastq_dir, f"{sample_id}*_2.*")))
        hits = [h for h in hits if h.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))]
        r2 = hits[0] if hits else None

    return r1, (r2 if r2 else "NA")

print("\nFinding FASTQ files...")
paths = merged.apply(get_fastq_paths, axis=1).tolist()
merged["fastq_r1"] = [p[0] for p in paths]
merged["fastq_r2"] = [p[1] for p in paths]

# ── Report ────────────────────────────────────────────────────────────────────
missing_r1 = merged[merged["fastq_r1"].isna()]
if len(missing_r1) > 0:
    print(f"\nWARNING: {len(missing_r1)} samples with no R1 FASTQ found:")
    print(missing_r1[["sample_id", "bam_path"]].to_string())
else:
    print(f"All {len(merged)} R1 FASTQs found.")

n_pe = (merged["fastq_r2"] != "NA").sum()
print(f"Paired-end: {n_pe}  |  Single-end: {len(merged) - n_pe}")

print("\nStrandedness distribution:")
print(merged["strandedness"].value_counts().to_string())

print("\nExample paths (first 5):")
print(merged[["sample_id", "fastq_r1", "fastq_r2", "strandedness"]].head().to_string())

# ── Condition and sf_group from sample_group ──────────────────────────────────
# sf_group: direct passthrough, WT groups collapsed
# condition: mutant | wt_heme | wt_tissue
SFGROUP_MAP = {
    "SF3B1_hotspot":        "SF3B1_hotspot",
    "SF3B1_other":          "SF3B1_other",
    "SRSF2_hotspot":        "SRSF2_hotspot",
    "SRSF2_other":          "SRSF2_other",
    "U2AF1_hotspot":        "U2AF1_hotspot",
    "U2AF1_other":          "U2AF1_other",
    "ZRSR2_mut":            "ZRSR2_mut",
    "MDS_AML_SFWT":         "MDS_AML_SFWT",       # excluded from comparisons
    "normal_hematopoietic": "normal_heme",
    "normal_tissue":        "normal_tissue",
}
CONDITION_MAP = {
    "SF3B1_hotspot":        "mutant",
    "SF3B1_other":          "mutant",
    "SRSF2_hotspot":        "mutant",
    "SRSF2_other":          "mutant",
    "U2AF1_hotspot":        "mutant",
    "U2AF1_other":          "mutant",
    "ZRSR2_mut":            "mutant",
    "MDS_AML_SFWT":         "mds_aml_wt",         # excluded from comparisons
    "normal_hematopoietic": "wt_heme",
    "normal_tissue":        "wt_tissue",
}

merged["sf_group"]  = merged["sample_group"].map(SFGROUP_MAP).fillna("other")
merged["condition"] = merged["sample_group"].map(CONDITION_MAP).fillna("other")

unmapped = merged[merged["sf_group"] == "other"]["sample_group"].unique()
if len(unmapped) > 0:
    print(f"\nWARNING: unmapped sample_group values: {unmapped}")

print("\nsf_group distribution:")
print(merged["sf_group"].value_counts().to_string())
print("\ncondition distribution:")
print(merged["condition"].value_counts().to_string())

# ── Write manifest ────────────────────────────────────────────────────────────
manifest = merged[merged["fastq_r1"].notna()][[
    "sample_id", "fastq_r1", "fastq_r2",
    "strandedness", "sf_group", "cohort", "condition",
    "read_length_median"
]].copy()

os.makedirs(OUT_DIR, exist_ok=True)
out_path = f"{OUT_DIR}/sample_manifest.tsv"
manifest.to_csv(out_path, sep="\t", index=False)
print(f"\nManifest written: {out_path}  ({len(manifest)} samples)")
if len(missing_r1) > 0:
    print(f"  ({len(missing_r1)} samples excluded — no R1 FASTQ found)")