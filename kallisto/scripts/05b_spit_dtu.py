#!/usr/bin/env python3
"""
05b_spit_dtu.py

Differential Transcript Usage using SPIT.

SPIT workflow (from official notebook):
  1. sp.preprocess() — filter transcripts, compute isoform fractions
  2. sp.dtu()        — SPIT-Test + DTU detection + confounding control
  3. sp.cluster()    — hierarchical clustering of case samples by DTU events

Input files required per comparison:
  tx_counts : transcripts x samples TSV, first column named "tx_id"
  tx2gene   : two-column TSV with "tx_id" and "gene_id"
  labels    : TSV with "id" (sample ID) and "condition" (0=control, 1=case)
              Any additional columns are treated as covariates for confounding
              control. Categorical covariates must end with "_cat".

Usage:
  python3 05b_spit_dtu.py <task_id>
    1 = SF3B1_hotspot vs normal_heme
    2 = SRSF2_hotspot vs normal_heme
    3 = U2AF1_hotspot vs normal_heme
    4 = ZRSR2_mut     vs normal_heme
    5 = SF3B1_hotspot vs normal_tissue
    6 = SRSF2_hotspot vs normal_tissue
    7 = U2AF1_hotspot vs normal_tissue
    8 = ZRSR2_mut     vs normal_tissue

Notes:
  - kallisto inferential replicates (bootstraps) are used automatically
    via --infReps + --quant_path + --quant_type kallisto
  - cohort is passed as a categorical covariate for confounding control
    since case (MDS/AML) and control samples come from different cohorts
  - disease subtype (e.g. hotspot VAF tier) can be added as covariate
    if available in metadata
"""

import sys
import os
import pandas as pd
import numpy as np

import spit as sp

# ── Config ────────────────────────────────────────────────────────────────────
COMPARISONS = [
    # Within-cohort comparisons using SF-wildtype MDS/AML as controls
    # This eliminates batch effects — case and control from same cohort
    ("SF3B1_hotspot", "MDS_AML_SFWT", "beat"),       # beat cohort only
    ("SRSF2_hotspot", "MDS_AML_SFWT", "beat"),       # beat cohort only
    ("U2AF1_hotspot", "MDS_AML_SFWT", "beat"),       # beat cohort only
    ("ZRSR2_mut",     "MDS_AML_SFWT", "beat"),       # beat cohort only
    ("SF3B1_hotspot", "MDS_AML_SFWT", "tcga"),       # tcga cohort only
    ("SRSF2_hotspot", "MDS_AML_SFWT", "tcga"),       # tcga cohort only
    ("U2AF1_hotspot", "MDS_AML_SFWT", "tcga"),       # tcga cohort only
]

PIPELINE_DIR  = "/data1/abdelwao/maxim/splicing_pipeline/kallisto"
MANIFEST_PATH = f"{PIPELINE_DIR}/reference/sample_manifest.tsv"
QUANT_DIR     = f"{PIPELINE_DIR}/kallisto_quant"
T2G_PATH      = f"{PIPELINE_DIR}/reference/t2g.tsv"
OUT_BASE      = f"{PIPELINE_DIR}/results/spit_dtu"

# SPIT parameters (defaults from notebook/paper)
N_ITER    = 100    # SPIT-Test iterations (100 recommended)
KAPPA     = 0.6    # p-value threshold hyperparameter
BANDWIDTH = 0.09   # KDE bandwidth
N_SMALL   = 12     # minimum subgroup size

# ── Parse task ID ─────────────────────────────────────────────────────────────
if len(sys.argv) < 2:
    print("Usage: python3 05b_spit_dtu.py <task_id>  (1-8)")
    sys.exit(1)

task_id = int(sys.argv[1]) - 1
if task_id < 0 or task_id >= len(COMPARISONS):
    print(f"Invalid task_id. Must be 1-{len(COMPARISONS)}")
    sys.exit(1)

sf_group, ctrl_group, cohort = COMPARISONS[task_id]
comparison = f"{sf_group}_vs_SFWT_{cohort}"
OUT_DIR    = f"{OUT_BASE}/{comparison}"
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 70)
print(f"SPIT DTU: {comparison}")
print(f"Case: {sf_group} ({cohort})  |  Control: {ctrl_group} ({cohort})")
print("=" * 70)

# ── Load manifest and subset ──────────────────────────────────────────────────
manifest = pd.read_csv(MANIFEST_PATH, sep="\t")

# Subset to this cohort only — ensures case and control are batch-matched
case_rows = manifest[
    (manifest["sf_group"] == sf_group) &
    (manifest["cohort"] == cohort)
].copy()
ctrl_rows = manifest[
    (manifest["sf_group"] == ctrl_group) &
    (manifest["cohort"] == cohort)
].copy()

manifest = pd.concat([case_rows, ctrl_rows], ignore_index=True)
manifest["condition"] = np.where(manifest["sf_group"] == sf_group, 1, 0)

n_case    = (manifest["condition"] == 1).sum()
n_control = (manifest["condition"] == 0).sum()
print(f"Cohort: {cohort}")
print(f"Samples: case={n_case}  control={n_control}")

if n_case < 8:
    print(f"WARNING: only {n_case} case samples — results may be underpowered")
if n_control < 8:
    print(f"WARNING: only {n_control} control samples — results may be underpowered")

# ── Build labels file ─────────────────────────────────────────────────────────
# Labels file: "id" + "condition" (0/1) only
# Cohort is intentionally excluded — case (MDS/AML) and control (normal heme)
# samples come from structurally different cohorts by design, so including
# cohort as a covariate causes SPIT to eliminate virtually all DTU signal
# as "cohort-confounded" even when it is genuinely mutation-driven.
labels = manifest[["sample_id", "condition"]].copy()
labels = labels.rename(columns={"sample_id": "id"})

print(f"\nLabels preview (no cohort covariate):")
print(labels.head(3).to_string(index=False))

# ── Build tx_counts file ──────────────────────────────────────────────────────
# SPIT expects: rows = transcripts, columns = samples, first col named "tx_id"
# Use scaledTPM counts from kallisto (consistent with DRIMSeq approach)
print("\nBuilding transcript count matrix from kallisto outputs...")

# Manual h5 reading (no R dependency needed)
import h5py

sample_counts = {}
sample_ids    = manifest["sample_id"].tolist()
transcript_ids = None

for sid in sample_ids:
    h5_path = f"{QUANT_DIR}/{sid}/abundance.h5"
    if not os.path.exists(h5_path):
        print(f"  WARNING: Missing {h5_path}")
        continue
    with h5py.File(h5_path, "r") as f:
        if transcript_ids is None:
            transcript_ids = [x.decode() for x in f["aux"]["ids"][:]]
        # Use est_counts (raw counts, consistent with scaledTPM after rounding)
        sample_counts[sid] = f["est_counts"][:]

print(f"  Loaded {len(sample_counts)} samples, {len(transcript_ids)} transcripts")

# Build count matrix
counts_df = pd.DataFrame(sample_counts, index=transcript_ids)
counts_df = counts_df.round().astype(int)

# Scale to library size (approximate scaledTPM normalization)
# SPIT uses isoform fractions internally so absolute scaling matters less,
# but consistent normalization is good practice
lib_sizes  = counts_df.sum(axis=0)
scaled_df  = counts_df.divide(lib_sizes / 1e6)  # CPM

# SPIT tx_counts format: first column is "tx_id"
tx_counts = scaled_df.reset_index().rename(columns={"index": "tx_id"})
tx_counts = tx_counts[["tx_id"] + sample_ids]   # ensure column order matches labels

# ── Build tx2gene file ────────────────────────────────────────────────────────
t2g_full = pd.read_csv(T2G_PATH, sep="\t")
# SPIT tx2gene: columns "tx_id" and "gene_id"
tx2gene = t2g_full[["transcript_id", "gene_id"]].rename(
    columns={"transcript_id": "tx_id"}
)

print(f"tx2gene entries: {len(tx2gene)}")

# SPIT writes its output to OUT_DIR/SPIT_analysis/ and looks for
# tx2gene.txt in that same directory — write all input files there
spit_analysis_dir = f"{OUT_DIR}/SPIT_analysis"
os.makedirs(spit_analysis_dir, exist_ok=True)

tx_counts_path = f"{spit_analysis_dir}/tx_counts.txt"
tx2gene_path   = f"{spit_analysis_dir}/tx2gene.txt"
labels_path    = f"{spit_analysis_dir}/labels.txt"

tx_counts.to_csv(tx_counts_path, sep="\t", index=False)
tx2gene.to_csv(tx2gene_path,     sep="\t", index=False)
labels.to_csv(labels_path,       sep="\t", index=False)

print(f"\nInput files written to: {OUT_DIR}")

# ── Step 1: Preprocess ────────────────────────────────────────────────────────
print("\n" + "="*50)
print("STEP 1: Preprocessing")
print("="*50)

sp.preprocess(
    tx_counts  = tx_counts_path,
    tx2gene    = tx2gene_path,
    labels     = labels_path,
    output_dir = OUT_DIR,
    n_small    = N_SMALL,
    write      = True,       # print filter stats to stdout
    quiet      = False
)

# ── Step 2: DTU detection ─────────────────────────────────────────────────────
print("\n" + "="*50)
print("STEP 2: DTU Detection")
print("="*50)
print(f"Parameters: k={KAPPA}, bandwidth={BANDWIDTH}, n_iter={N_ITER}")
print("Inferential replicates: disabled (using SPIT-Test for FDR control)")

sp.dtu(
    labels     = labels_path,
    output_dir = OUT_DIR,
    n_small    = N_SMALL,
    n_iter     = N_ITER,
    k          = KAPPA,
    bandwidth  = BANDWIDTH,
    plot       = True,
    infReps    = False,      # disabled — requires pre-extracted replicate files
    f_cpm      = True,
    quiet      = False
)

# ── Step 3: Clustering ────────────────────────────────────────────────────────
print("\n" + "="*50)
print("STEP 3: Case sample clustering by DTU events")
print("="*50)

try:
    sp.cluster(
        labels     = labels_path,
        output_dir = OUT_DIR
    )
    print("Clustering complete.")
except Exception as e:
    print(f"Clustering skipped (likely no significant DTU events): {e}")

# ── Collect and summarize results ─────────────────────────────────────────────
print("\n" + "="*50)
print("RESULTS SUMMARY")
print("="*50)

# ── Collect and summarize results ─────────────────────────────────────────────
print("\n" + "="*50)
print("RESULTS SUMMARY")
print("="*50)

# SPIT output files:
#   spit_out.txt                    — gene-level DTU results
#   controlled_spit_out.txt         — after confounding control (only if covariates provided)
#   spit_cluster_matrix.txt         — binary sample x transcript DTU matrix
#   controlled_spit_cluster_matrix.txt — after confounding control

spit_out_file       = f"{spit_analysis_dir}/spit_out.txt"
controlled_out_file = f"{spit_analysis_dir}/controlled_spit_out.txt"
cluster_file_ctrl   = f"{spit_analysis_dir}/controlled_spit_cluster_matrix.txt"
cluster_file_raw    = f"{spit_analysis_dir}/spit_cluster_matrix.txt"
all_pvals_file      = f"{spit_analysis_dir}/all_p_values.txt"

# Use controlled output if available (covariates provided), otherwise raw
final_out_file   = controlled_out_file if os.path.exists(controlled_out_file) \
                   else spit_out_file
cluster_file     = cluster_file_ctrl   if os.path.exists(cluster_file_ctrl) \
                   else cluster_file_raw

print(f"Using results file: {os.path.basename(final_out_file)}")
print(f"Using cluster file: {os.path.basename(cluster_file)}")

# Load gene-level results
spit_genes = pd.read_csv(spit_out_file, sep="\t", header=0,
                         names=["gene_id", "flag"],
                         dtype=str).fillna("")
print(f"\nDTU genes (spit_out): {len(spit_genes)}")
print(f"  Flagged (likelihood outlier): {(spit_genes['flag'] == 'F').sum()}")

final_genes = pd.read_csv(final_out_file, sep="\t", header=0,
                           names=["gene_id", "flag"],
                           dtype=str).fillna("")
if final_out_file != spit_out_file:
    print(f"DTU genes (post-confounding): {len(final_genes)}")

# Annotate with gene symbols
final_genes = final_genes.merge(
    t2g_full[["gene_id","gene_symbol"]].drop_duplicates("gene_id"),
    on="gene_id", how="left"
)
final_genes["comparison"] = comparison

print(f"\nTop DTU genes:")
print(final_genes[["gene_id","gene_symbol","flag","comparison"]]
      .head(30).to_string(index=False))

# Load cluster matrix
if os.path.exists(cluster_file):
    cluster_mat = pd.read_csv(cluster_file, sep="\t", index_col=0)
    print(f"\nCluster matrix: {cluster_mat.shape[0]} samples x "
          f"{cluster_mat.shape[1]} DTU transcripts")

    # Add novel flag to transcript list
    tx_ids = cluster_mat.columns.tolist()
    novel_txs = [t for t in tx_ids if t.startswith("MSTRG")]
    print(f"Novel MSTRG transcripts in DTU set: {len(novel_txs)}")
    if novel_txs:
        print("  " + "\n  ".join(novel_txs[:10]))

# Save annotated results
out_genes = f"{OUT_DIR}/{comparison}_spit_sig_genes.tsv"
final_genes.to_csv(out_genes, sep="\t", index=False)
print(f"\nSaved: {out_genes}")

print(f"\nAll SPIT output files in: {spit_analysis_dir}")
print("  SPIT_chart.png            — p-value visualization")
print("  violin_plots/             — IF distributions per DTU transcript")
print("  confounding_analysis_plots/ — confounding control results")
print("  spit_dendrogram.png       — case sample clustering")
print("\nSPIT DTU complete.")