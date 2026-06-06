#!/usr/bin/env python3
"""
00_build_surface_gene_list.py

Build a unified surface protein gene list (ENSG IDs + gene symbols) from:
  1. SURFY (Surfaceome Label == 'surface')
  2. CSPA (all entries — already curated surface proteins)
  3. UniProt GO:0009986 reviewed

Outputs:
  reference/surface_genes.tsv  — columns: ensg_id, gene_symbol, source
  reference/surface_ensg.txt   — plain list of ENSG IDs for GTF filtering
"""

import pandas as pd
import re
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANNOT_DIR = "/data1/abdelwao/maxim/splicing_pipeline/kallisto/reference"

# ---------- 1. SURFY ----------
surfy = pd.read_excel(
    f"{ANNOT_DIR}/SURFY_table_S3_surfaceome.xlsx", header=1
)
surfy_surface = surfy[surfy["Surfaceome Label"] == "surface"].copy()
surfy_genes = surfy_surface[["Ensembl gene", "UniProt gene"]].rename(
    columns={"Ensembl gene": "ensg_id", "UniProt gene": "gene_symbol"}
)
surfy_genes["source"] = "SURFY"
surfy_genes = surfy_genes.dropna(subset=["ensg_id"])
# Some rows have semicolon-separated ENSG IDs — explode
surfy_genes["ensg_id"] = surfy_genes["ensg_id"].astype(str).str.split(";")
surfy_genes = surfy_genes.explode("ensg_id")
surfy_genes["ensg_id"] = surfy_genes["ensg_id"].str.strip()
surfy_genes = surfy_genes[surfy_genes["ensg_id"].str.startswith("ENSG")]

# ---------- 2. CSPA ----------
cspa = pd.read_excel(
    f"{ANNOT_DIR}/CSPA_validated_surfaceome_proteins.xlsx", sheet_name="Table A"
)
# CSPA has UniProt accession as ID_link — need to map to ENSG via SURFY lookup
# Also build a gene symbol set directly
cspa_genes = cspa[["ENTREZ gene symbol"]].rename(
    columns={"ENTREZ gene symbol": "gene_symbol"}
)
cspa_genes["ensg_id"] = None
cspa_genes["source"] = "CSPA"
# Enrich CSPA with ENSG from SURFY where possible
surfy_sym2ensg = surfy[["UniProt gene", "Ensembl gene"]].dropna()
surfy_sym2ensg.columns = ["gene_symbol", "ensg_id"]
surfy_sym2ensg["ensg_id"] = surfy_sym2ensg["ensg_id"].astype(str).str.split(";")
surfy_sym2ensg = surfy_sym2ensg.explode("ensg_id")
surfy_sym2ensg["ensg_id"] = surfy_sym2ensg["ensg_id"].str.strip()
sym2ensg = dict(zip(surfy_sym2ensg["gene_symbol"], surfy_sym2ensg["ensg_id"]))
cspa_genes["ensg_id"] = cspa_genes["gene_symbol"].map(sym2ensg)
cspa_genes = cspa_genes.dropna(subset=["gene_symbol"])

# ---------- 3. UniProt GO:0009986 ----------
uniprot = pd.read_csv(
    f"{ANNOT_DIR}/uniprotkb_go_0009986_AND_reviewed_true_2025_07_26.tsv.txt",
    sep="\t"
)
# Gene Names column: first token is primary symbol
uniprot["gene_symbol"] = uniprot["Gene Names"].astype(str).str.split().str[0]
uniprot_genes = uniprot[["gene_symbol"]].copy()
uniprot_genes["ensg_id"] = uniprot_genes["gene_symbol"].map(sym2ensg)
uniprot_genes["source"] = "UniProt_GO0009986"

# ---------- Merge ----------
all_genes = pd.concat([surfy_genes, cspa_genes, uniprot_genes], ignore_index=True)
all_genes = all_genes.dropna(subset=["gene_symbol"])
all_genes["ensg_id"] = all_genes["ensg_id"].fillna("NA")
all_genes["gene_symbol"] = all_genes["gene_symbol"].str.strip()

# Deduplicate: keep all sources per gene for provenance
all_genes = all_genes.drop_duplicates(subset=["ensg_id", "gene_symbol", "source"])

# Pivot source to single row per gene
def agg_sources(grp):
    return pd.Series({
        "gene_symbol": grp["gene_symbol"].iloc[0],
        "ensg_id": grp["ensg_id"].iloc[0],
        "sources": ";".join(sorted(grp["source"].unique()))
    })

merged = (
    all_genes.groupby("gene_symbol", sort=False)
    .apply(agg_sources)
    .reset_index(drop=True)
)

os.makedirs(f"{BASE}/reference", exist_ok=True)
merged.to_csv(f"{BASE}/reference/surface_genes.tsv", sep="\t", index=False)

# Write ENSG list (drop NA)
ensg_list = merged[merged["ensg_id"] != "NA"]["ensg_id"].unique()
with open(f"{BASE}/reference/surface_ensg.txt", "w") as f:
    for e in ensg_list:
        f.write(e + "\n")

print(f"Total surface genes: {len(merged)}")
print(f"  With ENSG ID:      {len(ensg_list)}")
print(f"  SURFY only:        {(merged['sources'] == 'SURFY').sum()}")
print(f"  CSPA only:         {(merged['sources'] == 'CSPA').sum()}")
print(f"  Multi-source:      {merged['sources'].str.contains(';').sum()}")
print(f"\nOutput: {BASE}/reference/surface_genes.tsv")