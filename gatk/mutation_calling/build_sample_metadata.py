#!/usr/bin/env python3
"""
build_final_metadata.py

Builds a final metadata sheet combining:
  - SF-mutant samples with granular hotspot classification:
      SRSF2_hotspot    : P95 missense (p.P95*), P95_R102 deletions, and R94dup
      SRSF2_other      : all other SRSF2 missense
      SF3B1_hotspot    : R625, K666, K700 missense mutations
      SF3B1_other      : all other SF3B1 missense
      U2AF1_hotspot    : S34 or Q157 missense mutations
      U2AF1_other      : all other U2AF1 missense
      ZRSR2_mut        : any ZRSR2 nonsense/frameshift
  - SF-wildtype MDS/AML samples       → labelled 'MDS_AML_SFWT'
  - Normal hematopoietic samples      → labelled 'normal_hematopoietic'
  - Normal tissue samples             → labelled 'normal_tissue'

Multi-mutation priority rules (applied when a sample has >1 qualifying SF hit):
  1. Any hotspot beats any non-hotspot → assign to the hotspot group
  2. Two hotspots → assign to the higher-VAF hotspot; flag co_mutation
  3. Two non-hotspots → assign to the higher-VAF non-hotspot; flag co_mutation
  All co-mutations are recorded in a 'co_mutation' column for transparency.

Usage:
    python3 build_final_metadata.py \
        --manifest  sample_manifest.tsv \
        --metadata  sf_mutation_metadata.tsv \
        --variants  sf_all_variants.tsv \
        --outfile   sf_final_metadata.tsv
"""

import re
import csv
import argparse
from collections import defaultdict

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--manifest', required=True)
parser.add_argument('--metadata', required=True)
parser.add_argument('--variants', required=True)
parser.add_argument('--outfile',  required=True)
args = parser.parse_args()

# ── Variant class definitions ─────────────────────────────────────────────────
MISSENSE = {'MISSENSE'}
NONSENSE = {'NONSENSE', 'NONSTOP', 'FRAME_SHIFT_INS', 'FRAME_SHIFT_DEL'}

# ── Disease label mapping ─────────────────────────────────────────────────────
DISEASE_LABELS = {
    'MDS_AML':              'MDS_AML_SFWT',
    'normal_hematopoietic': 'normal_hematopoietic',
    'normal_tissue':        'normal_tissue',
}

# ── Group tier: lower number = higher priority when VAF is tied ───────────────
# Used only as tiebreaker; VAF is the primary comparator
GROUP_TIER = {
    'SRSF2_hotspot': 1,   # canonical hotspots always beat ZRSR2_mut
    'SF3B1_hotspot': 1,
    'U2AF1_hotspot': 1,
    'ZRSR2_mut':     2,   # tier 2 — loses to hotspots, beats non-hotspot missense
    'SRSF2_other':   3,
    'SF3B1_other':   3,
    'U2AF1_other':   3,
}

HOTSPOT_GROUPS = {'SRSF2_hotspot', 'SF3B1_hotspot', 'U2AF1_hotspot', 'ZRSR2_mut'}

# ── Hotspot classifiers ───────────────────────────────────────────────────────
def classify_srsf2(protein_change, variant_classes):
    pc = protein_change or ''
    # Indel hotspots — check before MISSENSE guard
    if (re.search(r'p\.P95_R102del', pc, re.IGNORECASE)
            or re.search(r'p\.R94dup',    pc, re.IGNORECASE)
            or re.search(r'p\.94_95insR', pc, re.IGNORECASE)):
        return 'SRSF2_hotspot'
    if not (variant_classes & MISSENSE):
        return None
    if re.search(r'p\.P95[A-Z]', pc):
        return 'SRSF2_hotspot'
    return 'SRSF2_other'


def classify_sf3b1(protein_change, variant_classes):
    if not (variant_classes & MISSENSE):
        return None
    pc = protein_change or ''
    if re.search(r'p\.(R625|K666|K700)[A-Z]', pc):
        return 'SF3B1_hotspot'
    return 'SF3B1_other'


def classify_u2af1(protein_change, variant_classes):
    if not (variant_classes & MISSENSE):
        return None
    pc = protein_change or ''
    if re.search(r'p\.(S34|Q157)[A-Z]', pc):
        return 'U2AF1_hotspot'
    return 'U2AF1_other'


def classify_zrsr2(variant_classes):
    """
    ZRSR2_mut : any nonsense/frameshift present in variant_classes.
    The MISSENSE exclusion is removed — low-depth missense calls are
    already filtered upstream by MIN_DEPTH=10 in the aggregation step.
    """
    if variant_classes & NONSENSE:
        return 'ZRSR2_mut'
    return None


# ── Load manifest ─────────────────────────────────────────────────────────────
print(f"Loading manifest from {args.manifest}...")
manifest = {}
with open(args.manifest) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        manifest[row['sample_id']] = row
print(f"  {len(manifest)} samples in manifest")

# ── Load variant classes ──────────────────────────────────────────────────────
print(f"Loading variant classes from {args.variants}...")
variant_classes = defaultdict(set)
all_variants    = defaultdict(list)   # (sample_id, gene) -> list of full rows
with open(args.variants) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        key = (row['sample_id'], row['gene'].upper())
        variant_classes[key].add(row['variant_class'].upper())
        all_variants[key].append({
            'variant_class':  row['variant_class'].upper(),
            'protein_change': row.get('protein_change', '.'),
            'vaf':            row.get('VAF', '.'),
        })
print(f"  {len(variant_classes)} sample-gene entries loaded")

# ── Load mutation metadata ────────────────────────────────────────────────────
print(f"Loading mutation metadata from {args.metadata}...")
metadata = {}
with open(args.metadata) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        metadata[row['sample_id']] = row
print(f"  {len(metadata)} samples in metadata")


# ── Collect ALL qualifying SF hits per sample ─────────────────────────────────
def get_all_sf_hits(sample_id, meta_row):
    """
    Returns a list of (group, gene, protein_change, vaf_float) for every
    qualifying SF mutation found in this sample. List may be empty (no hits)
    or contain multiple entries (co-mutations).
    """
    hits = []

    classifiers = {
        'SRSF2': classify_srsf2,
        'SF3B1': classify_sf3b1,
        'U2AF1': classify_u2af1,
    }

    for gene, classifier in classifiers.items():
        call = meta_row.get(f'{gene}_call', '')
        pc   = meta_row.get(f'{gene}_change', '')
        vaf  = meta_row.get(f'{gene}_vaf', '.')
        if call != f'{gene}_mut':
            continue
        classes = variant_classes.get((sample_id, gene), set())
        group = classifier(pc, classes)
        if group:
            try:
                vaf_f = float(vaf)
            except (ValueError, TypeError):
                vaf_f = 0.0
            hits.append((group, gene, pc, vaf, vaf_f))

    # ZRSR2 — require nonsense/frameshift, VAF > 0.5.
    # The aggregation script may have selected a missense as the primary call
    # (e.g. high-VAF recurrent artifact on chrX) even when a truncating
    # mutation is also present. Here we look directly at variant_classes:
    # if any truncating class exists for this sample, we use that call
    # regardless of what ZRSR2_change says.
    call = meta_row.get('ZRSR2_call', '')
    pc   = meta_row.get('ZRSR2_change', '')
    vaf  = meta_row.get('ZRSR2_vaf', '.')
    if call == 'ZRSR2_mut':
        classes = variant_classes.get((sample_id, 'ZRSR2'), set())
        if classes & NONSENSE:
            # A truncating variant exists — find its protein change and VAF
            # directly from sf_all_variants rather than trusting ZRSR2_change
            trunc = None
            for v in all_variants.get((sample_id, 'ZRSR2'), []):
                if v['variant_class'].upper() in NONSENSE:
                    try:
                        v_vaf = float(v['vaf'])
                    except (ValueError, TypeError):
                        v_vaf = 0.0
                    if trunc is None or v_vaf > float(trunc['vaf']):
                        trunc = v
            if trunc and float(trunc['vaf']) > 0.5:
                hits.append(('ZRSR2_mut', 'ZRSR2',
                              trunc['protein_change'], trunc['vaf'],
                              float(trunc['vaf'])))
        else:
            # No truncating variant — fall back to metadata call with VAF gate
            try:
                vaf_f = float(vaf)
            except (ValueError, TypeError):
                vaf_f = 0.0
            if vaf_f > 0.5:
                group = classify_zrsr2(classes)
                if group:
                    hits.append((group, 'ZRSR2', pc, vaf, vaf_f))

    return hits


def resolve_hits(hits):
    """
    Given a list of hits from get_all_sf_hits(), apply priority rules and
    return (winning_group, gene, protein_change, vaf, co_mutation_str).

    Priority rules:
      1. Any hotspot beats any non-hotspot
      2. Among equals (both hotspot or both non-hotspot): higher VAF wins
      3. Ties in VAF: lower GROUP_TIER wins (SRSF2>SF3B1>U2AF1 for hotspots)

    co_mutation_str is '.' for single hits, or a semicolon-separated list of
    the non-winning mutations (gene:protein_change:VAF) for multi-hit samples.
    """
    if not hits:
        return None

    if len(hits) == 1:
        group, gene, pc, vaf, _ = hits[0]
        return (group, gene, pc, vaf, '.')

    hotspots    = [h for h in hits if h[0] in HOTSPOT_GROUPS]
    non_hotspot = [h for h in hits if h[0] not in HOTSPOT_GROUPS]

    # Rule 1: hotspot beats non-hotspot
    candidates = hotspots if hotspots else non_hotspot

    # Rule 2+3: sort by GROUP_TIER ascending first, then VAF descending
    # within the same tier. This ensures SF3B1/SRSF2/U2AF1 hotspots (tier 1)
    # always beat ZRSR2_mut (tier 2) regardless of VAF.
    candidates.sort(key=lambda h: (GROUP_TIER.get(h[0], 99), -h[4]))
    winner = candidates[0]

    # Build co_mutation string from all non-winning hits
    losers = [h for h in hits if h is not winner]
    co_str = ';'.join(f"{h[1]}:{h[2]}:VAF={h[3]}" for h in losers)

    group, gene, pc, vaf, _ = winner
    return (group, gene, pc, vaf, co_str)


# ── Build final metadata ──────────────────────────────────────────────────────
print("Building final metadata...")

FIELDNAMES = [
    'sample_id', 'cohort', 'tier', 'disease', 'sample_group',
    'mutated_gene', 'protein_change', 'vaf', 'co_mutation', 'bam_path',
    'SF3B1_call', 'SF3B1_change', 'SF3B1_vaf',
    'SRSF2_call',  'SRSF2_change',  'SRSF2_vaf',
    'U2AF1_call',  'U2AF1_change',  'U2AF1_vaf',
    'ZRSR2_call',  'ZRSR2_change',  'ZRSR2_vaf',
]

rows  = []
counts       = defaultdict(int)
co_mut_count = 0

for sample_id, mrow in manifest.items():
    disease = mrow.get('disease', '')
    meta    = metadata.get(sample_id, {})

    resolved = None
    if meta and disease == 'MDS_AML':
        hits     = get_all_sf_hits(sample_id, meta)
        resolved = resolve_hits(hits)
        if len(hits) > 1:
            co_mut_count += 1

    if resolved:
        sample_group, mutated_gene, protein_change, vaf, co_mutation = resolved
    else:
        sample_group   = DISEASE_LABELS.get(disease, disease)
        mutated_gene   = '.'
        protein_change = '.'
        vaf            = '.'
        co_mutation    = '.'

    counts[sample_group] += 1

    rows.append({
        'sample_id':      sample_id,
        'cohort':         mrow.get('cohort', '.'),
        'tier':           mrow.get('tier', '.'),
        'disease':        disease,
        'sample_group':   sample_group,
        'mutated_gene':   mutated_gene,
        'protein_change': protein_change,
        'vaf':            vaf,
        'co_mutation':    co_mutation,
        'bam_path':       meta.get('bam_path', '.'),
        'SF3B1_call':     meta.get('SF3B1_call',   'no_data'),
        'SF3B1_change':   meta.get('SF3B1_change',  '.'),
        'SF3B1_vaf':      meta.get('SF3B1_vaf',     '.'),
        'SRSF2_call':     meta.get('SRSF2_call',   'no_data'),
        'SRSF2_change':   meta.get('SRSF2_change',  '.'),
        'SRSF2_vaf':      meta.get('SRSF2_vaf',     '.'),
        'U2AF1_call':     meta.get('U2AF1_call',   'no_data'),
        'U2AF1_change':   meta.get('U2AF1_change',  '.'),
        'U2AF1_vaf':      meta.get('U2AF1_vaf',     '.'),
        'ZRSR2_call':     meta.get('ZRSR2_call',   'no_data'),
        'ZRSR2_change':   meta.get('ZRSR2_change',  '.'),
        'ZRSR2_vaf':      meta.get('ZRSR2_vaf',     '.'),
    })

# ── Write output ──────────────────────────────────────────────────────────────
with open(args.outfile, 'w', newline='') as fh:
    writer = csv.DictWriter(fh, delimiter='\t', fieldnames=FIELDNAMES)
    writer.writeheader()
    writer.writerows(rows)

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\nSample group breakdown:")
for group, count in sorted(counts.items()):
    print(f"  {group:<25} {count}")
print(f"  {'TOTAL':<25} {len(rows)}")
print(f"\n  Samples with co-mutations : {co_mut_count}")
print(f"  Output written to         : {args.outfile}")