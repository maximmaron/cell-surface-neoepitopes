#!/usr/bin/env Rscript
# =============================================================================
# 05_drimseq_dtu.R
#
# Differential Transcript Usage (DTU) using DRIMSeq + stageR.
# Runs as a SLURM array (1 job per SF group comparison).
#
# Usage: Rscript 05_drimseq_dtu.R <array_task_id>
#   1 = SF3B1_hotspot vs WT
#   2 = SRSF2_hotspot vs WT
#   3 = U2AF1_hotspot vs WT
#   4 = ZRSR2_mut vs WT
#
# Outputs per comparison (all under results/drimseq_dtu/):
#   {SF}_vs_WT_DTU.tsv                  - full results table
#   {SF}_vs_WT_DTU_sig.tsv              - significant hits only (stageR padj<0.05)
#   plots/{SF}_vs_WT/
#     01_sample_counts_density.pdf      - count distribution QC per sample
#     02_precision_dispersion.pdf       - dispersion vs expression
#     03_pvalue_histograms.pdf          - gene + transcript p-value distributions
#     04_top_dtu_genes.pdf              - isoform proportion boxplots, top 20 genes
#     05_sig_tx_heatmap.pdf             - row-scaled heatmap of significant transcripts
#     06_volcano_tx.pdf                 - proportion log2FC vs -log10(stageR padj)
# =============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DRIMSeq)
  library(stageR)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(BiocParallel)
})

# ── Config ────────────────────────────────────────────────────────────────────
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
BPPARAM <- MulticoreParam(N_CORES)
message(sprintf("Using %d cores", N_CORES))

# 8 comparisons: 4 SF groups x 2 control types
# Separating heme and tissue controls avoids confounding from pooling
# biologically heterogeneous tissues as a single reference
COMPARISONS <- data.frame(
  sf_group = rep(c("SF3B1_hotspot", "SRSF2_hotspot", "U2AF1_hotspot", "ZRSR2_mut"), 2),
  wt_group = c(rep("normal_heme", 4), rep("normal_tissue", 4)),
  stringsAsFactors = FALSE
)

PIPELINE_DIR <- "/data1/abdelwao/maxim/splicing_pipeline/kallisto"
T2G          <- file.path(PIPELINE_DIR, "reference/t2g.tsv")
MANIFEST     <- file.path(PIPELINE_DIR, "reference/sample_manifest.tsv")
QUANT_DIR    <- file.path(PIPELINE_DIR, "kallisto_quant")
OUT_DIR      <- file.path(PIPELINE_DIR, "results/drimseq_dtu")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PIPELINE_DIR, "logs"), showWarnings = FALSE)

# ── Parse array task ID ───────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
if (is.na(task_id) || task_id < 1 || task_id > nrow(COMPARISONS)) {
  stop(sprintf("Invalid task_id: %s. Must be 1-%d.", args[1], nrow(COMPARISONS)))
}
sf_group   <- COMPARISONS$sf_group[task_id]
wt_group   <- COMPARISONS$wt_group[task_id]
wt_label   <- ifelse(wt_group == "normal_heme", "heme", "tissue")
comparison <- sprintf("%s_vs_%s", sf_group, wt_label)
PLOT_DIR   <- file.path(OUT_DIR, "plots", comparison)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
message(sprintf("Task %d: %s (control: %s)", task_id, comparison, wt_group))

# ── Load data ─────────────────────────────────────────────────────────────────
t2g <- read_tsv(T2G, show_col_types = FALSE) %>%
  dplyr::select(transcript_id, gene_id, gene_symbol) %>%
  as.data.frame()

manifest <- read_tsv(MANIFEST, show_col_types = FALSE) %>%
  filter(sf_group %in% c(!!sf_group, !!wt_group)) %>%
  mutate(
    files = file.path(QUANT_DIR, sample_id, "abundance.h5"),
    group = factor(
      ifelse(sf_group == !!sf_group, "case", "control"),
      levels = c("control", "case")
    )
  )

n_case    <- sum(manifest$group == "case")
n_control <- sum(manifest$group == "control")
n_min     <- min(n_case, n_control)
message(sprintf("Samples: case=%d  control=%d", n_case, n_control))
if (n_case < 3 || n_control < 3) stop("Too few samples.")

txi <- tximport(
  files               = setNames(manifest$files, manifest$sample_id),
  type                = "kallisto",
  txOut               = TRUE,
  countsFromAbundance = "scaledTPM"
)
counts_mat <- round(txi$counts)
message(sprintf("Transcripts loaded: %d", nrow(counts_mat)))

# ── QC Plot 1: Count distributions ───────────────────────────────────────────
message("QC Plot 1: count distributions...")
pdf(file.path(PLOT_DIR, "01_sample_counts_density.pdf"), width = 14, height = 10)

# Total counts per sample bar plot
total_df <- data.frame(
  sample_id = colnames(counts_mat),
  total     = colSums(counts_mat)
) %>% left_join(manifest %>% dplyr::select(sample_id, group), by = "sample_id")

p1 <- ggplot(total_df, aes(x = reorder(sample_id, total), y = total / 1e6,
                            fill = group)) +
  geom_col() +
  scale_fill_manual(values = c(control = "#4393C3", case = "#D6604D")) +
  labs(title = paste0(comparison, " - Total scaled counts per sample"),
       x = "Sample", y = "Total counts (millions)", fill = "Group") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(p1)

# Log count density by group
count_long <- log1p(counts_mat) %>%
  as.data.frame() %>%
  rownames_to_column("transcript_id") %>%
  pivot_longer(-transcript_id, names_to = "sample_id", values_to = "log_count") %>%
  left_join(manifest %>% dplyr::select(sample_id, group), by = "sample_id")

p2 <- ggplot(count_long %>% filter(log_count > 0),
             aes(x = log_count, color = group, group = sample_id)) +
  geom_density(linewidth = 0.3, alpha = 0.4) +
  scale_color_manual(values = c(control = "#4393C3", case = "#D6604D")) +
  labs(title = "Log1p count density per sample",
       x = "log1p(scaled counts)", y = "Density", color = "Group") +
  theme_bw(base_size = 11)
print(p2)

dev.off()
message("  -> 01_sample_counts_density.pdf")

# ── Build DRIMSeq object and filter ──────────────────────────────────────────
counts_df <- as.data.frame(counts_mat) %>%
  rownames_to_column("feature_id") %>%
  left_join(t2g %>% dplyr::select(transcript_id, gene_id),
            by = c("feature_id" = "transcript_id")) %>%
  filter(!is.na(gene_id)) %>%
  dplyr::select(gene_id, feature_id, everything())

samp_df <- manifest %>%
  dplyr::select(sample_id, group) %>%
  as.data.frame()

d <- dmDSdata(counts = counts_df, samples = samp_df)

n_pre <- nrow(counts(d))
d <- dmFilter(d,
  min_samps_gene_expr    = floor(n_min * 0.5),
  min_samps_feature_expr = floor(n_min * 0.25),
  min_gene_expr          = 10,
  min_feature_expr       = 3
)
n_genes <- length(unique(counts(d)$gene_id))
n_tx    <- nrow(counts(d))
message(sprintf("After filter: %d genes, %d transcripts (removed %d)",
                n_genes, n_tx, n_pre - n_tx))

# ── Fit model ─────────────────────────────────────────────────────────────────
message("Estimating precision...")
set.seed(42)

# Build design matrix explicitly - required by dmPrecision and dmFit
# The full model includes group; the null model has only intercept for dmTest
design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))
design_null <- model.matrix(~ 1,     data = DRIMSeq::samples(d))

d <- dmPrecision(d, design = design_full, BPPARAM = BPPARAM)

# ── QC Plot 2: Dispersion vs expression (built-in DRIMSeq plot) ──────────────
message("QC Plot 2: dispersion...")
pdf(file.path(PLOT_DIR, "02_precision_dispersion.pdf"), width = 8, height = 5)
plotPrecision(d)
dev.off()
message("  -> 02_precision_dispersion.pdf")

message("Fitting and testing...")
d <- dmFit(d,  design = design_full, BPPARAM = BPPARAM)
d <- dmTest(d, design = design_null, BPPARAM = BPPARAM)

# ── Extract results ───────────────────────────────────────────────────────────
res_gene <- results(d, level = "gene") %>%
  dplyr::rename(gene_pvalue = pvalue, gene_adj_pvalue = adj_pvalue)

res_tx <- results(d, level = "feature") %>%
  dplyr::rename(tx_pvalue = pvalue, tx_adj_pvalue = adj_pvalue)

# ── QC Plot 3: P-value histograms (built-in DRIMSeq plots) ───────────────────
message("QC Plot 3: p-value histograms...")
pdf(file.path(PLOT_DIR, "03_pvalue_histograms.pdf"), width = 10, height = 8)
plotPValues(d, level = "gene")
plotPValues(d, level = "feature")
dev.off()
message("  -> 03_pvalue_histograms.pdf")

# ── stageR correction ─────────────────────────────────────────────────────────
message("Applying stageR correction...")

pconf <- res_tx %>%
  filter(!is.na(tx_pvalue)) %>%
  dplyr::select(feature_id, gene_id, tx_pvalue)

# Only keep genes in pscreen that have at least one transcript in pconf
# Mismatch occurs when gene-level results include genes whose transcripts
# were all filtered out at the transcript level
genes_with_tx <- unique(pconf$gene_id)
pscreen <- res_gene %>%
  filter(gene_id %in% genes_with_tx, !is.na(gene_pvalue)) %>%
  { setNames(.$gene_pvalue, .$gene_id) }

message(sprintf("  Genes in pscreen: %d  |  Transcripts in pconf: %d",
                length(pscreen), nrow(pconf)))

pconf_mat <- matrix(pconf$tx_pvalue, ncol = 1,
                    dimnames = list(pconf$feature_id, NULL))

stageRObj <- stageRTx(
  pScreen         = pscreen,
  pConfirmation   = pconf_mat,
  pScreenAdjusted = FALSE,
  tx2gene         = data.frame(tx = pconf$feature_id, gene = pconf$gene_id)
)
stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)
padj_stage <- getAdjustedPValues(stageRObj, order = FALSE,
                                 onlySignificantGenes = FALSE)

# ── Assemble results ──────────────────────────────────────────────────────────
result <- res_tx %>%
  left_join(res_gene %>% dplyr::select(gene_id, gene_pvalue, gene_adj_pvalue),
            by = "gene_id") %>%
  left_join(
    as.data.frame(padj_stage) %>%
      rownames_to_column("feature_id") %>%
      dplyr::rename(gene_padj_stageR = gene, tx_padj_stageR = transcript),
    by = "feature_id"
  ) %>%
  left_join(t2g, by = c("feature_id" = "transcript_id", "gene_id")) %>%
  dplyr::rename(transcript_id = feature_id) %>%
  mutate(
    comparison = comparison,
    is_novel   = grepl("^MSTRG", transcript_id)
  ) %>%
  arrange(gene_padj_stageR, tx_padj_stageR)

sig_genes <- result %>%
  filter(!is.na(gene_padj_stageR), gene_padj_stageR < 0.05) %>%
  pull(gene_id) %>% unique() %>% length()
sig_tx_n <- result %>%
  filter(!is.na(tx_padj_stageR), tx_padj_stageR < 0.05) %>% nrow()
message(sprintf("DTU genes (stageR padj<0.05):       %d", sig_genes))
message(sprintf("DTU transcripts (stageR padj<0.05): %d", sig_tx_n))

# ── QC Plot 4: Proportion plots for top 20 DTU genes (built-in) ──────────────
message("QC Plot 4: top DTU gene proportion plots...")
top_genes <- result %>%
  filter(!is.na(gene_padj_stageR)) %>%
  arrange(gene_padj_stageR) %>%
  pull(gene_id) %>%
  unique() %>%
  head(20)

pdf(file.path(PLOT_DIR, "04_top_dtu_genes.pdf"), width = 12, height = 5)
for (gid in top_genes) {
  tryCatch({
    print(plotProportions(d, gene_id = gid, group_variable = "group"))
  }, error = function(e) NULL)
}
dev.off()
message("  -> 04_top_dtu_genes.pdf")

# ── QC Plot 5: Heatmap ────────────────────────────────────────────────────────
message("QC Plot 5: heatmap...")
sig_tx_ids <- result %>%
  filter(!is.na(tx_padj_stageR), tx_padj_stageR < 0.05) %>%
  pull(transcript_id)

if (length(sig_tx_ids) >= 2) {
  tpm_sig <- log1p(txi$abundance[sig_tx_ids, , drop = FALSE])
  if (nrow(tpm_sig) > 200) {
    top200 <- result %>%
      filter(transcript_id %in% sig_tx_ids) %>%
      arrange(tx_padj_stageR) %>%
      slice_head(n = 200) %>%
      pull(transcript_id)
    tpm_sig <- tpm_sig[top200, ]
  }

  ann_col <- manifest %>%
    dplyr::select(sample_id, group) %>%
    column_to_rownames("sample_id") %>%
    as.data.frame()
  ann_colors <- list(group = c(control = "#4393C3", case = "#D6604D"))

  row_labels <- result %>%
    filter(transcript_id %in% rownames(tpm_sig)) %>%
    dplyr::select(transcript_id, gene_symbol) %>%
    distinct() %>%
    mutate(label = ifelse(is.na(gene_symbol), transcript_id,
                          sprintf("%s|%s", gene_symbol, transcript_id)))
  rlv <- setNames(row_labels$label, row_labels$transcript_id)[rownames(tpm_sig)]

  pdf(file.path(PLOT_DIR, "05_sig_tx_heatmap.pdf"),
      width = 14, height = max(8, nrow(tpm_sig) * 0.15 + 4))
  pheatmap(tpm_sig,
    annotation_col    = ann_col,
    annotation_colors = ann_colors,
    labels_row        = rlv,
    show_colnames     = FALSE,
    scale             = "row",
    clustering_method = "ward.D2",
    color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    fontsize_row      = 6,
    main              = sprintf("%s - Significant DTU transcripts (n=%d)",
                                comparison, length(sig_tx_ids))
  )
  dev.off()
  message(sprintf("  -> 05_sig_tx_heatmap.pdf  (%d transcripts shown)", nrow(tpm_sig)))
} else {
  message("  Skipping heatmap - fewer than 2 significant transcripts")
}

# ── QC Plot 6: Volcano ────────────────────────────────────────────────────────
message("QC Plot 6: volcano...")

# Compute mean proportion per group per transcript from DRIMSeq counts
prop_summary <- counts(d) %>%
  pivot_longer(cols = -c(gene_id, feature_id),
               names_to = "sample_id", values_to = "count") %>%
  left_join(manifest %>% dplyr::select(sample_id, group), by = "sample_id") %>%
  group_by(gene_id, feature_id, sample_id) %>%
  mutate(gene_total = sum(count)) %>%
  ungroup() %>%
  mutate(prop = count / pmax(gene_total, 1)) %>%
  group_by(feature_id, group) %>%
  summarise(mean_prop = mean(prop), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_prop) %>%
  mutate(prop_lfc = log2((case + 0.01) / (control + 0.01)))

volcano_df <- result %>%
  left_join(prop_summary, by = c("transcript_id" = "feature_id")) %>%
  filter(!is.na(prop_lfc), !is.na(tx_padj_stageR)) %>%
  mutate(
    neg_log10_p = -log10(pmax(tx_padj_stageR, 1e-300)),
    sig_label   = factor(case_when(
      tx_padj_stageR < 0.05 & prop_lfc >  0.5 ~ "Up in mutant",
      tx_padj_stageR < 0.05 & prop_lfc < -0.5 ~ "Down in mutant",
      TRUE ~ "NS"
    ), levels = c("Up in mutant", "Down in mutant", "NS")),
    label = ifelse(tx_padj_stageR < 0.01 & abs(prop_lfc) > 1 & !is.na(gene_symbol),
                   gene_symbol, NA_character_)
  )

pdf(file.path(PLOT_DIR, "06_volcano_tx.pdf"), width = 10, height = 8)
p6 <- ggplot(volcano_df, aes(x = prop_lfc, y = neg_log10_p, color = sig_label)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  scale_color_manual(
    values = c(
      "Up in mutant"   = "#D6604D",
      "Down in mutant" = "#4393C3",
      "NS"             = "#CCCCCC"
    ),
    drop = FALSE
  ) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 25,
                  segment.color = "grey60", na.rm = TRUE) +
  labs(
    title    = sprintf("%s - DTU transcript volcano", comparison),
    subtitle = sprintf("Up: %d  Down: %d  (stageR padj<0.05, |prop_lfc|>0.5)",
                       sum(volcano_df$sig_label == "Up in mutant"),
                       sum(volcano_df$sig_label == "Down in mutant")),
    x     = "log2(proportion: case / control)",
    y     = "-log10(stageR transcript padj)",
    color = NULL
  ) +
  theme_bw(base_size = 11)
print(p6)
dev.off()
message("  -> 06_volcano_tx.pdf")

# ── Save tables ───────────────────────────────────────────────────────────────
out_full <- file.path(OUT_DIR, sprintf("%s_DTU.tsv", comparison))
out_sig  <- file.path(OUT_DIR, sprintf("%s_DTU_sig.tsv", comparison))
write_tsv(result, out_full)
write_tsv(result %>% filter(!is.na(tx_padj_stageR), tx_padj_stageR < 0.05), out_sig)

message(sprintf("\nSaved: %s  (%d rows)", out_full, nrow(result)))
message(sprintf("Saved: %s  (%d significant)", out_sig, sig_tx_n))
message("DRIMSeq DTU complete.")