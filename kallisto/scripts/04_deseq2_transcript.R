#!/usr/bin/env Rscript
#SBATCH --job-name=deseq2_surface
#SBATCH --output=logs/04_deseq2_surface_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
# =============================================================================
# 04_deseq2_transcript_surface.R
#
# Transcript-level differential expression using kallisto abundances.
# Focused on surface-protein-encoding transcripts with strict biological
# filters to reduce false positives and focus on actionable candidates.
#
# Strategy to control false positives:
#
#   1. STATISTICAL: lfcThreshold = 1 changes the null hypothesis from
#      "LFC != 0" to "LFC <= 1", meaning padj reflects evidence that
#      the fold change EXCEEDS 2-fold — not merely that it is nonzero.
#      This is far more stringent than post-hoc LFC filtering.
#
#   2. MINIMUM EXPRESSION IN CASES: transcript must have median TPM >= 1
#      in the case group. Eliminates transcripts that are significant only
#      because they drop to zero in controls (not useful for targeting).
#
#   3. CONSISTENCY: transcript must be detectably expressed (TPM >= 0.5)
#      in >= 25% of case samples. Removes outlier-driven signals.
#
#   4. LOW EXPRESSION IN NORMALS: median TPM < 0.5 in BOTH normal_heme
#      and normal_tissue controls (evaluated across all normal samples,
#      not just those in the current contrast). This is the key filter for
#      tumour-selective surface targets — you want transcripts that are
#      absent or very low in healthy tissue.
#
#   5. BATCH CORRECTION: sequencing_batch (collapsed cohort) included in
#      design where estimable.
#
# Outputs per comparison:
#   results/deseq2_surface/
#     {SF}_vs_WT_surface_transcript_DE.tsv   — full results with filter flags
#     {SF}_vs_WT_volcano.pdf                 — volcano plot with top hits labelled
#     candidate_surface_transcripts.tsv      — priority candidates (all comparisons)
#     igv_regions.txt                        — IGV-ready locus strings
#     igv_batch.txt                          — IGV batch screenshot script
# =============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(BiocParallel)
})

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
message(sprintf("Using %d cores", N_CORES))

# ── Config ────────────────────────────────────────────────────────────────────
TEST_MODE    <- FALSE
TEST_N_CASE  <- 5
TEST_N_CTRL  <- 10

LFC_THRESHOLD      <- 1       # null hypothesis: |LFC| <= 1 (tests for >2-fold)
MIN_TPM_CASE       <- 1.0     # min median TPM in case samples
MIN_FRAC_EXPRESSED <- 0.25    # min fraction of cases with TPM >= 0.5
MAX_TPM_NORMAL     <- 0.5     # max median TPM allowed in normals
PADJ_THRESHOLD     <- 0.05
IGV_PADDING        <- 2000    # bp padding around transcript for IGV view
TOP_N_LABEL        <- 20      # top N transcripts labelled on volcano
TOP_N_IGV          <- 100     # top N candidates for IGV batch script

PIPELINE_DIR  <- "/data1/abdelwao/maxim/splicing_pipeline/kallisto"
GTF_MERGED    <- "/data1/abdelwao/maxim/splicing_pipeline/stringtie/merge/merged_final.gtf"
T2G           <- file.path(PIPELINE_DIR, "reference/t2g.tsv")
MANIFEST      <- file.path(PIPELINE_DIR, "reference/sample_manifest.tsv")
QUANT_DIR     <- file.path(PIPELINE_DIR, "kallisto_quant")
SURFACE_GENES <- file.path(PIPELINE_DIR, "reference/surface_genes.tsv")
OUT_DIR       <- file.path(PIPELINE_DIR,
                           if (TEST_MODE) "results/deseq2_surface_test"
                           else           "results/deseq2_surface")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PIPELINE_DIR, "logs"), showWarnings = FALSE)

# ── Cohort to sequencing batch ────────────────────────────────────────────────
cohort_to_batch <- function(cohort) {
  dplyr::case_when(
    cohort %in% c("tcga", "bodymap2") ~ "unstranded_50bp",
    cohort %in% c("hpa_2014")         ~ "unstranded_101bp",
    TRUE                              ~ "stranded_100bp"
  )
}

# ── Load reference tables ─────────────────────────────────────────────────────
t2g_df <- read_tsv(T2G, show_col_types = FALSE) %>%
  mutate(gene_id_base = sub("\\.\\d+$", "", gene_id))

surf_df <- read_tsv(SURFACE_GENES, show_col_types = FALSE)
names(surf_df) <- trimws(names(surf_df))
surf_ensg <- unique(sub("\\.\\d+$", "", surf_df$ensg_id))
message(sprintf("Surface gene list: %d unique ENSG IDs", length(surf_ensg)))

# ── Load GTF for genomic coordinates (used for IGV) ───────────────────────────
message("Parsing GTF for transcript coordinates...")
gtf_coords <- local({
  con  <- file(GTF_MERGED, "r")
  rows <- list()
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (grepl("^\t", line) || !grepl("\ttranscript\t", line)) next
    tid  <- regmatches(line, regexpr('(?<=transcript_id ")([^"]+)', line, perl=TRUE))
    if (length(tid) == 0) next
    f    <- strsplit(line, "\t")[[1]]
    rows[[length(rows)+1]] <- list(
      transcript_id = tid,
      chrom  = f[1], start = as.integer(f[4]),
      end    = as.integer(f[5]), strand = f[7],
      class_code = {
        m <- regmatches(line, regexpr('(?<=class_code ")([^"]+)', line, perl=TRUE))
        if (length(m)) m else NA_character_
      }
    )
  }
  close(con)
  bind_rows(rows)
})
message(sprintf("  GTF transcripts loaded: %d", nrow(gtf_coords)))

# ── Load manifest ─────────────────────────────────────────────────────────────
manifest <- read_tsv(MANIFEST, show_col_types = FALSE) %>%
  mutate(
    files            = file.path(QUANT_DIR, sample_id, "abundance.h5"),
    sequencing_batch = cohort_to_batch(cohort)
  )

missing <- manifest$files[!file.exists(manifest$files)]
if (length(missing) > 0)
  stop(sprintf("%d abundance.h5 files missing. First: %s",
               length(missing), missing[1]))
message(sprintf("Manifest loaded: %d samples", nrow(manifest)))

# ── Define groups ─────────────────────────────────────────────────────────────
SF_GROUPS <- c("SF3B1_hotspot", "SRSF2_hotspot", "U2AF1_hotspot", "ZRSR2_mut")
WT_GROUPS <- c("normal_heme", "normal_tissue")

# ── Subset manifest ───────────────────────────────────────────────────────────
if (TEST_MODE) {
  case_sub <- manifest %>%
    filter(sf_group %in% SF_GROUPS) %>%
    group_by(sf_group) %>% slice_head(n = TEST_N_CASE) %>% ungroup()
  ctrl_sub <- manifest %>%
    filter(sf_group %in% WT_GROUPS) %>% slice_head(n = TEST_N_CTRL)
  manifest <- bind_rows(case_sub, ctrl_sub)
  message(sprintf("TEST MODE: %d samples", nrow(manifest)))
} else {
  manifest <- manifest %>%
    filter(sf_group %in% c(SF_GROUPS, WT_GROUPS, "MDS_AML_SFWT"))
  message(sprintf("Full run: %d samples", nrow(manifest)))
}

# ── Import all quantifications once ──────────────────────────────────────────
message("Importing kallisto quantifications with tximport...")
txi <- tximport(
  files               = setNames(manifest$files, manifest$sample_id),
  type                = "kallisto",
  txOut               = TRUE,
  countsFromAbundance = "scaledTPM"
)
message(sprintf("Total transcripts imported: %d", nrow(txi$counts)))

# ── Filter to surface transcripts ────────────────────────────────────────────
surf_tx <- t2g_df %>%
  filter(gene_id_base %in% surf_ensg) %>%
  pull(transcript_id) %>% unique()
surf_tx_present <- intersect(surf_tx, rownames(txi$counts))
message(sprintf("Surface transcripts in quantification: %d / %d",
                length(surf_tx_present), length(surf_tx)))
if (length(surf_tx_present) == 0)
  stop("No surface transcripts found — check ENSG ID format.")

txi_surf <- list(
  counts              = txi$counts[surf_tx_present, ],
  length              = txi$length[surf_tx_present, ],
  abundance           = txi$abundance[surf_tx_present, ],
  countsFromAbundance = txi$countsFromAbundance
)
message(sprintf("Proceeding with %d surface transcripts", nrow(txi_surf$counts)))

# ── Pre-compute normal TPM summaries ─────────────────────────────────────────
heme_samples   <- manifest %>% filter(sf_group == "normal_heme")   %>% pull(sample_id)
tissue_samples <- manifest %>% filter(sf_group == "normal_tissue") %>% pull(sample_id)

tpm_normal_heme   <- txi_surf$abundance[, heme_samples,   drop = FALSE]
tpm_normal_tissue <- txi_surf$abundance[, tissue_samples, drop = FALSE]

median_tpm_heme   <- apply(tpm_normal_heme,   1, median)
median_tpm_tissue <- apply(tpm_normal_tissue, 1, median)

message(sprintf("Normal expression: %d heme + %d tissue samples",
                length(heme_samples), length(tissue_samples)))

# ── Helper: is batch term estimable? ─────────────────────────────────────────
batch_is_estimable <- function(col_df) {
  tab <- table(col_df$sequencing_batch, col_df$group)
  any(apply(tab, 1, function(r) all(r > 0)))
}

# ── Helper: biological filters ───────────────────────────────────────────────
apply_biological_filters <- function(res_df, case_sid) {
  tpm_case       <- txi_surf$abundance[, case_sid, drop = FALSE]
  median_tpm_case <- apply(tpm_case, 1, median)
  frac_expressed  <- rowMeans(tpm_case >= 0.5)

  filter_df <- tibble(
    transcript_id     = rownames(tpm_case),
    median_tpm_case   = median_tpm_case,
    frac_expressed    = frac_expressed,
    median_tpm_heme   = median_tpm_heme  [rownames(tpm_case)],
    median_tpm_tissue = median_tpm_tissue[rownames(tpm_case)]
  )

  res_df %>%
    left_join(filter_df, by = "transcript_id") %>%
    mutate(
      pass_expr_case   = median_tpm_case   >= MIN_TPM_CASE,
      pass_consistency = frac_expressed    >= MIN_FRAC_EXPRESSED,
      pass_low_heme    = median_tpm_heme   <  MAX_TPM_NORMAL,
      pass_low_tissue  = median_tpm_tissue <  MAX_TPM_NORMAL,
      pass_all_filters = pass_expr_case & pass_consistency &
                         pass_low_heme  & pass_low_tissue,
      is_novel         = grepl("^MSTRG", transcript_id)
    )
}

# ── Helper: volcano plot ──────────────────────────────────────────────────────
make_volcano <- function(res_df, contrast_name, out_pdf) {

  # Categorise each transcript
  df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10_p = -log10(pmax(padj, 1e-300)),
      category    = case_when(
        padj < PADJ_THRESHOLD & log2FoldChange > 0 & mutation_specific & is_novel
                            ~ "Mutation-specific: novel",
        padj < PADJ_THRESHOLD & log2FoldChange > 0 & mutation_specific
                            ~ "Mutation-specific: annotated",
        padj < PADJ_THRESHOLD & log2FoldChange > 0 & pass_all_filters & in_disease_blacklist
                            ~ "Disease effect (not specific)",
        padj < PADJ_THRESHOLD & log2FoldChange > 0 & !pass_all_filters
                            ~ "Sig (filtered out)",
        padj < PADJ_THRESHOLD & log2FoldChange < 0
                            ~ "Sig (downregulated)",
        TRUE                ~ "NS"
      ),
      category = factor(category, levels = c(
        "Mutation-specific: novel", "Mutation-specific: annotated",
        "Disease effect (not specific)",
        "Sig (filtered out)", "Sig (downregulated)", "NS"
      ))
    )

  # Top candidates to label: priority hits sorted by padj
  top_label <- df %>%
    filter(category %in% c("Mutation-specific: novel", "Mutation-specific: annotated")) %>%
    arrange(padj) %>%
    slice_head(n = TOP_N_LABEL) %>%
    mutate(label = ifelse(!is.na(gene_symbol), gene_symbol,
                          sub("MSTRG\\.(\\d+)\\..*", "MSTRG.\\1", transcript_id)))

  colour_map <- c(
    "Mutation-specific: novel"      = "#B2182B",
    "Mutation-specific: annotated"  = "#D6604D",
    "Disease effect (not specific)" = "#F4A582",
    "Sig (filtered out)"            = "#92C5DE",
    "Sig (downregulated)"           = "#4393C3",
    "NS"                            = "#CCCCCC"
  )
  size_map <- c(
    "Mutation-specific: novel"      = 2.5,
    "Mutation-specific: annotated"  = 2.0,
    "Disease effect (not specific)" = 1.2,
    "Sig (filtered out)"            = 0.8,
    "Sig (downregulated)"           = 0.8,
    "NS"                            = 0.4
  )
  alpha_map <- c(
    "Mutation-specific: novel"      = 1.0,
    "Mutation-specific: annotated"  = 0.9,
    "Disease effect (not specific)" = 0.5,
    "Sig (filtered out)"            = 0.3,
    "Sig (downregulated)"           = 0.3,
    "NS"                            = 0.2
  )

  n_priority <- sum(df$category %in%
                    c("Mutation-specific: novel", "Mutation-specific: annotated"))
  n_novel    <- sum(df$category == "Mutation-specific: novel")
  n_disease  <- sum(df$category == "Disease effect (not specific)")

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_p,
                      colour = category, size = category, alpha = category)) +
    geom_point() +
    scale_colour_manual(values = colour_map, drop = FALSE) +
    scale_size_manual(  values = size_map,   drop = FALSE) +
    scale_alpha_manual( values = alpha_map,  drop = FALSE) +
    geom_hline(yintercept = -log10(PADJ_THRESHOLD),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0,
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_text_repel(
      data          = top_label,
      aes(label     = label),
      size          = 2.8,
      colour        = "black",
      fontface      = ifelse(top_label$is_novel, "bold.italic", "plain"),
      max.overlaps  = 30,
      segment.color = "grey60",
      segment.size  = 0.3,
      box.padding   = 0.4
    ) +
    labs(
      title    = contrast_name,
      subtitle = sprintf(
        "Mutation-specific: %d  |  Novel: %d  |  Disease effect: %d  |  Filtered: %d",
        n_priority, n_novel, n_disease,
        sum(df$category == "Sig (filtered out)")
      ),
      x        = paste0("log2 Fold Change (apeglm shrunken, null: |LFC| <= ",
                         LFC_THRESHOLD, ")"),
      y        = paste0("-log10(padj)"),
      colour   = NULL, size = NULL, alpha = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 9),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank()
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
      size   = "none",
      alpha  = "none"
    )

  ggsave(out_pdf, p, width = 10, height = 8, device = "pdf")
  message(sprintf("    Volcano -> %s", out_pdf))
  invisible(p)
}

# ── Helper: run DESeq2 for one contrast ──────────────────────────────────────
run_deseq2 <- function(sub_manifest, contrast_name) {

  n_case <- sum(sub_manifest$group == "case")
  n_ctrl <- sum(sub_manifest$group == "control")
  message(sprintf("\n  Running DESeq2: %s  (n_case=%d, n_ctrl=%d)",
                  contrast_name, n_case, n_ctrl))

  sid      <- sub_manifest$sample_id
  case_sid <- sub_manifest %>% filter(group == "case") %>% pull(sample_id)

  sub_txi <- list(
    counts              = txi_surf$counts[, sid],
    length              = txi_surf$length[, sid],
    abundance           = txi_surf$abundance[, sid],
    countsFromAbundance = txi_surf$countsFromAbundance
  )

  col_df <- sub_manifest %>%
    dplyr::select(sample_id, group, cohort, sequencing_batch) %>%
    mutate(
      group            = relevel(factor(group),            ref = "control"),
      sequencing_batch = relevel(factor(sequencing_batch), ref = "stranded_100bp")
    ) %>%
    as.data.frame()
  rownames(col_df) <- col_df$sample_id

  message("    Batch x group:")
  tab <- table(col_df$sequencing_batch, col_df$group)
  message(paste0("    ", capture.output(print(tab)), collapse = "\n"))

  design_form <- if (batch_is_estimable(col_df)) {
    message("    Design: ~ sequencing_batch + group")
    ~ sequencing_batch + group
  } else {
    message("    WARNING: batch confounded — using ~ group only")
    ~ group
  }

  dds <- DESeqDataSetFromTximport(sub_txi, col_df, design = design_form)
  n_min <- max(3, floor(ncol(dds) * 0.1))
  keep  <- rowSums(counts(dds) >= 10) >= n_min
  dds   <- dds[keep, ]
  message(sprintf("    Transcripts after count filter: %d", nrow(dds)))

  dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(N_CORES))

  res <- lfcShrink(dds,
    coef         = "group_case_vs_control",
    type         = "apeglm",
    lfcThreshold = LFC_THRESHOLD,
    parallel     = TRUE,
    BPPARAM      = MulticoreParam(N_CORES)
  )

  # When lfcThreshold > 0, apeglm returns svalue instead of padj.
  # Rename handled in res_df annotation block below.
  res_cols <- names(as.data.frame(res))
  message(sprintf("    lfcShrink columns: %s", paste(res_cols, collapse = ", ")))

  # Annotate
  res_df <- as.data.frame(res) %>%
    rownames_to_column("transcript_id") %>%
    # Rename svalue to padj if lfcThreshold was used (apeglm returns svalue not padj)
    { if ("svalue" %in% names(.) && !"padj" %in% names(.))
        dplyr::rename(., padj = svalue) else . } %>%
    left_join(t2g_df %>% dplyr::select(transcript_id, gene_id, gene_id_base),
              by = "transcript_id") %>%
    left_join(
      surf_df %>%
        dplyr::select(gene_symbol, ensg_id) %>%
        mutate(ensg_id_base = sub("\\.\\d+$", "", ensg_id)) %>%
        distinct(),
      by = c("gene_id_base" = "ensg_id_base")
    ) %>%
    # Add genomic coordinates from GTF
    left_join(
      gtf_coords %>% dplyr::select(transcript_id, chrom, start, end,
                                    strand, class_code),
      by = "transcript_id"
    ) %>%
    mutate(
      comparison  = contrast_name,
      igv_locus   = ifelse(
        !is.na(chrom),
        sprintf("%s:%d-%d", chrom,
                pmax(1L, start - IGV_PADDING),
                end + IGV_PADDING),
        NA_character_
      )
    ) %>%
    arrange(padj)

  # Apply biological filters
  res_df <- apply_biological_filters(res_df, case_sid)

  # Report filter cascade
  stat_sig <- sum(res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > 0,
                  na.rm = TRUE)
  bio_pass <- sum(res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > 0 &
                  res_df$pass_all_filters, na.rm = TRUE)

  message(sprintf("    Total tested:                    %d", nrow(res_df)))
  message(sprintf("    Statistically significant (UP):  %d  (padj/svalue<%.2f, LFC>%g)",
                  stat_sig, PADJ_THRESHOLD, LFC_THRESHOLD))
  message(sprintf("    After biological filters:        %d", bio_pass))
  message(sprintf("      pass expr (TPM>=%.1f):          %d",
                  MIN_TPM_CASE,
                  sum(res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > 0 &
                      res_df$pass_expr_case, na.rm = TRUE)))
  message(sprintf("      pass consistency (>=%.0f%%):     %d",
                  MIN_FRAC_EXPRESSED * 100,
                  sum(res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > 0 &
                      res_df$pass_consistency, na.rm = TRUE)))
  message(sprintf("      pass low heme (<%.1f):          %d",
                  MAX_TPM_NORMAL,
                  sum(res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > 0 &
                      res_df$pass_low_heme, na.rm = TRUE)))
  message(sprintf("      pass low tissue (<%.1f):        %d",
                  MAX_TPM_NORMAL,
                  sum(res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > 0 &
                      res_df$pass_low_tissue, na.rm = TRUE)))

  # Volcano plot
  volcano_path <- file.path(OUT_DIR,
                             sprintf("%s_volcano.pdf", contrast_name))
  tryCatch(
    make_volcano(res_df, contrast_name, volcano_path),
    error = function(e) message(sprintf("    Volcano plot failed: %s", e$message))
  )

  return(res_df)
}

# ── Run SFWT vs WT to identify disease-effect transcripts (blacklist) ─────────
# Transcripts upregulated in SFWT MDS/AML vs normal are driven by disease,
# not by the SF mutation. These must be removed from all mutant comparisons.
message("\n=== MDS_AML_SFWT vs WT (disease blacklist) ===")

sfwt_rows <- manifest %>% filter(sf_group == "MDS_AML_SFWT")
ctrl_rows  <- manifest %>% filter(sf_group %in% WT_GROUPS)

disease_blacklist <- character(0)

if (nrow(sfwt_rows) >= 3 && nrow(ctrl_rows) >= 3) {
  sfwt_manifest <- bind_rows(
    sfwt_rows %>% mutate(group = "case"),
    ctrl_rows %>% mutate(group = "control")
  )
  tryCatch({
    sfwt_res <- run_deseq2(sfwt_manifest, "MDS_AML_SFWT_vs_WT")
    write_tsv(sfwt_res,
              file.path(OUT_DIR, "MDS_AML_SFWT_vs_WT_surface_transcript_DE.tsv"))

    # Blacklist: transcripts significantly upregulated in SFWT
    disease_blacklist <- sfwt_res %>%
      filter(padj < PADJ_THRESHOLD, log2FoldChange > 0) %>%
      pull(transcript_id)

    message(sprintf("  Disease-effect blacklist: %d transcripts",
                    length(disease_blacklist)))
  }, error = function(e) {
    message(sprintf("  ERROR in SFWT blacklist: %s", e$message))
  })
} else {
  message("  Skipping — insufficient SFWT samples in manifest")
  message("  Add MDS_AML_SFWT to WT_GROUPS filter to include these samples")
}

# ── Run all comparisons ───────────────────────────────────────────────────────
all_results <- list()

for (sf in SF_GROUPS) {
  message(sprintf("\n=== %s vs pooled WT ===", sf))

  case_rows <- manifest %>% filter(sf_group == sf)
  ctrl_rows <- manifest %>% filter(sf_group %in% WT_GROUPS)

  if (nrow(case_rows) < 3 || nrow(ctrl_rows) < 3) {
    message("  Skipping — too few samples.")
    next
  }

  sub_manifest <- bind_rows(
    case_rows %>% mutate(group = "case"),
    ctrl_rows %>% mutate(group = "control")
  )

  contrast_name <- sprintf("%s_vs_WT", sf)

  tryCatch({
    res_df   <- run_deseq2(sub_manifest, contrast_name)

    # Apply disease-effect blacklist
    # Transcripts also upregulated in SFWT are disease markers, not mutation-specific
    res_df <- res_df %>%
      mutate(
        in_disease_blacklist = transcript_id %in% disease_blacklist,
        mutation_specific    = pass_all_filters & !in_disease_blacklist
      )

    n_blacklisted <- sum(res_df$padj < PADJ_THRESHOLD &
                         res_df$log2FoldChange > 0 &
                         res_df$pass_all_filters &
                         res_df$in_disease_blacklist, na.rm = TRUE)
    message(sprintf("    Removed by disease blacklist: %d", n_blacklisted))
    message(sprintf("    Mutation-specific candidates: %d",
                    sum(res_df$mutation_specific, na.rm = TRUE)))

    out_path <- file.path(OUT_DIR,
                          sprintf("%s_surface_transcript_DE.tsv", contrast_name))
    write_tsv(res_df, out_path)
    all_results[[contrast_name]] <- res_df
  }, error = function(e) {
    message(sprintf("  ERROR in %s: %s", contrast_name, e$message))
  })
}

# ── Combine results and write candidate files ─────────────────────────────────
message("\n", strrep("=", 60))
message("SUMMARY")
message(strrep("-", 60))

for (nm in names(all_results)) {
  df   <- all_results[[nm]]
  stat <- sum(df$padj < PADJ_THRESHOLD & df$log2FoldChange > 0, na.rm = TRUE)
  bio  <- sum(df$padj < PADJ_THRESHOLD & df$log2FoldChange > 0 &
              df$pass_all_filters, na.rm = TRUE)
  mut  <- sum(df$mutation_specific, na.rm = TRUE)
  nov  <- sum(df$mutation_specific & df$is_novel, na.rm = TRUE)
  message(sprintf("  %-40s  stat=%d  bio=%d  mutation_specific=%d  novel=%d",
                  nm, stat, bio, mut, nov))
}

if (length(all_results) > 0) {
  combined <- bind_rows(all_results)
  write_tsv(combined,
            file.path(OUT_DIR, "all_comparisons_surface_transcript_DE.tsv"))

  # Priority candidates: mutation-specific (pass bio filters + not in disease blacklist)
  candidates <- combined %>%
    filter(padj < PADJ_THRESHOLD, log2FoldChange > 0, mutation_specific) %>%
    dplyr::select(any_of(c(
      "transcript_id", "gene_symbol", "gene_id",
      "chrom", "start", "end", "strand", "class_code",
      "log2FoldChange", "lfcSE", "pvalue", "padj", "svalue",
      "median_tpm_case", "frac_expressed",
      "median_tpm_heme", "median_tpm_tissue",
      "is_novel", "in_disease_blacklist", "mutation_specific",
      "igv_locus", "comparison"
    ))) %>%
    arrange(comparison, padj)

  write_tsv(candidates,
            file.path(OUT_DIR, "candidate_surface_transcripts.tsv"))
  message(sprintf("\nPriority candidates: %d", nrow(candidates)))

  # ── IGV files ─────────────────────────────────────────────────────────────
  # Top candidates sorted by padj across all comparisons
  top_igv <- candidates %>%
    arrange(padj) %>%
    filter(!is.na(igv_locus)) %>%
    # Deduplicate by transcript (may appear in multiple comparisons)
    distinct(transcript_id, .keep_all = TRUE) %>%
    slice_head(n = TOP_N_IGV)

  # igv_regions.txt — paste into IGV search box
  writeLines(top_igv$igv_locus,
             file.path(OUT_DIR, "igv_regions.txt"))

  # igv_batch.txt — automated screenshot script
  igv_batch_lines <- c(
    "new",
    "genome hg38",
    "snapshotDirectory igv_snapshots",
    ""
  )
  for (i in seq_len(nrow(top_igv))) {
    row   <- top_igv[i, ]
    gsym  <- if (!is.na(row$gene_symbol)) row$gene_symbol else row$gene_id
    tid   <- row$transcript_id
    ccode <- if (!is.na(row$class_code)) row$class_code else "?"
    igv_batch_lines <- c(
      igv_batch_lines,
      sprintf("goto %s", row$igv_locus),
      sprintf("snapshot %s_%s_%s.png", gsym, tid, ccode),
      ""
    )
  }
  writeLines(igv_batch_lines, file.path(OUT_DIR, "igv_batch.txt"))

  message(sprintf("IGV regions:  %s/igv_regions.txt  (%d loci)",
                  OUT_DIR, nrow(top_igv)))
  message(sprintf("IGV batch:    %s/igv_batch.txt", OUT_DIR))
}

message(sprintf("\nDone (%s).", if (TEST_MODE) "test run" else "full run"))