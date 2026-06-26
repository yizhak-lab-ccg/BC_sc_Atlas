# ============================================================
# TCGA-BRCA MYC CNV analysis (LOW vs HIGH)
# - Downloads gene-level CNV (GISTIC2) from UCSC Xena
# - Extracts MYC CNV
# - Matches by TCGA case (first 12 chars)
# - Collapses multiple aliquots per case (median)
# - Continuous test: Wilcoxon (HIGH vs LOW)
# - Discrete enrichment: Fisher (CNV>=+1 and CNV>=+2)
# - Saves tables + PDF plots (box, violin, heatmap, barplot)
#
# For publication: cite UCSC Xena + UCSCXenaTools
# ============================================================
library(UCSCXenaTools)
library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(grid)
library(scales)
library(tidyr)

# ----------------------------
# 0) USER SETTINGS
# ----------------------------
base_dir <- "path to folder/TCGA VALIDATION"

low_path  <- file.path(base_dir, "ER_ratio_LOW_samples.csv")#list of the ER low ratio
high_path <- file.path(base_dir, "ER_ratio_HIGH_samples.csv")# list of ER high ratio

# If your sample-id column isn't named "sample_id", change here:
sample_id_col <- "sample_id"

# Gene to analyze (easy to reuse for other genes)
gene_of_interest <- "MYC"

# Thresholds for discrete CNV enrichment
gain_threshold <- 1  # CNV >= 1 (GAIN/AMP)
amp_threshold  <- 2  # CNV >= 2 (AMP)

# Output subfolder for downloads
cnv_destdir <- file.path(base_dir, "TCGA_CNV_download")
dir.create(cnv_destdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) HELPERS
# ----------------------------
read_sample_list <- function(path, sample_col) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!sample_col %in% colnames(df)) {
    stop(paste("Missing column", sample_col, "in", path))
  }
  unique(na.omit(as.character(df[[sample_col]])))
}

pick_brca_gistic_by_genes_dataset <- function(prefer_thresholded = TRUE) {
  data("XenaData", package = "UCSCXenaTools")
  
  brca_ds <- XenaData %>%
    filter(XenaCohorts == "TCGA Breast Cancer (BRCA)") %>%
    as.data.frame()
  
  if (nrow(brca_ds) == 0) {
    stop("No rows found for 'TCGA Breast Cancer (BRCA)' in XenaData. Try: XenaDataUpdate() then restart R.")
  }
  
  cnv_candidates <- brca_ds %>%
    filter(grepl("Gistic2|gistic2", XenaDatasets, ignore.case = TRUE)) %>%
    filter(grepl("CopyNumber|copynumber|cnv", XenaDatasets, ignore.case = TRUE)) %>%
    filter(grepl("by_genes|genes", XenaDatasets, ignore.case = TRUE)) %>%
    arrange(XenaDatasets)
  
  if (nrow(cnv_candidates) == 0) {
    stop("Could not find a GISTIC2-by-genes CNV dataset in XenaData for BRCA.")
  }
  
  if (prefer_thresholded) {
    idx_best <- which(grepl("threshold", cnv_candidates$XenaDatasets, ignore.case = TRUE))
    if (length(idx_best) >= 1) return(cnv_candidates$XenaDatasets[idx_best[1]])
  }
  
  cnv_candidates$XenaDatasets[1]
}

download_prepare_cnv <- function(dataset_name, destdir) {
  cnv_file <- XenaGenerate() %>%
    XenaFilter(
      filterCohorts  = "TCGA Breast Cancer (BRCA)",
      filterDatasets = dataset_name
    ) %>%
    XenaQuery() %>%
    XenaDownload(destdir = destdir)
  
  XenaPrepare(cnv_file)
}

detect_gene_column <- function(cnv_mat) {
  gene_col <- intersect(
    c("Gene Symbol", "gene", "Gene", "Hugo_Symbol", "symbol"),
    colnames(cnv_mat)
  )
  if (length(gene_col) == 0) stop("No gene symbol column found. Check colnames(cnv_mat).")
  gene_col[1]
}

extract_gene_vector <- function(cnv_mat, gene_symbol, gene_col) {
  gene_row <- cnv_mat %>%
    filter(.data[[gene_col]] == gene_symbol)
  
  if (nrow(gene_row) == 0) {
    gene_row <- cnv_mat %>%
      filter(grepl(paste0("^", gene_symbol, "($|\\|)"), .data[[gene_col]]))
  }
  
  if (nrow(gene_row) == 0) {
    stop(paste0(gene_symbol, " not found in CNV table. Inspect cnv_mat[[gene_col]]."))
  }
  
  gene_row <- gene_row[1, ]
  sample_cols <- setdiff(colnames(cnv_mat), gene_col)
  
  vec <- as.numeric(gene_row[1, sample_cols] %>% unlist(use.names = FALSE))
  names(vec) <- sample_cols
  vec
}

build_case_table <- function(gene_vec, high_cases, low_cases, gene_name) {
  data.frame(
    sample = names(gene_vec),
    case   = substr(names(gene_vec), 1, 12),
    CNV    = as.numeric(gene_vec),
    stringsAsFactors = FALSE
  ) %>%
    group_by(case) %>%
    summarise(
      CNV = median(CNV, na.rm = TRUE),
      n_aliquots = n(),
      .groups = "drop"
    ) %>%
    mutate(
      group = case_when(
        case %in% high_cases ~ "HIGH",
        case %in% low_cases  ~ "LOW",
        TRUE ~ NA_character_
      ),
      gene = gene_name
    ) %>%
    filter(!is.na(group))
}

summarize_groups <- function(df) {
  df %>%
    group_by(group) %>%
    summarise(
      n      = n(),
      mean   = mean(CNV, na.rm = TRUE),
      sd     = sd(CNV, na.rm = TRUE),
      median = median(CNV, na.rm = TRUE),
      IQR    = IQR(CNV, na.rm = TRUE),
      .groups = "drop"
    )
}

wilcox_compare <- function(df) {
  wilcox.test(CNV ~ group, data = df, exact = FALSE)
}

fisher_enrichment <- function(df, threshold) {
  tab <- df %>%
    mutate(AMP = CNV >= threshold) %>%
    count(group, AMP) %>%
    tidyr::pivot_wider(names_from = AMP, values_from = n, values_fill = 0)
  
  if (!("FALSE" %in% colnames(tab))) tab$`FALSE` <- 0
  if (!("TRUE"  %in% colnames(tab))) tab$`TRUE`  <- 0
  
  m <- as.matrix(tab[, c("FALSE","TRUE")])
  rownames(m) <- tab$group
  fisher.test(m)
}

save_boxplot <- function(df, wt, outpath, title_txt, subtitle_txt = NULL) {
  pval_txt <- paste0("Wilcoxon P = ", signif(wt$p.value, 3))
  
  p <- ggplot(df, aes(x = group, y = CNV)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.7) +
    annotate(
      "text",
      x = 1.5,
      y = max(df$CNV, na.rm = TRUE) * 1.05,
      label = pval_txt,
      size = 5
    ) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "",
      y = "CNV (GISTIC2 call)"
    ) +
    theme_classic()
  
  ggsave(outpath, plot = p, width = 5, height = 5)
  p
}

save_violin <- function(df, wt, outpath, title_txt, subtitle_txt = NULL) {
  pval_txt <- paste0("Wilcoxon P = ", signif(wt$p.value, 3))
  
  p <- ggplot(df, aes(x = group, y = CNV, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.7) +
    annotate(
      "text",
      x = 1.5,
      y = max(df$CNV, na.rm = TRUE) * 1.05,
      label = pval_txt,
      size = 5
    ) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "",
      y = "CNV (GISTIC2 call)"
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  ggsave(outpath, plot = p, width = 5, height = 5)
  p
}

save_heatmap <- function(df, outpath, main_txt, low_first = TRUE) {
  heat_df <- df %>%
    filter(!is.na(group), !is.na(CNV), !is.na(case)) %>%
    mutate(group = factor(group, levels = if (low_first) c("LOW","HIGH") else c("HIGH","LOW"))) %>%
    arrange(group, desc(CNV))
  
  stopifnot(nrow(heat_df) > 0)
  
  case_ids <- make.unique(as.character(heat_df$case))
  
  heat_mat <- matrix(as.numeric(heat_df$CNV), nrow = 1)
  rownames(heat_mat) <- unique(heat_df$gene)[1]
  colnames(heat_mat) <- case_ids
  
  annotation_col <- data.frame(Group = heat_df$group)
  rownames(annotation_col) <- case_ids
  
  cnv_colors <- colorRampPalette(c("#2166ac", "#d1e5f0", "#f7f7f7", "#fddbc7", "#b2182b"))(100)
  breaks <- seq(-2, 2, length.out = 101)
  
  p <- pheatmap(
    heat_mat,
    color = cnv_colors,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    show_colnames = FALSE,
    border_color = NA,
    fontsize_row = 14,
    main = main_txt
  )
  
  pdf(outpath, width = 10, height = 3)
  grid.newpage()
  grid.draw(p$gtable)
  dev.off()
  
  p
}

save_amp_fraction_bar <- function(df, threshold, outpath, title_txt) {
  amp_df <- df %>%
    filter(!is.na(group), !is.na(CNV)) %>%
    mutate(
      group = factor(group, levels = c("LOW","HIGH")),
      AMP = CNV >= threshold
    ) %>%
    group_by(group) %>%
    summarise(
      n = n(),
      n_amp = sum(AMP, na.rm = TRUE),
      frac_amp = mean(AMP, na.rm = TRUE),
      .groups = "drop"
    )
  
  f <- fisher_enrichment(df, threshold)
  
  label_txt <- paste0(
    "CNV≥", threshold, ": Fisher P = ", signif(f$p.value, 3),
    " | OR = ", signif(unname(f$estimate), 3)
  )
  
  p <- ggplot(amp_df, aes(x = group, y = frac_amp, fill = group)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(n_amp, "/", n)), vjust = -0.6, size = 4) +
    annotate(
      "text",
      x = 1.5,
      y = max(amp_df$frac_amp) * 1.15,
      label = label_txt,
      size = 4.5
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, max(amp_df$frac_amp) * 1.25)
    ) +
    labs(title = title_txt, x = "", y = paste0("% CNV≥", threshold)) +
    theme_classic() +
    theme(legend.position = "none")
  
  ggsave(outpath, plot = p, width = 5.5, height = 4.5)
  
  list(plot = p, amp_rates = amp_df, fisher = f)
}

# ----------------------------
# 2) LOAD LOW/HIGH LISTS
# ----------------------------
low_samples  <- read_sample_list(low_path,  sample_id_col)
high_samples <- read_sample_list(high_path, sample_id_col)

low_cases  <- unique(substr(low_samples,  1, 12))
high_cases <- unique(substr(high_samples, 1, 12))

cat("LOW cases :", length(low_cases),  "\n")
cat("HIGH cases:", length(high_cases), "\n")

# ----------------------------
# 3) SELECT DATASET + DOWNLOAD CNV
# ----------------------------
cnv_dataset <- pick_brca_gistic_by_genes_dataset(prefer_thresholded = TRUE)
cat("\nUsing CNV dataset:\n", cnv_dataset, "\n")

cnv_mat <- download_prepare_cnv(cnv_dataset, cnv_destdir)
cat("\nCNV matrix dim:", paste(dim(cnv_mat), collapse = " x "), "\n")

gene_col <- detect_gene_column(cnv_mat)
cat("Using gene column:", gene_col, "\n")

# ----------------------------
# 4) EXTRACT GENE CNV + BUILD CASE TABLE
# ----------------------------
gene_vec <- extract_gene_vector(cnv_mat, gene_of_interest, gene_col)

case_df <- build_case_table(gene_vec, high_cases, low_cases, gene_of_interest) %>%
  mutate(group = factor(group, levels = c("LOW","HIGH")))  # consistent ordering

cat("\nMatched cases in CNV data:", nrow(case_df), "\n")
cat("HIGH matched:", sum(case_df$group == "HIGH"), "\n")
cat("LOW matched :", sum(case_df$group == "LOW"), "\n")

# ----------------------------
# 5) STATS
# ----------------------------
wt <- wilcox_compare(case_df)
print(wt)

summary_tbl <- summarize_groups(case_df)
cat("\nGroup summary:\n")
print(summary_tbl)

f_gain <- fisher_enrichment(case_df, gain_threshold)
f_amp  <- fisher_enrichment(case_df, amp_threshold)

cat("\nFisher enrichment:\n")
cat("CNV>=", gain_threshold, " p=", signif(f_gain$p.value, 3), " OR=", signif(unname(f_gain$estimate), 3), "\n")
cat("CNV>=", amp_threshold,  " p=", signif(f_amp$p.value, 3),  " OR=", signif(unname(f_amp$estimate), 3),  "\n")

# ----------------------------
# 6) SAVE TABLES
# ----------------------------
out_case_csv   <- file.path(base_dir, paste0(gene_of_interest, "_CNV_per_case_LOW_HIGH.csv"))
out_summ_csv   <- file.path(base_dir, paste0(gene_of_interest, "_CNV_summary_LOW_HIGH.csv"))

write.csv(case_df,    out_case_csv, row.names = FALSE)
write.csv(summary_tbl, out_summ_csv, row.names = FALSE)

# Also save discrete breakdown table (GAIN/AMP/OTHER)
breakdown_tbl <- case_df %>%
  mutate(
    status = case_when(
      CNV >= 2 ~ "AMP (+2)",
      CNV == 1 ~ "GAIN (+1)",
      TRUE ~ "NEUTRAL/LOSS (<=0)"
    )
  ) %>%
  count(group, status) %>%
  group_by(group) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

out_breakdown_csv <- file.path(base_dir, paste0(gene_of_interest, "_CNV_gain_amp_breakdown.csv"))
write.csv(breakdown_tbl, out_breakdown_csv, row.names = FALSE)

# ----------------------------
# 7) PLOTS (PDF)
# ----------------------------
title_core <- paste0(gene_of_interest, " CNV (TCGA-BRCA; GISTIC2 gene-level)")
subtitle_core <- "LOW vs HIGH ER+ signature"

# Boxplot
out_box <- file.path(base_dir, paste0(gene_of_interest, "_CNV_boxplot_HIGH_vs_LOW.pdf"))
p_box <- save_boxplot(case_df, wt, out_box, title_core, subtitle_core)

# Violin
out_violin <- file.path(base_dir, paste0(gene_of_interest, "_CNV_violin_HIGH_vs_LOW.pdf"))
p_violin <- save_violin(case_df, wt, out_violin, title_core, subtitle_core)

# Heatmap
out_hm <- file.path(base_dir, paste0(gene_of_interest, "_CNV_heatmap_LOW_vs_HIGH.pdf"))
p_hm <- save_heatmap(case_df, out_hm, paste0(gene_of_interest, " CNV (GISTIC2) — LOW vs HIGH"))

# Barplot: fraction CNV>=+1
out_bar_gain <- file.path(base_dir, paste0(gene_of_interest, "_CNV_gain_fraction_barplot.pdf"))
bar_gain <- save_amp_fraction_bar(
  case_df, gain_threshold, out_bar_gain,
  paste0("Fraction of ", gene_of_interest, " gain/amplification (CNV≥", gain_threshold, ")")
)

# Barplot: fraction CNV>=+2
out_bar_amp <- file.path(base_dir, paste0(gene_of_interest, "_CNV_amp_fraction_barplot.pdf"))
bar_amp <- save_amp_fraction_bar(
  case_df, amp_threshold, out_bar_amp,
  paste0("Fraction of ", gene_of_interest, " amplification (CNV≥", amp_threshold, ")")
)

# ----------------------------
# 8) RUN LOG (helpful for supplement)
# ----------------------------
cat("\nSaved outputs:\n")
cat(out_case_csv, "\n")
cat(out_summ_csv, "\n")
cat(out_breakdown_csv, "\n")
cat(out_box, "\n")
cat(out_violin, "\n")
cat(out_hm, "\n")
cat(out_bar_gain, "\n")
cat(out_bar_amp, "\n")

cat("\nSession info:\n")
print(sessionInfo())
