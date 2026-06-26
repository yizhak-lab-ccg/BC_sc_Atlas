# ============================================================
# TCGA-BRCA ATAC-seq (UCSC Xena ATAC Hub):
# Download + HIGH vs LOW Differential Accessibility analysis
#
# Includes:
#   A) Discover datasets on atacseq.xenahubs.net
#   B) Download BRCA "All peaks" log2counts matrix
#   C) Download probeMap (peak coordinates) and attach to peaks
#   D) Match LOW/HIGH sample lists by TCGA case (first 12 chars)
#   E) Optional: collapse multiple aliquots per case (median)
#   F) Differential accessibility per peak (Wilcoxon + BH FDR)
#   G) Save tables + plots (volcano, heatmap, PCA)
#
# Outputs in base_dir:
#   - TCGA_BRCA_ATAC_all_peaks_Log2Counts.rds
#   - brca_peak.probeMap  (downloaded annotation)
#   - TCGA_BRCA_ATAC_peak_coordinates.csv
#   - ATAC_LOW_HIGH_column_mapping.csv
#   - ATAC_DA_peaks_High_vs_Low.csv
#   - ATAC_DA_peaks_High_vs_Low_with_coords.csv
#   - ATAC_DA_peaks_High_vs_Low_significant.csv
#   - ATAC_volcano_HIGH_vs_LOW.pdf
#   - ATAC_heatmap_topDA.pdf
#   - ATAC_PCA_LOW_vs_HIGH.pdf
# ============================================================
library(UCSCXenaTools)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(matrixStats)
library(ggplot2)
library(pheatmap)
library(grid)
library(scales)

# ----------------------------
# 0) USER SETTINGS
# ----------------------------
base_dir <- "FOLDER PATH/TCGA VALIDATION"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# LOW/HIGH sample lists
low_path  <- file.path(base_dir, "ER_ratio_LOW_samples.csv")
high_path <- file.path(base_dir, "ER_ratio_HIGH_samples.csv")
sample_id_col <- "sample_id"

# ATAC hub
atac_host <- "https://atacseq.xenahubs.net"

# Dataset selection rules
# - We'll pick BRCA cohort + Label starts with "All" (All peaks)
# - If you want enhancer/promoter-only, change label_prefix to "Enha" or "Prom"
label_prefix <- "All"

# Download settings
download_dir <- file.path(base_dir, "xena_atac_brca_download")
dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)

# Save downloaded matrix here
atac_rds_out <- file.path(base_dir, "TCGA_BRCA_ATAC_all_peaks_Log2Counts.rds")

# ProbeMap output
probemap_out <- file.path(base_dir, "brca_peak.probeMap")

# DA analysis parameters
collapse_by_case <- TRUE
collapse_fun <- "median"   # "median" recommended, or "mean"
min_non_na_per_group <- 10
fdr_cutoff <- 0.05
top_n_heatmap <- 200
set.seed(1)

# ----------------------------
# 1) ATAC-SEQ
# ----------------------------
data("XenaData", package = "UCSCXenaTools")

atac_hits <- XenaData %>%
  filter(XenaHosts == atac_host) %>%
  filter(str_detect(tolower(XenaDatasets), "peak")) %>%
  select(XenaHosts, XenaCohorts, XenaDatasets, Label, DataSubtype, SampleCount) %>%
  arrange(desc(SampleCount))

cat("\nTop ATAC peak datasets on the ATAC hub:\n")
print(atac_hits, n = 30)

# ----------------------------
# 2) PICK BRCA 'All peaks' DATASET
# ----------------------------
brca_ds <- XenaData %>%
  filter(XenaHosts == atac_host) %>%
  filter(str_detect(XenaCohorts, "Breast Cancer \\(BRCA\\)")) %>%
  filter(str_detect(Label, paste0("^", label_prefix))) %>%
  arrange(desc(SampleCount)) %>%
  slice(1)

if (nrow(brca_ds) == 0) {
  stop("No BRCA dataset matched your label_prefix. Inspect atac_hits and adjust label_prefix.")
}

target_dataset <- brca_ds$XenaDatasets[1]
cat("\nSelected dataset:\n")
print(brca_ds)

cat("\nUsing Xena dataset id:\n", target_dataset, "\n", sep="")

# ----------------------------
# 3) DOWNLOAD + PREPARE MATRIX
# ----------------------------
# Note: We use curl -L to avoid redirects.
xe <- XenaGenerate(subset = XenaHosts == atac_host & XenaDatasets == target_dataset)
xq <- XenaQuery(xe)

xd <- XenaDownload(
  xq,
  destdir = download_dir,
  method  = "curl",
  extra   = "-L",
  force   = TRUE
)

obj <- XenaPrepare(xd)

cat("\nPrepared object class:", class(obj), "\n")
cat("Prepared object dim:", paste(dim(obj), collapse=" x "), "\n")

# XenaPrepare usually returns a data.frame/tibble:
# For this hub, first column often stores peak id; sample columns follow.
# Your earlier code used obj$sample; we handle both common cases safely.
first_col_name <- colnames(obj)[1]
cat("First column name appears to be:", first_col_name, "\n")

peak_ids <- obj[[first_col_name]]
mat <- as.matrix(obj[, -1, drop = FALSE])
rownames(mat) <- peak_ids

cat("\nATAC matrix dim (peaks x samples):", paste(dim(mat), collapse=" x "), "\n")
print(mat[1:min(3,nrow(mat)), 1:min(3,ncol(mat)), drop=FALSE])

saveRDS(mat, file = atac_rds_out)
cat("\nSaved ATAC matrix RDS:\n", atac_rds_out, "\n", sep="")

# ----------------------------
# 4) DOWNLOAD PROBEMAP (PEAK COORDINATES)
# ----------------------------
# Many ATAC hub datasets use a probeMap file at:
#   https://atacseq.xenahubs.net/download/<cohort>/<cohort>_peak.probeMap
#
# We can ask UCSCXenaTools for the probeMap URL:
pm_url <- has_probeMap(atac_host, target_dataset, return_url = TRUE)
cat("\nprobeMap URL reported by has_probeMap():\n")
print(pm_url)

# If pm_url is NA/NULL, fall back to a common default guess.
if (length(pm_url) == 0 || all(is.na(pm_url))) {
  # typical fallback for BRCA
  pm_url <- paste0(atac_host, "/download/brca/brca_peak.probeMap")
  cat("probeMap URL fallback used:\n", pm_url, "\n")
}

download.file(pm_url, destfile = probemap_out, mode = "wb")
cat("Downloaded probeMap to:\n", probemap_out, "\n", sep="")

annot <- read_tsv(probemap_out, show_col_types = FALSE)
cat("\nprobeMap annotation dim:", paste(dim(annot), collapse=" x "), "\n")
cat("probeMap columns:\n")
print(colnames(annot))

# Expecting columns like: id, chrom, chromStart, chromEnd (common)
needed <- c("id","chrom","chromStart","chromEnd")
if (!all(needed %in% colnames(annot))) {
  stop("probeMap missing required columns (id, chrom, chromStart, chromEnd). Inspect annot.")
}

# Create a peak coordinate table matching matrix rownames
peak_annot <- annot[match(rownames(mat), annot$id), c("chrom","chromStart","chromEnd")]
peak_annot$start <- peak_annot$chromStart + 1  # BED->1-based
peak_annot$end   <- peak_annot$chromEnd
peak_annot <- peak_annot[, c("chrom","start","end")]
peak_annot$peak <- rownames(mat)

coord_out <- file.path(base_dir, "TCGA_BRCA_ATAC_peak_coordinates.csv")
write.csv(peak_annot, coord_out, row.names = FALSE)
cat("Saved peak coordinates:\n", coord_out, "\n", sep="")

# ----------------------------
# 5) LOAD LOW/HIGH LISTS
# ----------------------------
low_df  <- read.csv(low_path,  stringsAsFactors = FALSE)
high_df <- read.csv(high_path, stringsAsFactors = FALSE)

if (!sample_id_col %in% colnames(low_df))  stop(paste("Missing column", sample_id_col, "in", low_path))
if (!sample_id_col %in% colnames(high_df)) stop(paste("Missing column", sample_id_col, "in", high_path))

low_samples  <- unique(na.omit(as.character(low_df[[sample_id_col]])))
high_samples <- unique(na.omit(as.character(high_df[[sample_id_col]])))

low_cases  <- unique(substr(low_samples,  1, 12))
high_cases <- unique(substr(high_samples, 1, 12))

cat("\nLOW cases :", length(low_cases),  "\n")
cat("HIGH cases:", length(high_cases), "\n")

# ----------------------------
# 6) MATCH ATAC COLUMNS TO LOW/HIGH BY CASE
# ----------------------------
atac_cols  <- colnames(mat)
atac_cases <- substr(atac_cols, 1, 12)

keep_low_idx  <- which(atac_cases %in% low_cases)
keep_high_idx <- which(atac_cases %in% high_cases)

cat("\nATAC columns total:", length(atac_cols), "\n")
cat("Matched LOW columns :", length(keep_low_idx),  "\n")
cat("Matched HIGH columns:", length(keep_high_idx), "\n")

if (length(keep_low_idx) < 2 || length(keep_high_idx) < 2) {
  stop("Too few matched samples in one/both groups. Check sample IDs and TCGA barcode formatting.")
}

mat_low  <- mat[, keep_low_idx,  drop = FALSE]
mat_high <- mat[, keep_high_idx, drop = FALSE]

# Save mapping (audit trail)
mapping <- bind_rows(
  data.frame(
    atac_col  = colnames(mat_low),
    atac_case = substr(colnames(mat_low), 1, 12),
    group     = "LOW",
    stringsAsFactors = FALSE
  ),
  data.frame(
    atac_col  = colnames(mat_high),
    atac_case = substr(colnames(mat_high), 1, 12),
    group     = "HIGH",
    stringsAsFactors = FALSE
  )
) %>% arrange(group, atac_case, atac_col)

map_path <- file.path(base_dir, "ATAC_LOW_HIGH_column_mapping.csv")
write.csv(mapping, map_path, row.names = FALSE)
cat("Saved mapping:", map_path, "\n")

# ----------------------------
# 7) (OPTIONAL/RECOMMENDED) COLLAPSE MULTIPLE ALIQUOTS PER CASE
# ----------------------------
collapse_matrix_by_case <- function(m, fun=c("median","mean")) {
  fun <- match.arg(fun)
  cases <- substr(colnames(m), 1, 12)
  split_cols <- split(seq_along(cases), cases)
  
  collapsed <- sapply(names(split_cols), function(cs) {
    sub <- m[, split_cols[[cs]], drop=FALSE]
    if (fun == "median") matrixStats::rowMedians(sub, na.rm = TRUE) else rowMeans(sub, na.rm = TRUE)
  })
  
  collapsed <- as.matrix(collapsed)
  rownames(collapsed) <- rownames(m)
  collapsed
}

if (collapse_by_case) {
  mat_low_use  <- collapse_matrix_by_case(mat_low,  fun = collapse_fun)
  mat_high_use <- collapse_matrix_by_case(mat_high, fun = collapse_fun)
  
  saveRDS(mat_low_use,  file.path(base_dir, "TCGA_BRCA_ATAC_LOW_byCase.rds"))
  saveRDS(mat_high_use, file.path(base_dir, "TCGA_BRCA_ATAC_HIGH_byCase.rds"))
  
  cat("\nCollapsed LOW  dim:", paste(dim(mat_low_use),  collapse=" x "), "\n")
  cat("Collapsed HIGH dim:", paste(dim(mat_high_use), collapse=" x "), "\n")
} else {
  mat_low_use  <- mat_low
  mat_high_use <- mat_high
}

# ----------------------------
# 8) DIFFERENTIAL ACCESSIBILITY (Wilcoxon per peak + BH FDR)
# ----------------------------
run_da_wilcox <- function(low_mat, high_mat, min_n = 10) {
  stopifnot(identical(rownames(low_mat), rownames(high_mat)))
  
  n_low  <- rowSums(!is.na(low_mat))
  n_high <- rowSums(!is.na(high_mat))
  keep <- (n_low >= min_n) & (n_high >= min_n)
  
  peaks <- rownames(low_mat)
  
  low_mean  <- rowMeans(low_mat,  na.rm = TRUE)
  high_mean <- rowMeans(high_mat, na.rm = TRUE)
  low_median  <- matrixStats::rowMedians(low_mat,  na.rm = TRUE)
  high_median <- matrixStats::rowMedians(high_mat, na.rm = TRUE)
  
  logFC <- high_mean - low_mean
  delta_median <- high_median - low_median
  
  pvals <- rep(NA_real_, length(peaks))
  idx <- which(keep)
  
  cat("\nPeaks passing min_n filter:", length(idx), " / ", length(peaks), "\n")
  
  for (i in idx) {
    x <- as.numeric(high_mat[i, ]); x <- x[!is.na(x)]
    y <- as.numeric(low_mat[i, ]);  y <- y[!is.na(y)]
    pvals[i] <- suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
  }
  
  df <- data.frame(
    peak = peaks,
    n_low = n_low,
    n_high = n_high,
    mean_low = low_mean,
    mean_high = high_mean,
    median_low = low_median,
    median_high = high_median,
    logFC = logFC,
    delta_median = delta_median,
    p_value = pvals,
    stringsAsFactors = FALSE
  )
  df$fdr <- p.adjust(df$p_value, method = "BH")
  df
}

da_tbl <- run_da_wilcox(mat_low_use, mat_high_use, min_n = min_non_na_per_group) %>%
  mutate(direction = case_when(
    is.na(logFC) ~ NA_character_,
    logFC > 0 ~ "HIGH_up",
    logFC < 0 ~ "LOW_up",
    TRUE ~ "no_change"
  )) %>%
  arrange(fdr, p_value)

da_out <- file.path(base_dir, "ATAC_DA_peaks_High_vs_Low.csv")
write.csv(da_tbl, da_out, row.names = FALSE)
cat("\nSaved DA table:", da_out, "\n")

# Attach coordinates
da_with_coords <- da_tbl %>%
  left_join(peak_annot, by = "peak") %>%
  relocate(chrom, start, end, .after = peak)

da_coord_out <- file.path(base_dir, "ATAC_DA_peaks_High_vs_Low_with_coords.csv")
write.csv(da_with_coords, da_coord_out, row.names = FALSE)
cat("Saved DA table with coords:", da_coord_out, "\n")

da_sig <- da_with_coords %>% filter(!is.na(fdr), fdr < fdr_cutoff)
da_sig_out <- file.path(base_dir, "ATAC_DA_peaks_High_vs_Low_significant.csv")
write.csv(da_sig, da_sig_out, row.names = FALSE)
cat("Saved significant DA peaks:", da_sig_out, "\n")
cat("Significant peaks (FDR<", fdr_cutoff, "): ", nrow(da_sig), "\n", sep="")

# ----------------------------
# 9) PLOTS (PDF)
# ----------------------------
# A) Volcano
volcano_out <- file.path(base_dir, "ATAC_volcano_HIGH_vs_LOW.pdf")

vol <- da_tbl %>%
  mutate(minus_log10_fdr = -log10(pmax(fdr, 1e-300)),
         sig = !is.na(fdr) & fdr < fdr_cutoff)

p_vol <- ggplot(vol, aes(x = logFC, y = minus_log10_fdr)) +
  geom_point(aes(alpha = sig), size = 1.2) +
  labs(
    title = "TCGA-BRCA ATAC-seq Differential Accessibility",
    subtitle = paste0("HIGH vs LOW (Wilcoxon per peak; BH-FDR<", fdr_cutoff, ")"),
    x = "log2(Accessibility) difference (HIGH - LOW)",
    y = "-log10(FDR)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(volcano_out, p_vol, width = 6, height = 5)
cat("Saved volcano:", volcano_out, "\n")

# B) Heatmap of top DA peaks
heat_out <- file.path(base_dir, "ATAC_heatmap_topDA.pdf")

top_peaks <- da_tbl %>%
  filter(!is.na(fdr)) %>%
  arrange(fdr) %>%
  slice_head(n = min(top_n_heatmap, n())) %>%
  pull(peak)

if (length(top_peaks) >= 2) {
  mat_all <- cbind(mat_low_use, mat_high_use)
  col_group <- c(rep("LOW", ncol(mat_low_use)), rep("HIGH", ncol(mat_high_use)))
  
  ann <- data.frame(Group = factor(col_group, levels = c("LOW","HIGH")))
  rownames(ann) <- colnames(mat_all)
  
  hm_mat <- mat_all[top_peaks, , drop=FALSE]
  hm_mat_z <- t(scale(t(hm_mat)))
  
  ph <- pheatmap(
    hm_mat_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    annotation_col = ann,
    border_color = NA,
    main = paste0("Top DA peaks (n=", length(top_peaks), ") — HIGH vs LOW")
  )
  
  pdf(heat_out, width = 10, height = 7)
  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()
  
  cat("Saved heatmap:", heat_out, "\n")
} else {
  cat("Skipped heatmap: not enough peaks with valid FDR.\n")
}

# C) PCA
pca_out <- file.path(base_dir, "ATAC_PCA_LOW_vs_HIGH.pdf")

mat_all <- cbind(mat_low_use, mat_high_use)
col_group <- c(rep("LOW", ncol(mat_low_use)), rep("HIGH", ncol(mat_high_use)))

vars <- matrixStats::rowVars(mat_all, na.rm = TRUE)
keep_var <- order(vars, decreasing = TRUE)[1:min(5000, length(vars))]

pca <- prcomp(t(mat_all[keep_var, , drop=FALSE]), center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = factor(col_group, levels = c("LOW","HIGH"))
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  labs(
    title = "ATAC-seq PCA (TCGA-BRCA matched samples)",
    subtitle = if (collapse_by_case) "Collapsed by case (median across aliquots)" else "Aliquot-level",
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic()

ggsave(pca_out, p_pca, width = 6, height = 5)
cat("Saved PCA:", pca_out, "\n")

# ----------------------------
# 10) SESSION INFO
# ----------------------------
cat("\nDONE.\n")
print(sessionInfo())
