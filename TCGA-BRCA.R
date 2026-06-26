library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(colorRamps)
library(ggplot2)
library(NbClust)
library(gridExtra)
library(circlize)
library(Rtsne)
library(gplots)
library(Exact)
library(ClusterR)
library(stats)
library(amap)
library(Hmisc)
library(ggpubr)
library(data.table)
library(plotrix)
library(infercnv)
library(anndata)
library(survival)
library(survminer)
library(grid)
options(scipen = 100)

# Load the transposed transcriptomic data from the CSV file
merged_data <- read.csv("downloaded_normalized_breast_cancer_TCGA_Project_with_Metadata_annotations")


# Print the first few rows to verify
View(merged_data[1:10,1:10])

#######robustness analysis ER negative DEG#####
# Full ER-negative gene list
er_negative_genes_full <- c(
  "PTMA","FTH1","TMSB10","FAU","TMSB4X","EEF1A1","UBA52","FTL","NACA",
  "H3F3A","SERF2","EIF1","TPT1","GAPDH","MYL6","EEF1D","MIF","B2M","OAZ1",
  "ATP5F1E","COX7C","H3F3B","HSPB1","FOS","COX6C","CALM2","PFDN5","RACK1",
  "SRP14","PPDPF","TMA7","BTF3","JUN","CHCHD2","CFL1","S100A11","HSP90AB1",
  "UBB","NDUFA4","COX6A1","TOMM7","UQCR11","PFN1","PPIA","COX6B1","UBL5",
  "UQCRQ","COX4I1","UBC","UQCRB","PSMA7","OST4","IFITM3","COX7A2","JUNB",
  "HSP90AA1","DYNLL1","EEF1B2","S100A10","LGALS1","CD63","TXN","HINT1",
  "CRIP1","PRDX1","S100A6","SEC61G","COX8A","POLR2L","KRT18","ELOB",
  "TMEM258","HSPA1A","COMMD6","COX5B","SUMO2","HSPE1","SOD1","SLC25A6",
  "UQCR10","HMGB1","ATP5MC2","SSR4","MT2A","SNRPD2","GPX4","HNRNPA1",
  "IER2","NDUFA11","PRDX2","DUSP1","TPI1","YBX1","NEDD8","NDUFA13",
  "WDR83OS","ATP5MD","CSTB","EDF1","NDUFC2","GUK1","ATP5ME","SH3BGRL3",
  "UQCRH","ZFP36","LY6E","NDUFB4","TSPO","FOSB","SKP1","IFITM2","NDUFS5",
  "NPM1","HSPA1B","NDUFB2","EEF2","GNG5","BTG1","SAT1","MYL12B","ZFAS1",
  "NDUFA1","TSTD1","TRMT112","CALM1","ATP5MG","ATP5MF","PABPC1","H2AFZ",
  "POMP","ATP5MPL","KRT19","DBI","ADIRF","NUPR1","KRT8","MZT2B","PRDX5",
  "SEC61B","EGR1","COX7B","MYL12A","C19orf70","ROMO1","NDUFB1","NDUFB11",
  "BRK1","CLIC1","NPC2","CTSD","GADD45B","RHOB","PSME1","EIF3K","FKBP1A",
  "ATP6V1G1","MARCKSL1","RAN","NOP10","SMIM22","SNRPE","ATP5PF","SPINT2",
  "PPIB","ANAPC11","KLF6","TUBA1B","NDUFB9","SELENOH","SLC25A5","NDUFC1",
  "LAMTOR5","TAGLN2","ITM2B","ERH","C19orf53","CD9","COX17","SAP18",
  "PET100","CUTA","HSPA8","NDUFB10","PPP1R15A","NHP2","KRTCAP2","BST2",
  "ARL6IP1","VAMP8","KDELR2","C12orf57","ATP6V1F","RABAC1","EIF3E","SNRPG",
  "C8orf59","HLA-B","MINOS1","POLR2I","S100A14","LDHA","HLA-C","CRABP2",
  "HMGN1","RBX1","NDUFA2","ATP5PO","ATP5MC3","RARRES3","ATOX1","NDUFB7",
  "DNAJB1","CLDN4","COPS9","FXYD3","PLA2G16","CD59","C6orf48","S100A16",
  "COX5A","PRELID1","HMGN2","TFF3","ATF3","BANF1","LSM4","TMED2","NDUFA3",
  "CYB5A","VIM","NCL","LAPTM4A","LGALS3","PSME2","RNF181","EIF3F","ALDOA",
  "ATP6V0B","HIST1H1C","CST3","EIF5A","MT1X","C4orf3","SF3B6","IFI27",
  "FAM96B","APRT","SCGB2A2","MTRNR2L12","HIGD2A","SLC25A3","EIF3G",
  "HSBP1","FKBP2","AURKAIP1","AP2S1","POLR2K","ARPC3","TIMM13","SNU13",
  "RHOC","PCBP1","GADD45GIP1","SLIRP","SNRPF","CD74","DDIT4","CRIP2",
  "CD24","IFITM1","NDUFB8","ENY2","LAMTOR2","DDT","ID2","HLA-A","DRAP1",
  "MRPL51","KRT10","ANXA2","CSNK2B","HLA-E","CNPY2","CYCS","GNB2L1",
  "EIF6","SNHG8","NFKBIA","SPCS1","SNRPB","NDUFAF3","TMEM147","IFI6",
  "SSR2","PEBP1","TUBB","RBM8A","TMBIM4","FXYD5","MUCL1"
)

results <- data.frame(
  Num_Genes = integer(),
  Cutoff_Type = character(),
  Cutoff_Value = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# To track which sizes have no significant p-values
no_significant_sets <- c()

# Loop through gene set sizes
for (i in 50:length(er_negative_genes_full)) {
  genes_subset <- er_negative_genes_full[1:i]
  valid_genes <- intersect(genes_subset, colnames(merged_data))
  if (length(valid_genes) < 10) next
  
  merged_data$er_neg_score <- rowMeans(merged_data[, valid_genes, drop = FALSE], na.rm = TRUE)
  luminal <- subset(merged_data, merged_data$er_status_by_ihc == "Positive")
  surv_object <- Surv(time = as.numeric(luminal$OS_time), event = as.numeric(luminal$OS_status))
  
  has_significant <- FALSE
  
  for (p in c(0.50, 0.60, 0.70, 0.75)) {
    cutoff <- quantile(luminal$er_neg_score, probs = p, na.rm = TRUE)
    luminal$group <- ifelse(luminal$er_neg_score > cutoff, "High", "Low")
    
    fit <- survfit(surv_object ~ group, data = luminal)
    surv_test <- survdiff(surv_object ~ group, data = luminal)
    pval <- 1 - pchisq(surv_test$chisq, df = 1)
    
    if (!is.na(pval) && pval < 0.05) has_significant <- TRUE
    
    results <- rbind(results, data.frame(
      Num_Genes = i,
      Cutoff_Type = paste0("Q", p * 100),
      Cutoff_Value = cutoff,
      P_Value = pval
    ))
  }
  
  if (!has_significant) {
    no_significant_sets <- c(no_significant_sets, i)
  }
}

# Save to CSV
write.csv(results, "KM_ER_negative_gene_score_results.csv", row.names = FALSE)

# Report gene counts with no significant p-values
cat("Gene set sizes with no significant p-values (< 0.05):\n")
print(no_significant_sets)


#######robustness analysis ER positive DEG#####
# === ER-positive genes ===
er_positive_genes_full <- c(
  "ESR1","TRPS1","THSD4","ERBB4","ANKRD30A","AFF3","PATJ","PBX1","LRBA","TTC6","MSI2",
  "FOXP1","RBM47","RERE","MED13L","ANKS1B","FTX","AGAP1","RORA","DACH1","NEBL","ASH1L",
  "APBB2","KIAA1217","ZFAND3","MAST4","MAGI1","SIPA1L3","MACF1","NFIA","IGF1R","NAALADL2",
  "LPP","AUTS2","TRIO","SHANK2","LDLRAD4","ZBTB20","MAGI3","BCL2","XIST","VPS13B","BMPR1B",
  "PVT1","TBC1D9","EXOC4","WWOX","SPIDR","TRAPPC9","WWP1","TBC1D5","MLPH","SRPK2","KDM4B",
  "MRTFB","PTK2","KAZN","HIVEP3","SMYD3","GALNT7","CERS6","ARID1B","NEDD4L","GSE1","PRLR",
  "RABEP1","NCOA2","ITPR2","GLI3","EVL","EIF4G3","NSMCE2","FAF1","IMMP2L","ADK","DOCK1",
  "JMJD1C","TCF12","KMT2C","SYTL2","USP34","GRHL2","SH3RF1","GOLGA4","SIK3","PTPRT","KAT6B",
  "FARP1","C1orf21","BIRC6","CASZ1","MYO6","PDXDC1","DNAH14","RAD51B","LIMCH1","SHROOM3",
  "ATXN1","FUT8","TULP4","PARD3","STAG1","ALCAM","MBD5","CACNA1D","ZNF609","ETV6","BRWD1",
  "TMCC1","RERG","FHIT","DNAJC1","TNRC6A","MIRLET7BHG","PPP6R3","RBM6","NF1","GATA3",
  "LDLRAD3","EHMT1","KIAA1324","FLNB","FBXL17","CELSR1","FAM172A","CUX1","UBR5","BCAS3",
  "ZMYND8","INPP4B","ARHGAP32","TTC3","KCNMA1","SAMD12","DENND1B","ANKRD17","PKP4","CTBP2",
  "ZCCHC7","KDM2A","DLEU2","TNRC18","ZNF91","SETD5","FBXL7","NIPBL","EXOC6B","DLG1","BTBD9",
  "GMDS-DT","ZNF638","FNDC3B","TANC2","KANSL1","EXT1","RALGPS2","BAZ2B","STRBP","DTNB",
  "ANK3","GTF2I","ERC1","FARS2","BTRC","MKLN1","SSBP2","VAV3","PLEKHA7","ARHGEF3","GPHN",
  "CSNK1A1","NCKAP5","ERBIN","AR","CBFA2T2","TBL1XR1","ARID4B","MAML3","POGZ","VWA8","CLASP2",
  "PAN3","ENAH","FAM214A","MAP7","MIR29B2CHG","IKBKB","MUC20-OT1","COP1","CADPS2",
  "C1QTNF3-AMACR","ATXN2","BPTF","NCOA1","PUM2","NSD1","TENT2","PTPRK","RALGAPA1","ZNF148",
  "GPATCH2","NTN4","SLC39A11","PACS1","GPATCH8","DIP2B","ATP8B1","CDK13","UBE2H","AMBRA1"
)

# === Results dataframe ===
results <- data.frame(
  Num_Genes = integer(),
  Cutoff_Type = character(),
  Cutoff_Value = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

no_significant_sets <- c()

# === Loop through gene set sizes ===
for (i in 50:length(er_positive_genes_full)) {
  genes_subset <- er_positive_genes_full[1:i]
  valid_genes <- intersect(genes_subset, colnames(merged_data))
  if (length(valid_genes) < 10) next
  
  merged_data$er_pos_score <- rowMeans(merged_data[, valid_genes, drop = FALSE], na.rm = TRUE)
  luminal <- subset(merged_data, merged_data$er_status_by_ihc == "Positive")
  surv_object <- Surv(time = as.numeric(luminal$OS_time), event = as.numeric(luminal$OS_status))
  
  has_significant <- FALSE
  
  for (p in c(0.50, 0.60, 0.70, 0.75)) {
    cutoff <- quantile(luminal$er_pos_score, probs = p, na.rm = TRUE)
    luminal$group <- ifelse(luminal$er_pos_score > cutoff, "High", "Low")
    
    fit <- survfit(surv_object ~ group, data = luminal)
    surv_test <- survdiff(surv_object ~ group, data = luminal)
    pval <- 1 - pchisq(surv_test$chisq, df = 1)
    
    if (!is.na(pval) && pval < 0.05) has_significant <- TRUE
    
    results <- rbind(results, data.frame(
      Num_Genes = i,
      Cutoff_Type = paste0("Q", p * 100),
      Cutoff_Value = cutoff,
      P_Value = pval
    ))
  }
  
  if (!has_significant) {
    no_significant_sets <- c(no_significant_sets, i)
  }
}

# Save to CSV
write.csv(results, "KM_ER_positive_gene_score_results.csv", row.names = FALSE)

# Report gene counts with no significant p-values
cat("Gene set sizes with no significant p-values (< 0.05):\n")
print(no_significant_sets)

### Cox proportional hazards regression-ER positive score###

er_positive_genes_score <- c("ESR1", "TRPS1", "THSD4", "ERBB4", "ANKRD30A", "AFF3", "PATJ", "PBX1", "LRBA", "TTC6",
                             "MSI2", "FOXP1", "RBM47", "RERE", "MED13L", "ANKS1B", "FTX", "AGAP1", "RORA", "DACH1",
                             "NEBL", "ASH1L", "APBB2", "KIAA1217", "ZFAND3", "MAST4", "MAGI1", "SIPA1L3", "MACF1", "NFIA",
                             "IGF1R", "NAALADL2", "LPP", "AUTS2", "TRIO", "SHANK2", "LDLRAD4", "ZBTB20", "MAGI3", "BCL2",
                             "XIST", "VPS13B", "BMPR1B", "PVT1", "TBC1D9", "EXOC4", "WWOX", "SPIDR", "TRAPPC9", "WWP1",
                             "TBC1D5", "MLPH", "SRPK2", "KDM4B", "MRTFB", "PTK2", "KAZN", "HIVEP3", "SMYD3", "GALNT7",
                             "CERS6", "ARID1B", "NEDD4L", "GSE1", "PRLR", "RABEP1", "NCOA2", "ITPR2", "GLI3", "EVL",
                             "EIF4G3", "NSMCE2", "FAF1", "IMMP2L", "ADK", "DOCK1", "JMJD1C", "TCF12", "KMT2C", "SYTL2",
                             "USP34", "GRHL2", "SH3RF1", "GOLGA4", "SIK3", "PTPRT", "KAT6B", "FARP1", "C1orf21", "BIRC6",
                             "CASZ1", "MYO6", "PDXDC1", "DNAH14", "RAD51B", "LIMCH1", "SHROOM3", "ATXN1", "FUT8", "TULP4",
                             "PARD3", "STAG1", "ALCAM", "MBD5", "CACNA1D", "ZNF609", "ETV6", "BRWD1", "TMCC1", "RERG",
                             "FHIT", "DNAJC1", "TNRC6A", "MIRLET7BHG", "PPP6R3", "RBM6", "NF1", "GATA3", "LDLRAD3",
                             "EHMT1", "KIAA1324", "FLNB", "FBXL17", "CELSR1", "FAM172A", "CUX1", "UBR5", "BCAS3", "ZMYND8",
                             "INPP4B", "ARHGAP32", "TTC3", "KCNMA1", "SAMD12", "DENND1B", "ANKRD17", "PKP4", "CTBP2",
                             "ZCCHC7", "KDM2A", "DLEU2", "TNRC18", "ZNF91", "SETD5", "ZNF638", "FNDC3B", "TANC2", "KANSL1", "EXT1", "RALGPS2", "BAZ2B")


# Make sure 'er_negative_genes_score' is numeric
luimnal$er_genes_score <- as.numeric(luimnal$er_positive_genes_score)

# Format covariates
luimnal$age_at_diagnosis <- as.numeric(luimnal$age_at_diagnosis)

luimnal$margin_status <- as.factor(luimnal$margin_status)
luimnal$ajcc_pathologic_tumor_stage <- as.factor(luimnal$ajcc_pathologic_tumor_stage)
luimnal$her2_fish_status <- as.factor(luimnal$her2_fish_status)
#luimnal$menapause_status <- as.factor(luimnal$menapause_status)
luimnal$ethnicity <- as.factor(luimnal$ethnicity)

# Ensure survival data is ready
luimnal$OS_time <- as.numeric(luimnal$OS_time)
luimnal$OS_status <- as.numeric(luimnal$OS_status)

# Create the survival object
surv_object <- Surv(time = luimnal$OS_time, event = luimnal$OS_status)

# Cox model
cox_model <- coxph(surv_object ~ er_positive_genes_score + 
                     age_at_diagnosis +
                     margin_status +
                     ajcc_pathologic_tumor_stage +
                     her2_fish_status +
                     ethnicity,
                   data = luimnal)

# View results
summary(cox_model)


# Create PDF file
pdf("cox_model_summary.pdf", width = 8.5, height = 11)

# Capture model summary as text
model_summary <- capture.output(summary(cox_model))

# Combine and render text in the PDF
grid.newpage()
grid.text(paste(model_summary, collapse = "\n"), x = 0.05, y = 0.95, just = c("left", "top"))

# Close the PDF device
dev.off()


### Cox proportional hazards regression-ER negative score###
er_negative_genes_score <- c("PTMA","FTH1","TMSB10","FAU","TMSB4X","EEF1A1","UBA52","FTL","NACA",
                             "H3F3A","SERF2","EIF1","TPT1","GAPDH","MYL6","EEF1D","MIF","B2M","OAZ1",
                             "ATP5F1E","COX7C","H3F3B","HSPB1","FOS","COX6C","CALM2","PFDN5","RACK1",
                             "SRP14","PPDPF","TMA7","BTF3","JUN","CHCHD2","CFL1","S100A11","HSP90AB1",
                             "UBB","NDUFA4","COX6A1","TOMM7","UQCR11","PFN1","PPIA","COX6B1","UBL5",
                             "UQCRQ","COX4I1","UBC","UQCRB")

# Make sure 'er_negative_genes_score' is numeric
luimnal$er_genes_score <- as.numeric(luimnal$er_negative_genes_score)

# Format covariates
luimnal$age_at_diagnosis <- as.numeric(luimnal$age_at_diagnosis)

luimnal$margin_status <- as.factor(luimnal$margin_status)
luimnal$ajcc_pathologic_tumor_stage <- as.factor(luimnal$ajcc_pathologic_tumor_stage)
luimnal$her2_fish_status <- as.factor(luimnal$her2_fish_status)
#luimnal$menapause_status <- as.factor(luimnal$menapause_status)
luimnal$ethnicity <- as.factor(luimnal$ethnicity)

# Ensure survival data is ready
luimnal$OS_time <- as.numeric(luimnal$OS_time)
luimnal$OS_status <- as.numeric(luimnal$OS_status)

# Create the survival object
surv_object <- Surv(time = luimnal$OS_time, event = luimnal$OS_status)

# Cox model
cox_model <- coxph(surv_object ~ er_negative_genes_score + 
                     age_at_diagnosis +
                     margin_status +
                     ajcc_pathologic_tumor_stage +
                     her2_fish_status +
                     ethnicity,
                   data = luimnal)

# View results
summary(cox_model)



# Create PDF file
pdf("cox_model_summary_er_negative.pdf", width = 8.5, height = 11)

# Capture model summary as text
model_summary <- capture.output(summary(cox_model))

# Combine and render text in the PDF
grid.newpage()
grid.text(paste(model_summary, collapse = "\n"), x = 0.05, y = 0.95, just = c("left", "top"))

# Close the PDF device
dev.off()

#######empirical p-value calculation#########
# Filter ER-positive samples
luimnal <- subset(merged_data, merged_data$er_status_by_ihc == "Positive")

# Ensure survival columns are numeric
luimnal$OS_time <- as.numeric(luimnal$OS_time)
luimnal$OS_status <- as.numeric(luimnal$OS_status)

# Genes available for analysis (exclude metadata columns)
all_genes <- colnames(merged_data)
exclude_cols <- c("OS_time", "OS_status", "er_status_by_ihc", "er_positive_genes", "er_pos_group")
gene_pool <- setdiff(intersect(all_genes, colnames(merged_data)), exclude_cols)

# Run 100 iterations
for (i in 1:n_iterations) {
  # 1. Randomly select 50 genes
  selected_genes <- sample(gene_pool, n_genes)
  # 2. Calculate gene score as mean expression (ensure numeric)
  gene_matrix <- sapply(luimnal[, selected_genes, drop = FALSE], as.numeric)
  luimnal$random_signature_score <- rowMeans(gene_matrix, na.rm = TRUE)
  
  # 3. Divide into high/low based on 25th percentile
  cutoff <- quantile(luimnal$random_signature_score, 0.25, na.rm = TRUE)
  luimnal$signature_group <- ifelse(luimnal$random_signature_score > cutoff, "High", "Low")
  
  # 4. Fit Kaplan-Meier model
  surv_obj <- Surv(luimnal$OS_time, luimnal$OS_status)
  fit <- survfit(surv_obj ~ signature_group, data = luimnal)
  
  # 5. Log-rank test p-value
  logrank_test <- survdiff(surv_obj ~ signature_group, data = luimnal)
  pval <- 1 - pchisq(logrank_test$chisq, df = 1)
  pvals[i] <- pval
  
  # 6. Count if significant
  if (pval < pval_cutoff) {
    significant_count <- significant_count + 1
  }
}

# Output results
cat("Number of significant p-values (< 0.05):", significant_count, "out of", n_iterations, "\n")
hist(pvals, breaks = 20, main = "Distribution of p-values from 100 Random Signatures", xlab = "p-value")

#######KM curves:
##################ER Negative score
er_negative_genes <- c("PTMA", "FTH1", "TMSB10", "FAU", "TMSB4X", "EEF1A1", "UBA52", "FTL", "NACA",
                       "H3F3A", "SERF2", "EIF1", "TPT1", "GAPDH", "MYL6", "EEF1D", "MIF", "B2M",
                       "OAZ1", "ATP5F1E", "COX7C", "H3F3B", "HSPB1", "FOS", "COX6C", "CALM2",
                       "PFDN5", "RACK1", "SRP14", "PPDPF", "TMA7", "BTF3", "JUN", "CHCHD2", "CFL1",
                       "S100A11", "HSP90AB1", "UBB", "NDUFA4", "COX6A1", "TOMM7", "UQCR11", "PFN1",
                       "PPIA", "COX6B1", "UBL5", "UQCRQ", "COX4I1", "UBC", "UQCRB")

#Ensure only genes present in merged_data are used for ER- score
valid_er_neg_genes <- intersect(er_negative_genes, colnames(merged_data))

#Compute `er_negative_genes_score`
merged_data$er_negative_genes_score <- rowMeans(merged_data[, valid_er_neg_genes, drop = FALSE], na.rm = TRUE)

#Ensure 'er_negative_genes_score' is numeric
luimnal$er_negative_genes_score <- as.numeric(luimnal$er_negative_genes_score)

#Define High vs. Low Groups (Threshold> 25%)
luimnal$er_neg_group <- ifelse(luimnal$er_negative_genes_score > 9.5, "High", "Low")

#Ensure survival data is numeric
luimnal$OS_time <- as.numeric(luimnal$OS_time)
luimnal$OS_status <- as.numeric(luimnal$OS_status)

#Create Kaplan-Meier Survival Object
surv_object <- Surv(time = luimnal$OS_time, event = luimnal$OS_status)

#Fit survival model
km_fit <- survfit(surv_object ~ er_neg_group, data = luimnal)

# Plot Kaplan-Meier curve
km_plot <- ggsurvplot(
  km_fit, data = luimnal,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE,
  xlab = "Time (Days)", ylab = "Survival Probability",
  title = "Kaplan-Meier Survival Curve: ER-Negative Gene Score",
  legend.title = "ER-Negative Score Group",
  legend.labs = c("Low", "High"),
  palette = custom_palette,
  ggtheme = theme_bw(base_size = 8)
)


er_positive_genes <- c("TRPS1", "THSD4", "ERBB4", "ANKRD30A", "AFF3", "PATJ", "PBX1", "LRBA", "TTC6",
                       "MSI2", "FOXP1", "RBM47", "RERE", "MED13L", "ANKS1B", "FTX", "AGAP1", "RORA", "DACH1",
                       "NEBL", "ASH1L", "APBB2", "KIAA1217", "ZFAND3", "MAST4", "MAGI1", "SIPA1L3", "MACF1",
                       "NFIA", "IGF1R", "NAALADL2", "LPP", "AUTS2", "TRIO", "SHANK2", "LDLRAD4", "ZBTB20",
                       "MAGI3", "BCL2", "XIST", "VPS13B", "BMPR1B", "PVT1", "TBC1D9", "EXOC4", "WWOX", "SPIDR",
                       "TRAPPC9", "WWP1", "TBC1D5", "MLPH", "SRPK2", "KDM4B", "MRTFB", "PTK2", "KAZN", "HIVEP3",
                       "SMYD3", "GALNT7", "CERS6", "ARID1B", "NEDD4L", "GSE1", "PRLR", "RABEP1", "NCOA2",
                       "ITPR2", "GLI3", "EVL", "EIF4G3", "NSMCE2", "FAF1", "IMMP2L", "ADK", "DOCK1", "JMJD1C",
                       "TCF12", "KMT2C", "SYTL2", "USP34", "GRHL2", "SH3RF1", "GOLGA4", "SIK3", "PTPRT", "KAT6B",
                       "FARP1", "C1orf21", "BIRC6", "CASZ1", "MYO6", "PDXDC1", "DNAH14", "RAD51B", "LIMCH1",
                       "SHROOM3", "ATXN1", "FUT8", "TULP4", "PARD3", "STAG1", "ALCAM", "MBD5", "CACNA1D",
                       "ZNF609", "ETV6", "BRWD1", "TMCC1", "RERG", "FHIT", "DNAJC1", "TNRC6A", "MIRLET7BHG",
                       "PPP6R3", "RBM6", "NF1", "GATA3", "LDLRAD3", "EHMT1", "KIAA1324", "FLNB", "FBXL17",
                       "CELSR1", "FAM172A", "CUX1", "UBR5", "BCAS3", "ZMYND8", "INPP4B", "ARHGAP32", "TTC3",
                       "KCNMA1", "SAMD12", "DENND1B", "ANKRD17", "PKP4", "CTBP2", "ZCCHC7", "KDM2A", "DLEU2")
valid_er_positive_genes <- intersect(er_positive_genes, colnames(merged_data))
merged_data$er_positive_genes <- rowMeans(merged_data[, valid_er_positive_genes, drop = FALSE], na.rm = TRUE)
luimnal <- subset(merged_data, merged_data$er_status_by_ihc == "Positive")
luimnal$er_positive_genes <- as.numeric(luimnal$er_positive_genes)
#Define High vs. Low Groups (Threshold> 25%)
luimnal$er_pos_group <- ifelse(luimnal$er_positive_genes > 9, "High", "Low")
luimnal$OS_time <- as.numeric(luimnal$OS_time)
luimnal$OS_status <- as.numeric(luimnal$OS_status)
surv_object <- Surv(time = luimnal$OS_time, event = luimnal$OS_status)
km_fit <- survfit(surv_object ~ er_pos_group, data = luimnal)
# 6️⃣ Plot Kaplan-Meier Curve
ggsurvplot(km_fit, data = luimnal,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, 
           xlab = "Time (Days)", ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve: ER-positive Gene Score",
           legend.title = "ER-positive Score Group",
           legend.labs = c("Low", "High"),
           palette = c("red", "blue"))

















