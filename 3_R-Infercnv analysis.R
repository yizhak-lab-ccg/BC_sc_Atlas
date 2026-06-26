#Infercnv code
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(colorRamps)
library(ggplot2)
library(NbClust)
library(gridExtra)
library(circlize)
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
options(scipen = 100)


link <- "#path to counts.csv"
p1 <- fread(link, header = TRUE) %>% as.matrix(rownames = TRUE)
p1 <- as.data.frame(p1)

# Define the genes of interest and the threshold
genes_of_interest <- c("PTPRC","CD3E", "CD8A", "CD8B", "CD19", "MS4A1", "CD14", "CD163")
threshold <- 0.5

# Create a new column indicating if any gene is above the threshold
p1 <- p1 %>%
  mutate(type = ifelse(rowSums(select(., any_of(genes_of_interest)) > threshold) > 0, "IMMUNE", "NOT_IMMUNE"))

f<-subset(p1,select=type)
p1<-p1[,-ncol(p1)]
p1<-t(p1)


#

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=p1,
                                    annotations_file=f,
                                    gene_order_file = read.table(#gene_order_file.txt"),
                                    ref_group_names=c("IMMUNE"))
# Set output path to 
out_dir <- "output file"
                                    
# Run infercnv
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              out_dir = out_dir,
                              cluster_by_groups = TRUE,
                              plot_steps = FALSE,
                              denoise = TRUE,
                              HMM = TRUE,
                              no_prelim_plot = TRUE,
                              png_res = 200,
                              analysis_mode = "subclusters",
                              write_expr_matrix = TRUE)

                                    
# Generate subcluster annotation plot
t <- plot_subclusters(infercnv_obj,
                      out_dir = out_dir,
                      output_filename = "subclusters_as_annotations")
                                    

#Following the subcluster INFERCNV possess each cell was assigned to a subclusters (based on the INFERCNV object):
#Non-immune samples without CNV were annotated as stroma cell.
#Immune cells with CNV were deleted.
#Non-immune cells with CNV were annotated as tumor cells.
#The processed h5ad files were subsetted to immune/tumor stroma based on the cells INFERCNV subcluster group.