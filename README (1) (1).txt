Breast Cancer Single-Cell Atlas: Analysis Code and Data

This repository contains the computational workflows used to generate and analyze the integrated breast cancer single-cell atlas described in this study. The repository includes preprocessing, cell-type integration and subclustering, tumor-cell analyses, copy-number variation analysis, cNMF analysis, ligand–receptor analysis, and validation-cohort analyses.

Datasets and data availability

The datasets included in the integrated breast cancer single-cell atlas are described in Supplementary Table 1 of the manuscript. For each dataset, Supplementary Table 1 provides the original study and source, links to the publicly available raw data, and the number of cells retained from each dataset after preprocessing and quality control and included in the atlas.

The individual datasets used to construct the atlas are available through the Breastscape website:

https://breastscape.twix.technion.ac.il/

Computational environment and package versions

Analyses were performed using the following Python and R environments. These package versions are provided to facilitate reproducibility of the computational workflow.

Python (version 3.13.5):

anndata 0.12.2
bbknn 1.6.0
cnmf 1.7.0
liana 1.6.1
lifelines 0.30.0
matplotlib 3.10.0
networkx 3.4.2
numpy 2.1.3
pandas 2.2.3
scanpy 1.11.4
scipy 1.15.3
scrublet 0.2.3
seaborn 0.13.2
statsmodels 0.14.5

Scrublet installation requires Microsoft C++ Build Tools on Windows. These are available from the Microsoft Visual C++ Build Tools website.

R (version 4.3.3; RStudio version 2025.5.1.513):

amap 0.8.20
anndata 0.8.0
circlize 0.4.16
ClusterR 1.3.3
colorRamps 2.3.4
data.table 1.17.8
dplyr 1.1.4
Exact 3.3
ggplot2 3.5.2
ggpubr 0.6.1
gplots 3.2.0
gridExtra 2.3
Hmisc 5.2.3
infercnv 1.20.0
NbClust 3.0.1
plotrix 3.8.4
RColorBrewer 1.1.3
survival 3.8.3
survminer 0.5.0
tidyverse 2.0.0

Additional installation instructions and package requirements are provided in 15_Instructions and Packages.txt.

Repository structure and analysis workflow:

1_Preprocessing of Single-cell RNA-seq Datasets.ipynb - Preprocessing of each single-cell RNA-seq dataset (Input- counts h5ad file; Output: QC processed h5ad file).

2_Converting h5ad files to csv.ipynb - Converting h5ad to CSV (for R **InferCNV** analysis. Input counts h5ad files; Output: csv file).
 
3_R-Infercnv analysis - InferCNV analysis (Input- csv file; Output tumor/immune/stroma h5ad files for each dataset).

4_Immune Integration and clustering.ipynb - Immune integration, clustering, and analysis (includes T cell and myeloid clustering and analysis). Input: all immune h5ad files; Output: Immune clustering, T and Myeloid subclustering, analysis process and graphs.

5_Stroma clustering and analysis.ipynb - Stroma integration, clustering, and analysis (Input: all stroma h5ad files; Output: Stroma clustering, analysis process and graphs).

6_Tumor.ipynb - Tumor integration (Input: all tumor h5ad files; Output: integrated tumor cells dataset).

7_cnmf_scripts.zip - cNMF analysis code (Input: tumor cells counts; Output: cNMF analysis including top 100 associated genes and usage for each program). Intended for execution on a multi-core server; the pipeline runs in parallel across available CPU cores (see README.txt inside the ZIP folder).

8_Integrated dataset analysis.ipynb - Integrated dataset analysis: integrating tumor, immune, and stromal cells into one dataset (Input: h5ad files of the tumor/immune/stroma populations; Output: integrated dataset ; immune abundance and ligand–receptor analysis).

9_Cohort characteristics.ipynb - Cohort characteristics graphs (Input: h5ad integrated dataset file; Output: Cohort characteristics graphs- Figure 1 ).

10_LA tumor analysis.ipynb - Luminal A tumor analysis (Input integrated tumor h5ad file, Output: Luminal A analysis graphs, ER signatures analysis and associated graphs).

11_infercnv_scripts.zip - Containing the **InferCNV** chromosome and gene amplification analysis for integrated samples, separated by breast cancer subtype and stage (Input: integrated csv counts file- see readme.txt inside; Output: txt and DAT files required for the following scripts analysis).

12_Luminal A ER signature-CNV analysis .ipynb - Containing the Luminal A chromosome and gene amplification analysis (Input integrated tumor h5ad file and INFERCNV txt and DAT output files from step 12).

13_TNBC Tumor analysis.ipynb - TNBC tumor analysis: cNMF programs, CNV chromosomal amplifications (Input integrated Tumor h5ad file, csv and excel files of the chromosome and gene amplification analysis from step 12; Output cNMF analysis and associated graphs).

14_Validation_Cohorts folder — Contains individual files for processing and analyzing each validation dataset. Each file specifies its own required input(s) in the code.

15_Instructions and Packages.txt - Necessary Python/R packages and installation information.

16_DEMO_INPUT.h5ad - h5ad DEMO input file for the first (Preprocessing) code (DEMO_1). Can be found on our website: https://breastscape.twix.technion.ac.il/

17_DEMO_after_preprocessing.h5ad - h5ad DEMO output file for the first (Preprocessing) code (DEMO_2). Can be found on our website: https://breastscape.twix.technion.ac.il/

