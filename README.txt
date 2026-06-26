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

16_DEMO_INPUT.h5ad - h5ad DEMO input file for the first (Preprocessing) code. Can be found on our website: https://breastscape.twix.technion.ac.il/

17_DEMO_after_preprocessing.h5ad - h5ad DEMO output file for the first (Preprocessing) code. Can be found on our website: https://breastscape.twix.technion.ac.il/

