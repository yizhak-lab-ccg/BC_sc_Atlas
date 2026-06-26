import scanpy as sc
import os
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt

#Upload the integrated tumor dataset

adata=sc.read_h5ad(r'path\Tumor.h5ad')
adata

#The following code include multiple steps combining the CNV results with the luminal A anndata objects.
#From the INFERCNV R pachage output it extracts for each cell the CNV cluster, then for each cluster the CNV(1-6) of all chromosomal regions.
#Next it matched the chromosome regions with the gene coded in the region.
#Following these steps the code calcualted (using Fisher exact test) the statistical significance of deleted vs amplyfed genes in ER+ enriched vs ER- enriched samples
#P values were then FDR corrected.
#INFERCNV output used:
    #'infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observation_groupings.txt'
    #'17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.DAT'

#C.TYPE=BC Subtype
#P_P=sample
#BIOPSY=Biopsy location
#Create a new column 'BIOPSY_GROUP' to categorize BIOPSY values
adata.obs['BIOPSY_GROUP'] = adata.obs['BIOPSY'].apply(
    lambda x: 'BREAST' if x == 'BREAST' else ('LN' if x == 'LN' else ('MET' if x not in ['BLOOD', 'NA'] else None))
)

#Filter out cells with 'BIOPSY_GROUP' as None
adata= adata[adata.obs['BIOPSY_GROUP'].notnull()]

# Load the CSV file
usage_data = pd.read_csv(r"CNMF Results folder\usage_30.csv")

# Add "Program+" prefix to each column name except 'CELL'
usage_data = usage_data.rename(columns={col: f"Program_{col}" for col in usage_data.columns if col != 'CELL'})

# Set 'CELL' as the index to match with `adata.obs_names`
usage_data.set_index('CELL', inplace=True)

# Join the data with `adata.obs`
adata.obs = adata.obs.join(usage_data, how="left")

# Identify the Program+ columns in adata.obs
program_columns = [col for col in adata.obs.columns if col.startswith("Program")]

# Calculate the Program with the maximum value for each cell
adata.obs['MAX_USAGE'] = adata.obs[program_columns].idxmax(axis=1)

# Print a sample of the new column to confirm
print(adata.obs[['MAX_USAGE']].head())

# Initialize column in adata.obs
adata.obs["CNV_GROUP"] = pd.NA

# Set base path
base_path = r"CNV results folder path"

# Get all subfolders
subfolders = sorted([f.path for f in os.scandir(base_path) if f.is_dir()])

# File name to look for (CNV cluster per cell)
file_name = 'infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observation_groupings.txt'

for i, folder in enumerate(subfolders):
    file_path = os.path.join(folder, file_name)

    if not os.path.exists(file_path):
        print(f"[{i}] Skipping folder (file not found): {folder}")
        continue

    # Read and clean
    rows = []
    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]  # Skip header

        for line in lines:
            parts = line.strip().split(' ')
            if len(parts) >= 2:
                cell_id = parts[0].strip().replace('"', '')
                group = parts[1].strip().replace('"', '')
                rows.append((cell_id, group))

    infercnv_df = pd.DataFrame(rows, columns=["cell_id", "group"])

    # Clean index
    adata.obs_names = adata.obs_names.astype(str).str.strip()
    infercnv_df["cell_id"] = infercnv_df["cell_id"].astype(str).str.strip()

    # Find matching cells
    matched_df = infercnv_df[infercnv_df["cell_id"].isin(adata.obs_names)]
    
    # Attach groups to adata.obs
    for _, row in matched_df.iterrows():
        cell = row["cell_id"]
        group = row["group"]
        if group.lower() not in ['nan', 'na', '']:
            adata.obs.at[cell, "CNV_GROUP"] = group

    # Report
    attached = matched_df["group"].str.lower().isin(["nan", "na", ""]).apply(lambda x: not x).sum()
    print(f"[{i}] Folder: {os.path.basename(folder)}")
    print(f"  Cells in file: {len(infercnv_df)}")
    print(f"  Matched in adata: {len(matched_df)}")
    print(f"  Groups attached (non-NaN): {attached}\n")


# === Setup ===
base_path = r"infercnv folder path"
subfolders = sorted([f.path for f in os.scandir(base_path) if f.is_dir()])
first_folder = subfolders[0]

groupings_file = os.path.join(first_folder, 'infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observation_groupings.txt')
cnv_genes_file = os.path.join(first_folder, '17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.DAT')

# === Step 1: Load CNV groupings ===
rows = []
with open(groupings_file, 'r') as f:
    lines = f.readlines()[1:]
    for line in lines:
        parts = line.strip().split(' ')
        if len(parts) >= 2:
            cell_id = parts[0].strip().replace('"', '')
            group = parts[1].strip().replace('"', '')
            rows.append((cell_id, group))

infercnv_df = pd.DataFrame(rows, columns=["cell_id", "group"])
adata.obs_names = adata.obs_names.astype(str).str.strip()
infercnv_df["cell_id"] = infercnv_df["cell_id"].astype(str).str.strip()

# Attach CNV_GROUP
adata.obs["CNV_GROUP"] = pd.NA
matched_df = infercnv_df[infercnv_df["cell_id"].isin(adata.obs_names)]
for _, row in matched_df.iterrows():
    if row["group"].lower() not in ['nan', 'na', '']:
        adata.obs.at[row["cell_id"], "CNV_GROUP"] = row["group"]

# === Step 2: Load CNV .DAT file and clean group names ===
cnv_genes_df = pd.read_csv(cnv_genes_file, sep='\t')
cnv_genes_df['cell_group_name_cleaned'] = cnv_genes_df['cell_group_name'].apply(lambda x: x.split('.')[-1] if isinstance(x, str) else x)

# === Step 3: Merge metadata (excluding MAX_USAGE) ===
base_columns = ['paper', 'ORIGIN', 'BIOPSY_GROUP', 'STAGE', 'ER', 'PR', 'HER2']
program_columns = [col for col in adata.obs.columns if col.startswith("Program")]
selected_columns = base_columns + program_columns

# Create metadata + CNV gene long format
merged_rows = []

for _, cnv_row in cnv_genes_df.iterrows():
    group_name = cnv_row['cell_group_name_cleaned']
    matching_obs = adata.obs[adata.obs['CNV_GROUP'] == group_name][selected_columns].copy()

    if not matching_obs.empty:
        for cell_id, obs_row in matching_obs.iterrows():
            combined = cnv_row.to_dict()
            combined['cell_id'] = cell_id
            for col in obs_row.index:
                combined[col] = obs_row[col]
            merged_rows.append(combined)

# Final merged table without MAX_USAGE
final_df = pd.DataFrame(merged_rows)

# === Step 4: Compute mean Program values for each CNV group ===
group_program_means = (
    adata.obs
    .groupby("CNV_GROUP")[program_columns]
    .mean()
    .reset_index()
    .rename(columns={"CNV_GROUP": "cell_group_name_cleaned"})
)

# === Step 5: Attach program means to the CNV .DAT file ===
cnv_genes_with_programs = cnv_genes_df.merge(group_program_means, on="cell_group_name_cleaned", how="left")

# === Step 6: Add 'C.TYPE', 'BIOPSY_GROUP', 'paper' to CNV summary ===
representative_cells = (
    adata.obs.reset_index()[["index", "CNV_GROUP", "C.TYPE", "BIOPSY_GROUP", "paper"]]
    .dropna(subset=["CNV_GROUP"])
    .drop_duplicates(subset=["CNV_GROUP"])
    .rename(columns={
        "CNV_GROUP": "cell_group_name_cleaned",
        "C.TYPE": "C_TYPE"
    })
)

# Merge metadata into CNV summary
cnv_genes_with_programs = cnv_genes_with_programs.merge(
    representative_cells,
    on="cell_group_name_cleaned",
    how="left"
)

# === Step 7: Save results ===
final_df.to_csv("merged_cell_level_metadata.csv", index=False)
cnv_genes_with_programs.to_csv("cnv_gene_table_with_program_means.csv", index=False)

# Preview
print("Processing complete.")
print("Merged CNV gene table with program means and metadata:")
print(cnv_genes_with_programs.head())

# Save as CSV instead of Excel (no size limit)
output_csv_path = os.path.join(base_path, "all_cnv_genes_with_programs.csv")
merged_all.to_csv(output_csv_path, index=False)

print("Saved as CSV (no Excel row limit):")
print(output_csv_path)
merged_all

# Load the merged CSV file
merged_all = pd.read_csv(r"folder path\all_cnv_genes_with_programs.csv")

# Filter for C.TYPE == "Luminal A"
luminal_a_only = merged_all[merged_all["C_TYPE"] == "Luminal A"].copy()

# Load the Luminal A subset
df = pd.read_csv(r"folder path\luminal_a_only.csv")

# Get unique genes and unique CNV cell groups
genes = df['gene'].unique()
cell_groups = df['cell_group_name'].unique()

# Initialize output matrix with genes as rows and CNV clusters as columns
matrix = pd.DataFrame(index=genes, columns=cell_groups)

# Fill matrix with 'state' values
for _, row in df.iterrows():
    gene = row['gene']
    group = row['cell_group_name']
    state = row['state']
    matrix.at[gene, group] = state

# Create BIOPSY_GROUP row from mapping
group_biopsy_map = df.drop_duplicates(subset=["cell_group_name"])[["cell_group_name", "BIOPSY_GROUP"]].set_index("cell_group_name")["BIOPSY_GROUP"]
matrix.loc["BIOPSY_GROUP"] = matrix.columns.map(group_biopsy_map)

# Remove columns where all rows are NaN
matrix = matrix.dropna(axis=1, how='all')

# Drop rows that are all NaNs
matrix = matrix.dropna(axis=0, how='all')

# Transpose (CNV group as rows, genes as columns)
matrix = matrix.T

# Save as CSV (not Excel, due to possible large size)
output_path = r"folder path\gene_state_matrix_luminal_a_with_biopsy_group.csv"
matrix.to_csv(output_path)

print("Matrix saved to:")
print(output_path)

# === Load the matrix ===
file_path = r"C:\Users\kerenylab\Documents\breast cancer\INFERCNV\gene_state_matrix_luminal_a_with_biopsy_group.csv"
matrix = pd.read_csv(file_path, index_col=0)

# === Extract BIOPSY_GROUP row safely ===
if 'BIOPSY_GROUP' in matrix.index:
    biopsy_row = matrix.loc['BIOPSY_GROUP']
    matrix = matrix.drop(index='BIOPSY_GROUP')
elif 'BIOPSY_GROUP' in matrix.columns:
    biopsy_row = matrix['BIOPSY_GROUP']
    matrix = matrix.drop(columns='BIOPSY_GROUP')
else:
    raise KeyError("BIOPSY_GROUP not found in index or columns.")

# === Define amplification/deletion state categories ===
deletion_states = [1, 2]
amplification_states = [4, 5, 6]
raw_results = []

# === Run Fisher's exact test on each gene ===
for gene in matrix.columns:
    values = matrix[gene].astype(int)
    biopsies = biopsy_row[values.index]

    is_deleted = values.isin(deletion_states)
    is_amplified = values.isin(amplification_states)
    is_breast = biopsies == "BREAST"
    is_met = biopsies == "MET"

    # --- Deleted in BREAST vs MET ---
    a = (is_deleted & is_breast).sum()
    b = (~is_deleted & is_breast).sum()
    c = (is_deleted & is_met).sum()
    d = (~is_deleted & is_met).sum()
    table = [[a, b], [c, d]]
    if all(sum(row) > 0 for row in table):
        _, pval = fisher_exact(table)
        raw_results.append({
            "gene": gene,
            "comparison": "Deleted_BREAST_vs_MET",
            "p_value": pval,
            "BREAST_deleted": a,
            "MET_deleted": c
        })

    # --- Deleted in MET vs BREAST ---
    table = [[c, d], [a, b]]
    if all(sum(row) > 0 for row in table):
        _, pval = fisher_exact(table)
        raw_results.append({
            "gene": gene,
            "comparison": "Deleted_MET_vs_BREAST",
            "p_value": pval,
            "MET_deleted": c,
            "BREAST_deleted": a
        })

    # --- Amplified in BREAST vs MET ---
    a = (is_amplified & is_breast).sum()
    b = (~is_amplified & is_breast).sum()
    c = (is_amplified & is_met).sum()
    d = (~is_amplified & is_met).sum()
    table = [[a, b], [c, d]]
    if all(sum(row) > 0 for row in table):
        _, pval = fisher_exact(table)
        raw_results.append({
            "gene": gene,
            "comparison": "Amplified_BREAST_vs_MET",
            "p_value": pval,
            "BREAST_amplified": a,
            "MET_amplified": c
        })

    # --- Amplified in MET vs BREAST ---
    table = [[c, d], [a, b]]
    if all(sum(row) > 0 for row in table):
        _, pval = fisher_exact(table)
        raw_results.append({
            "gene": gene,
            "comparison": "Amplified_MET_vs_BREAST",
            "p_value": pval,
            "MET_amplified": c,
            "BREAST_amplified": a
        })

# === Apply FDR correction and save results ===
results_df = pd.DataFrame(raw_results)

# Apply FDR correction
results_df["FDR"] = multipletests(results_df["p_value"], method="fdr_bh")[1]

# Filter significant results
results_sig = results_df[results_df["FDR"] < 0.05].copy()
results_sig.sort_values("FDR", inplace=True)

# Save results to CSV
output_path = r"folder path\Fisher_tests_luminal_a_fdr_sig_only.csv"
results_sig.to_csv(output_path, index=False)

print("Fisher test with FDR complete. Significant results (FDR < 0.05) saved to:")
print(output_path)

# === Load Fisher results with FDR from CSV ===
results_path = r"folder path\Fisher_tests_luminal_a_fdr_sig_only.csv"
results_df = pd.read_csv(results_path)

# === Define gene sets ===
genes_breast_vs_met = [  # your long Luminal A–like list
    "ESR1", "TRPS1", "THSD4", "ERBB4", "ANKRD30A", "AFF3", "PATJ", "PBX1", "LRBA", "TTC6",
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
    "ZCCHC7", "KDM2A", "DLEU2", "TNRC18", "ZNF91", "SETD5", "FBXL7", "NIPBL", "EXOC6B", "DLG1"]

genes_met_vs_breast = [  # your MET-like list
    "PTMA", "FTH1", "TMSB10", "FAU", "TMSB4X", "EEF1A1", "UBA52", "FTL", "NACA",
    "H3F3A", "SERF2", "EIF1", "TPT1", "GAPDH", "MYL6", "EEF1D", "MIF", "B2M",
    "OAZ1", "ATP5F1E", "COX7C", "H3F3B", "HSPB1", "FOS", "COX6C", "CALM2",
    "PFDN5", "RACK1", "SRP14", "PPDPF", "TMA7", "BTF3", "JUN", "CHCHD2", "CFL1",
    "S100A11", "HSP90AB1", "UBB", "NDUFA4", "COX6A1", "TOMM7", "UQCR11", "PFN1",
    "PPIA", "COX6B1", "UBL5", "UQCRQ", "COX4I1", "UBC", "UQCRB", "PSMA7", "OST4",
    "IFITM3", "COX7A2", "JUNB", "HSP90AA1", "DYNLL1", "EEF1B2", "S100A10", "LGALS1",
    "CD63", "TXN", "HINT1", "CRIP1", "PRDX1", "S100A6", "SEC61G", "COX8A", "POLR2L",
    "KRT18", "ELOB", "TMEM258", "HSPA1A", "COMMD6", "COX5B", "SUMO2", "HSPE1", "MYC","MYCBP2"
]

# === Filter Fisher results by gene and comparison ===
amp_breast_vs_met = results_df[
    (results_df["comparison"] == "Amplified_BREAST_vs_MET") &
    (results_df["gene"].isin(genes_breast_vs_met))
]

amp_met_vs_breast = results_df[
    (results_df["comparison"] == "Amplified_MET_vs_BREAST") &
    (results_df["gene"].isin(genes_met_vs_breast))
]

# === Display results ===
print("\n Genes from list significantly amplified in BREAST vs MET:")
print(amp_breast_vs_met[["gene", "p_value", "FDR"]])

print("\n Genes from list significantly amplified in MET vs BREAST:")
print(amp_met_vs_breast[["gene", "p_value", "FDR"]])

# === Save to CSV ===
amp_breast_vs_met.to_csv(
    r"folder path\Amplified_in_BREAST_vs_MET_selected.csv",
    index=False
)

amp_met_vs_breast.to_csv(
    r"folder path\Amplified_in_MET_vs_BREAST_selected.csv",
    index=False
)


# === Load matrix ===
matrix_path = r"folder path\gene_state_matrix_luminal_a_with_biopsy_group.csv"
df = pd.read_csv(matrix_path, index_col=0)

# === Load ESR1-based ER status file ===
er_path = r"folder path\CNV_GROUP_with_ESR1_percent.csv"
er_df = pd.read_csv(er_path)

# If ER_Status column is missing, compute it from mean ESR1%
if "ER_Status" not in er_df.columns:
    def parse_percent_list(pct_string):
        if pd.isna(pct_string): return []
        return [float(x.replace('%', '')) for x in pct_string.split(',') if x.strip()]
    
    er_df["ESR1_List"] = er_df["ESR1_Positive_Percent"].apply(parse_percent_list)
    er_df["ESR1_Mean_Percent"] = er_df["ESR1_List"].apply(lambda x: sum(x)/len(x) if x else 0)
    er_df["ER_Status"] = er_df["ESR1_Mean_Percent"].apply(lambda x: "POSITIVE" if x >= 50 else "NEGATIVE")

# === Build a mapping from CNV group name → ER status ===
cnv_to_er = dict(zip(er_df["CNV_GROUP"], er_df["ER_Status"]))

# === Extract CNV_GROUP from row names (index) ===
# Your rownames are like: "BREAST_TUMOR.BREAST_TUMOR_s3" → match "BREAST_TUMOR_s3"
df["CNV_GROUP"] = df.index.to_series().apply(lambda x: x.split('.')[-1])

# === Map ER_Status===
df["ER_Status"] = df["CNV_GROUP"].map(cnv_to_er).fillna("Unknown")

# === Move ER_Status next to BIOPSY_GROUP (for clarity)
cols = df.columns.tolist()
if "BIOPSY_GROUP" in cols:
    cols.insert(cols.index("BIOPSY_GROUP") + 1, cols.pop(cols.index("ER_Status")))
    df = df[cols]

# === Save updated file ===
output_path = r"folder path\gene_state_matrix_luminal_a_with_ER_column.csv"
df.to_csv(output_path)

print("ER_Status added and saved to:")
print(output_path)


# === Load matrix with ER column ===
file_path = r"folder path\gene_state_matrix_luminal_a_with_ER_column.csv"
matrix = pd.read_csv(file_path, index_col=0)

# === Extract ER_Status and drop from matrix ===
er_status = matrix["ER_Status"]
matrix = matrix.drop(columns=["ER_Status"])

# === Drop any non-numeric metadata columns ===
non_gene_cols = ["BIOPSY_GROUP", "C.TYPE", "CNV_GROUP"]
matrix = matrix.drop(columns=[col for col in non_gene_cols if col in matrix.columns])

# === Ensure all values are numeric and fill NaNs ===
matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(3)

# === Define state categories ===
deletion_states = [1, 2]
amplification_states = [4, 5, 6]
raw_results = []

# === Fisher's exact test ===
for gene in matrix.columns:
    values = matrix[gene].astype(int)
    status = er_status[values.index]

    is_deleted = values.isin(deletion_states)
    is_amplified = values.isin(amplification_states)
    is_er_pos = status == "POSITIVE"
    is_er_neg = status == "NEGATIVE"

    # --- Deleted in ER+ vs ER− ---
    a = (is_deleted & is_er_pos).sum()
    b = (~is_deleted & is_er_pos).sum()
    c = (is_deleted & is_er_neg).sum()
    d = (~is_deleted & is_er_neg).sum()
    table = [[a, b], [c, d]]
    if all(sum(row) > 0 for row in table):
        _, pval = fisher_exact(table)
        raw_results.append({
            "gene": gene,
            "comparison": "Deleted_ERpos_vs_ERneg",
            "p_value": pval,
            "ERpos_deleted": a,
            "ERneg_deleted": c
        })

    # --- Amplified in ER+ vs ER− ---
    a = (is_amplified & is_er_pos).sum()
    b = (~is_amplified & is_er_pos).sum()
    c = (is_amplified & is_er_neg).sum()
    d = (~is_amplified & is_er_neg).sum()
    table = [[a, b], [c, d]]
    if all(sum(row) > 0 for row in table):
        _, pval = fisher_exact(table)
        raw_results.append({
            "gene": gene,
            "comparison": "Amplified_ERpos_vs_ERneg",
            "p_value": pval,
            "ERpos_amplified": a,
            "ERneg_amplified": c
        })

# === FDR correction and save ===
results_df = pd.DataFrame(raw_results)
results_df["FDR"] = multipletests(results_df["p_value"], method="fdr_bh")[1]
results_sig = results_df[results_df["FDR"] < 0.05].copy()
results_sig = results_sig.sort_values("FDR")

# === Save ===
output_path = r"folder path\Fisher_tests_ERpos_vs_ERneg_fdr_sig_only.csv"
results_sig.to_csv(output_path, index=False)

print("Fisher test (ER+ vs ER−) complete. FDR < 0.05 results saved to:")
print(output_path)

# === Load FDR-corrected Fisher test results ===
fisher_file = r"folder path\Fisher_tests_ERpos_vs_ERneg_fdr_sig_only.csv"
fisher_df = pd.read_csv(fisher_file)

# === Define gene lists ===
er_positive_genes = [
    "ESR1", "TRPS1", "THSD4", "ERBB4", "ANKRD30A", "AFF3", "PATJ", "PBX1", "LRBA", "TTC6",
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
    "ZCCHC7", "KDM2A", "DLEU2", "TNRC18", "ZNF91", "SETD5", "FBXL7", "NIPBL", "EXOC6B", "DLG1",
    "BTBD9", "GMDS-DT", "ZNF638", "FNDC3B", "TANC2", "KANSL1", "EXT1", "RALGPS2", "BAZ2B"
]

er_negative_genes = [
    "PTMA", "FTH1", "TMSB10", "FAU", "TMSB4X", "EEF1A1", "UBA52", "FTL", "NACA",
    "H3F3A", "SERF2", "EIF1", "TPT1", "GAPDH", "MYL6", "EEF1D", "MIF", "B2M",
    "OAZ1", "ATP5F1E", "COX7C", "H3F3B", "HSPB1", "FOS", "COX6C", "CALM2",
    "PFDN5", "RACK1", "SRP14", "PPDPF", "TMA7", "BTF3", "JUN", "CHCHD2", "CFL1",
    "S100A11", "HSP90AB1", "UBB", "NDUFA4", "COX6A1", "TOMM7", "UQCR11", "PFN1",
    "PPIA", "COX6B1", "UBL5", "UQCRQ", "COX4I1", "UBC", "UQCRB", "PSMA7", "OST4",
    "IFITM3", "COX7A2", "JUNB", "HSP90AA1", "DYNLL1", "EEF1B2", "S100A10", "LGALS1",
    "CD63", "TXN", "HINT1", "CRIP1", "PRDX1", "S100A6", "SEC61G", "COX8A", "POLR2L",
    "KRT18", "ELOB", "TMEM258", "HSPA1A", "COMMD6", "COX5B", "SUMO2", "HSPE1", "MYC", "MYCBP2"
]

# === Identify amplified genes in ER+ and ER− ===
amplified_er_pos = fisher_df[
    (fisher_df["comparison"] == "Amplified_ERpos_vs_ERneg") &
    (fisher_df["gene"].isin(er_positive_genes))
]["gene"].unique().tolist()

amplified_er_neg = fisher_df[
    (fisher_df["comparison"] == "Amplified_ERneg_vs_ERpos") &
    (fisher_df["gene"].isin(er_negative_genes))
]["gene"].unique().tolist()

# === Print the results ===
print("Genes amplified in ER-positive group:")
print(amplified_er_pos)

print("Genes amplified in ER-negative group:")
print(amplified_er_neg)


##Final heatmap code##
# === Load expression matrix and ER status ===
matrix_path = r"folder path\gene_state_matrix_luminal_a_with_ER_column.csv"
matrix = pd.read_csv(matrix_path, index_col=0)

# Separate ER status
er_status = matrix["ER_Status"]
matrix = matrix.drop(columns=["ER_Status"], errors="ignore")
matrix = matrix.drop(columns=[col for col in ["BIOPSY_GROUP", "CNV_GROUP", "C.TYPE"] if col in matrix.columns], errors="ignore")
matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(3)

# === Load Fisher results and define gene list ===
fdr_path = r"folder path\Fisher_tests_ERpos_vs_ERneg_fdr_sig_only.csv"
fdr_df = pd.read_csv(fdr_path)

# Gene list with exclusions
amplified_er_pos_genes = [
    'SIK3', 'KMT2C', 'GALNT7', 'VPS13B', 'DNAH14', 'IMMP2L', 'UBR5', 'CUX1', 'SRPK2',
    'TRPS1','SAMD12','EXT1','TRAPPC9','PTK2','DENND1B','LPP','NSMCE2','SMYD3','ASH1L',
    'SH3RF1', 'BIRC6', 'DACH1', 'ANKRD17', 'FARP1', 'TULP4', 'BAZ2B', 'PKP4', 'IGF1R', 'LRBA',
    'GSE1', 'PDXDC1', 'ESR1', 'EVL', 'MBD5', 'ARHGAP32', 'ZCCHC7', 'RABEP1', 'ZMYND8', 'EXOC4',
    'ARID1B', 'NEBL', 'C1orf21', 'DLG1', 'AFF3'
]
genes_to_remove = {'EVL','ZMYND8','DLG1','C1orf21','PDXDC1','TRPS1', 'SAMD12', 'EXT1', 'TRAPPC9', 'PTK2', 'DENND1B', 'LPP', 'NSMCE2', 'SMYD3', 'ASH1L', 'UBR5', 'VPS13B', 'DNAH14'}
filtered_gene_list = [g for g in amplified_er_pos_genes if g not in genes_to_remove]

# === Order genes by mean expression ===
valid_genes = [g for g in filtered_gene_list if g in matrix.columns]
gene_means = matrix[valid_genes].mean(axis=0).sort_values(ascending=False)
sorted_genes = gene_means.index.tolist()

# === Subset matrix and sort samples ===
subset_matrix = matrix[sorted_genes].T
er_pos = er_status[er_status == "POSITIVE"].index.tolist()
er_neg = er_status[er_status == "NEGATIVE"].index.tolist()
ordered_samples = er_pos + er_neg
subset_matrix = subset_matrix[ordered_samples]

# ER color bar
er_color_map = {'POSITIVE': '#4575b4', 'NEGATIVE': '#d73027'}
col_colors = er_status.loc[ordered_samples].map(er_color_map)

# === Plot ===
plt.figure(figsize=(8.27 / 1.2, 11.69 / 2))  # PDF size ≈ (6.89, 9.74)
sns.heatmap(
    subset_matrix,
    cmap=sns.diverging_palette(240, 10, as_cmap=True),
    linewidths=0.3,
    linecolor="gray",
    cbar_kws={"label": "Gene State (1–6)"},
    xticklabels=False,
    yticklabels=True,
    vmin=1,
    vmax=6
)

# Add stripe for ER status
ax = plt.gca()
for idx, label in enumerate(ordered_samples):
    ax.add_patch(plt.Rectangle(
        (idx, -0.5), 1, 0.3,
        color=col_colors[label],
        transform=ax.transData,
        clip_on=False
    ))

# Add vertical black line to separate ER+ from ER−
ax.axvline(x=len(er_pos), color='black', linewidth=2)

# === Format ===
plt.title("ER+ Amplified Genes", fontsize=9)
plt.xlabel("CNV_Clusters", fontsize=8)
plt.ylabel("Genes", fontsize=8)
plt.yticks(fontsize=8)

# === Save ===
output_path = r"folder path\heatmap_ERpos_genes_sorted_by_mean.pdf"
plt.tight_layout()
plt.savefig(output_path, format="pdf", dpi=300)

print(f"Heatmap saved to:\n{output_path}")
plt.show()

