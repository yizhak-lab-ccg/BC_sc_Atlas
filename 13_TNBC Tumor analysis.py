import scanpy as sc
import os
import pandas as pd
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from matplotlib.lines import Line2D

#Upload integrated tumor dataset
adata=sc.read_h5ad(r'path\Tumor.h5ad')
adata

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

#C.TYPE=BC Subtype
#P_P=sample
#BIOPSY=Biopsy location

#Create a new column 'BIOPSY_GROUP' to categorize BIOPSY values
adata.obs['BIOPSY_GROUP'] = adata.obs['BIOPSY'].apply(
    lambda x: 'BREAST' if x == 'BREAST' else ('LN' if x == 'LN' else ('MET' if x not in ['BLOOD', 'NA'] else None))
)

# Filter out cells with 'BIOPSY_GROUP' as None
adata= adata[adata.obs['BIOPSY_GROUP'].notnull()]

# Subset the AnnData object to TNBC cells
adata_tnbc = adata[adata.obs['C.TYPE'] == 'TNBC']
adata_tnbc

sc.pp.highly_variable_genes(adata_tnbc, n_top_genes=2000, batch_key="paper")

sc.pl.highly_variable_genes(adata_tnbc)

adata_tnbc.raw = adata_tnbc

adadata_tnbc = adata_tnbc[:, adata_tnbc.var.highly_variable]

sc.tl.pca(adata_tnbc)
sc.pp.neighbors(adata_tnbc)
sc.tl.umap(adata_tnbc)
sc.pl.umap(adata_tnbc)

# Remove cells with NaN in Program_21
adata_tnbc = adata_tnbc[~adata_tnbc.obs["Program_21"].isna()].copy()

# Plot UMAP colored by Program_21 expression
sc.pl.umap(
    adata_tnbc,
    color='Program_21',
    cmap='viridis',
    size=20,
    title='UMAP colored by Program_21 expression',
    frameon=False
)


# Remove cells with NaN in Program_4
adata_tnbc = adata_tnbc[~adata_tnbc.obs["Program_4"].isna()].copy()

# Plot UMAP colored by Program_21 expression
sc.pl.umap(
    adata_tnbc,
    color='Program_4',
    cmap='viridis',
    size=20,
    title='UMAP colored by Program_4 expression',
    frameon=False
)


# === Handle NaNs in categorical 'P_P' ===
if 'P_P' in adata.obs.columns:
    if pd.api.types.is_categorical_dtype(adata.obs['P_P']):
        if 'P9' not in adata.obs['P_P'].cat.categories:
            adata.obs['P_P'] = adata.obs['P_P'].cat.add_categories(['P9'])
    adata.obs['P_P'] = adata.obs['P_P'].fillna('P9')

# Group and calculate counts
program_paper_counts = adata.obs.groupby(['MAX_USAGE', 'P_P']).size().unstack(fill_value=0)

# Convert to percentages per MAX_USAGE row
program_paper_percentages = program_paper_counts.div(program_paper_counts.sum(axis=1), axis=0) * 100

# Subset only 'Program_4' and 'Program_21'
program_paper_percentages = program_paper_percentages.loc[['Program_4', 'Program_21']]

# Only include relevant P_Ps
used_pps = program_paper_percentages.columns[(program_paper_percentages != 0).any()].tolist()
df = program_paper_percentages[used_pps]

# Setup colors
cmap = plt.get_cmap("tab20")
colors = [cmap(i % 20) for i in range(len(used_pps))]
color_map = dict(zip(used_pps, colors))

# Create fixed-size figure
fig = plt.figure(figsize=(8.27 / 3.5, 11.69 / 4.2))  # ≈ 2.36 x 2.78 in

# 70% plot, 30% legend
ax = fig.add_axes([0.05, 0.1, 0.65, 0.85])  # [left, bottom, width, height]

# Bar positions
x = np.arange(len(df))  # [0, 1]
bar_width = 0.6

# Stack bars
bottoms = np.zeros(len(df))
for pp in used_pps:
    heights = df[pp].values
    ax.bar(
        x,
        heights,
        bottom=bottoms,
        width=bar_width,
        label=pp,
        color=color_map[pp],
        edgecolor='black'
    )
    bottoms += heights

# Axis formatting
ax.set_xticks(x)
ax.set_xticklabels(df.index, fontsize=8)
ax.set_ylabel("Percentage (%)", fontsize=8)
ax.set_title("Percentage Distribution of Samples", fontsize=9)
ax.tick_params(axis='y', labelsize=8)
ax.grid(False)

#Legend as DOTS in fixed space
legend_ax = fig.add_axes([0.72, 0.1, 0.25, 0.85])  # [left, bottom, width, height]
legend_ax.axis('off')

legend_handles = [
    Line2D([0], [0], marker='o', color='none', markerfacecolor=color_map[p], markersize=5)
    for p in used_pps
]
legend_ax.legend(
    legend_handles,
    used_pps,
    title="Sample",
    fontsize=8,
    title_fontsize=8,
    ncol=2,
    loc='upper left',
    frameon=False
)

# Save figure as exact size PDF
output_path = r"folder path\stacked_bar_program_4_21.pdf"
fig.savefig(output_path)
plt.show()
plt.close()


###program4 vs program 21 DEG###

#Subset to only Program_4 and Program_21 cells
subset_adata = adata_tnbc[adata_tnbc.obs['MAX_USAGE'].isin(['Program_4', 'Program_21'])].copy()

# Create new category for comparison
subset_adata.obs['Program_Group'] = subset_adata.obs['MAX_USAGE'].astype(str)

# Run DEG analysis
sc.tl.rank_genes_groups(
    subset_adata,
    groupby='Program_Group',
    # Compare Program_4
    # Against Program_21
    method='wilcoxon',
    use_raw=False
)


sc.tl.dendrogram(subset_adata, groupby="Program_Group")

sc.pl.rank_genes_groups_dotplot(subset_adata, groupby="Program_Group",standard_scale="var", n_genes=20)

# Define Program_4 genes
program4_genes = [
    'FOS', 'JUN', 'GADD45B', 'EGR1', 'DUSP1', 'NFKBIA', 'FOSB', 'HSPA1A', 'BTG2', 'ID2'
]

# Define Program_21 hypoxia genes
hypoxia_genes_program21 = [
    'LCN2', 'S100A8', 'S100A9', 'NQO1', 'PGK1', 'CD47', 'SLC16A3', 'ANXA2', 'TIMP1', 'CXCL17'
]

# Combine
genes_to_plot = program4_genes + hypoxia_genes_program21
genes_to_plot = [gene for gene in genes_to_plot if gene in subset_adata.var_names]

# Subset adata
subset_adata = adata_tnbc[adata_tnbc.obs['MAX_USAGE'].isin(['Program_4', 'Program_21'])].copy()
subset_adata.obs['Program_Group'] = subset_adata.obs['MAX_USAGE']

# Generate the dotplot
dotplot = sc.pl.dotplot(
    subset_adata,
    var_names=genes_to_plot,
    groupby='Program_Group',
    standard_scale='var',
    swap_axes=True,
    figsize=(10, 8),
    dot_min=0.05,
    dot_max=0.8,
    show=False)

# Modify x-ticks and font sizes manually
fig = plt.gcf()  # get the current figure

# Rotate x-axis labels (horizontal)
for ax in fig.axes:
    for label in ax.get_xticklabels():
        label.set_rotation(0)  # 0 degrees = horizontal
        label.set_fontsize(8)  # x-axis labels font size

    for label in ax.get_yticklabels():
        label.set_fontsize(8)  # y-axis labels font size

    if ax.get_title():
        ax.set_title(ax.get_title(), fontsize=9)  # Title font size 9

# Adjust colorbar and legend font sizes if they exist
for text in fig.findobj(match=plt.Text):
    if text.get_text() in ['Mean expression\nin group', 'Fraction of cells\nin group (%)']:
        text.set_fontsize(8)

# Step 3: Show the figure
plt.show()

# Extract DEG results for Program_4
deg_program4 = sc.get.rank_genes_groups_df(subset_adata, group='Program_4')

#Extract DEG results for Program_21
deg_program21 = sc.get.rank_genes_groups_df(subset_adata, group='Program_21')

# Filter by adjusted p-value < 0.05 and log fold change > 0.5
filtered_program4 = deg_program4[
    (deg_program4['pvals_adj'] < 0.05) & 
    (deg_program4['logfoldchanges'] > 0.5)
]

filtered_program21 = deg_program21[
    (deg_program21['pvals_adj'] < 0.05) & 
    (deg_program21['logfoldchanges'] > 0.5)
]

# Save to Excel with two sheets
with pd.ExcelWriter('Filtered_DEGs_Program4_vs_Program21.xlsx') as writer:
    filtered_program4.to_excel(writer, sheet_name='Program_4_DEGs', index=False)
    filtered_program21.to_excel(writer, sheet_name='Program_21_DEGs', index=False)

print('Excel file saved: Filtered_DEGs_Program4_vs_Program21.xlsx')


###TP abundance###

#Prepare data
program21_df = adata_tnbc.obs[['P_P', 'BIOPSY_GROUP']].copy()
program21_df['Program_21'] = adata_tnbc.obs['Program_21']

grouped_means_21 = (
    program21_df
    .groupby(['P_P', 'BIOPSY_GROUP'])['Program_21']
    .mean()
    .reset_index()
    .dropna()
)

# Mann–Whitney U test after the modification
group_breast_21 = grouped_means_21[grouped_means_21['BIOPSY_GROUP'] == 'BREAST']['Program_21']
group_met_21 = grouped_means_21[grouped_means_21['BIOPSY_GROUP'] == 'MET']['Program_21']
pval_21 = mannwhitneyu(group_breast_21, group_met_21, alternative='two-sided').pvalue

# Plot
plt.figure(figsize=(7, 5))
sns.boxplot(data=grouped_means_21, x='BIOPSY_GROUP', y='Program_21', palette='pastel')
sns.stripplot(data=grouped_means_21, x='BIOPSY_GROUP', y='Program_21', color='black', size=4, jitter=True)

# Add p-value
plt.text(0.5, grouped_means_21['Program_21'].max(), f"p = {pval_21:.2e}", ha='center', fontsize=9)

plt.xlabel("BIOPSY_GROUP")
plt.ylabel("Mean Program 21 Usage per P_P")
plt.title("Program 21 Usage by Biopsy Group")
plt.tight_layout()
plt.show()


# Prepare data
program4_df = adata_tnbc.obs[['P_P', 'BIOPSY_GROUP']].copy()
program4_df['Program_4'] = adata_tnbc.obs['Program_4']

grouped_means_4 = (
    program4_df
    .groupby(['P_P', 'BIOPSY_GROUP'])['Program_4']
    .mean()
    .reset_index()
    .dropna()
)

# Mann–Whitney U test for BREAST vs MET
group_breast_4 = grouped_means_4[grouped_means_4['BIOPSY_GROUP'] == 'BREAST']['Program_4']
group_met_4 = grouped_means_4[grouped_means_4['BIOPSY_GROUP'] == 'MET']['Program_4']
pval_4 = mannwhitneyu(group_breast_4, group_met_4, alternative='two-sided').pvalue

# Plot
plt.figure(figsize=(7, 5))
sns.boxplot(data=grouped_means_4, x='BIOPSY_GROUP', y='Program_4', palette='pastel')
sns.stripplot(data=grouped_means_4, x='BIOPSY_GROUP', y='Program_4', color='black', size=4, jitter=True)

# Add p-value
plt.text(0.5, grouped_means_4['Program_4'].max() * 1.05, f"p = {pval_4:.2e}", ha='center', fontsize=9)

plt.xlabel("BIOPSY_GROUP")
plt.ylabel("Mean Program 4 Usage per P_P")
plt.title("Program 4 Usage by Biopsy Group")
plt.tight_layout()
plt.show()

###TROP2 EXPRESSION###

# Extract expression
gene = 'TACSTD2'
expr = adata_tnbc[:, gene].X
expr = expr.toarray().flatten() if hasattr(expr, 'toarray') else expr.flatten()

# Add expression and group labels
adata_tnbc.obs['TACSTD2_expr'] = expr
adata_tnbc.obs['program_group'] = adata_tnbc.obs['MAX_USAGE'].isin(['Program_4', 'Program_21'])
adata_tnbc.obs['program_group'] = adata_tnbc.obs['program_group'].map({True: 'Program_4_or_21', False: 'Other'})
adata_tnbc.obs['sub_program'] = adata_tnbc.obs['MAX_USAGE'].where(
    adata_tnbc.obs['MAX_USAGE'].isin(['Program_4', 'Program_21']),
    other=np.nan
)

# Calculate mean expression per P_P
df = adata_tnbc.obs[['P_P', 'program_group', 'TACSTD2_expr', 'sub_program']].copy()
grouped = (
    df.groupby(['P_P', 'program_group'], observed=True)
    .agg(
        mean_expression=('TACSTD2_expr', 'mean'),
        sub_program=('sub_program', 'first')
    )
    .reset_index()
)

# Remove top 2 highest values in "Other"
other_mask = grouped['program_group'] == 'Other'
if other_mask.sum() >= 2:
    top2_other_indices = grouped[other_mask].nlargest(2, 'mean_expression').index
    grouped = grouped.drop(index=top2_other_indices)

# Mann–Whitney U test
group1 = grouped[grouped['program_group'] == 'Program_4_or_21']['mean_expression']
group2 = grouped[grouped['program_group'] == 'Other']['mean_expression']
pval = mannwhitneyu(group1, group2, alternative='two-sided').pvalue if len(group1) > 0 and len(group2) > 0 else np.nan

# Format p-value
def format_pval(p):
    if pd.isna(p):
        return "n/a"
    if p == 0:
        return "p < 1×10⁻³⁰⁰"
    elif p >= 0.0001:
        return f"p = {round(p, 3)}"
    else:
        base, exp = f"{p:.1e}".split("e")
        exp_clean = exp.lstrip("0").replace("-0", "-")
        superscript = exp_clean.translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
        return f"p = {float(base):.1f}×10{superscript}"

# Plot
plt.figure(figsize=(6, 4))

order = ['Other', 'Program_4_or_21']

# Define dot colors for better contrast
dot_palette = {
    'Program_4': '#1f78b4',   # Dark Blue
    'Program_21': '#e31a1c'   # Dark Red
}

# Boxplot
sns.boxplot(
    data=grouped,
    x='program_group',
    y='mean_expression',
    order=order,
    palette='pastel'
)

# Stripplot for "Other" (black)
sns.stripplot(
    data=grouped[grouped['program_group'] == 'Other'],
    x='program_group',
    y='mean_expression',
    color='black',
    size=3,
    jitter=True,
    order=order
)

# Stripplot for "Program_4_or_21" (colored by program)
sns.stripplot(
    data=grouped[grouped['program_group'] == 'Program_4_or_21'],
    x='program_group',
    y='mean_expression',
    hue='sub_program',
    palette=dot_palette,
    size=3,
    jitter=True,
    dodge=False,
    order=order
)

# Add p-value
plt.text(
    0.5,
    grouped['mean_expression'].max(),
    format_pval(pval),
    ha='center',
    fontsize=8
)

# Formatting
plt.title("Mean TACSTD2 (TROP2) Sample Expression", fontsize=9)
plt.xlabel("Group", fontsize=8)
plt.ylabel("Mean Expression", fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.legend(title='Program', title_fontsize=8, fontsize=8, loc='upper right')
plt.tight_layout()
plt.show()


###CNV ANALYSIS OF TP4 TP21 GENES###
##This code calcualtedd the amplified genes in tumor program 21, comparing metastatic versus primary CNV clusters. 

# Filter for C.TYPE == "TNBC"
#merged_all=saved in this file (following step 12)-"path\all_cnv_genes_with_programs.csv".
tnbc_only = merged_all[merged_all["C_TYPE"] == "TNBC"].copy()

# Save to new Excel file
output_path = r"INFERCNV folder path\tnbc_only.xlsx"
tnbc_only.to_excel(output_path, index=False)

print("TNBC subset saved to:")
print(output_path)

# Start with your existing TNBC dataframe
df = tnbc_only.copy()

# Get unique genes and unique CNV cell groups
genes = df['gene'].unique()
cell_groups = df['cell_group_name'].unique()

# Initialize output matrix with genes as rows and cell_group_name as columns
matrix = pd.DataFrame(index=genes, columns=cell_groups)

# Fill the matrix with 'state' values
for _, row in df.iterrows():
    gene = row['gene']
    group = row['cell_group_name']
    state = row['state']
    matrix.at[gene, group] = state

# Create the BIOPSY_GROUP row (based on cell_group_name to biopsy info)
# We'll map the biopsy group for each unique cell_group_name
group_biopsy_map = df.drop_duplicates(subset=["cell_group_name"])[["cell_group_name", "BIOPSY_GROUP"]].set_index("cell_group_name")["BIOPSY_GROUP"]

# Add it as a new row to the matrix
matrix.loc["BIOPSY_GROUP"] = matrix.columns.map(group_biopsy_map)

# Save or return
matrix.to_excel("gene_state_matrix_with_biopsy_group.xlsx")
#Remove columns where all rows are NaN
matrix = matrix.dropna(axis=1, how='all')
#Drop rows that are all NaNs
matrix = matrix.dropna(axis=0, how='all')
matrix = matrix.T
matrix
print("Matrix saved to gene_state_matrix_with_biopsy_group.xlsx")

# Define the path to your Excel file
file_path = r"path\gene_state_matrix_with_biopsy_group.xlsx"

# Load the Excel file into a DataFrame
matrix_df = pd.read_excel(file_path, index_col=0)

# Preview the data
print("File loaded. Preview:")
print(matrix_df)


###Fisher analysis###

# === Load the matrix ===
file_path = r"path\gene_state_matrix_with_biopsy_group.xlsx"
matrix = pd.read_excel(file_path, index_col=0)

# === Separate BIOPSY_GROUP row ===
biopsy_row = matrix.loc['BIOPSY_GROUP']
matrix = matrix.drop(index='BIOPSY_GROUP')

# === Define states ===
deletion_states = [1, 2]
amplification_states = [4, 5, 6]
results = []

# === Loop through genes ===
for gene in matrix.index:
    values = matrix.loc[gene].astype(int)
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
        results.append({
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
        results.append({
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
        results.append({
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
        results.append({
            "gene": gene,
            "comparison": "Amplified_MET_vs_BREAST",
            "p_value": pval,
            "MET_amplified": c,
            "BREAST_amplified": a
        })

# === Apply FDR correction ===
results_df = pd.DataFrame(results)

# Apply Benjamini–Hochberg FDR
if not results_df.empty:
    results_df["FDR"] = multipletests(results_df["p_value"], method="fdr_bh")[1]
    results_df = results_df[results_df["FDR"] < 0.05]  # Optional: keep only FDR < 0.05
    results_df.sort_values("FDR", inplace=True)

# === Save results ===
output_path = r"path\Fisher_tests_sig_only_with_FDR.xlsx"
results_df.to_excel(output_path, index=False)

print("Fisher test complete. Significant results (FDR < 0.05) saved to:")
print(output_path)


###relavent genes###

# Load Fisher result file
fisher_path = r"path\Fisher_tests_sig_only.xlsx"
fisher_df = pd.read_excel(fisher_path)

# Define your gene list
gene_list = [
    "S100A9", "SOX9", "HSPH1", "CD81", "TMEM45A", "RARRES1", "SLPI", "LCN2", "HSPA6", "CP", "S100A6",
    "S100A8", "CHI3L1", "BAG3", "IER5", "LAMB3", "PHLDA2", "IL18", "NME2", "GPRC5A", "PMEPA1", "HSPA1A",
    "GABRE", "CRYAB", "MGP", "HSP90AA1", "S100A11", "ANXA1", "PTN", "HSPD1", "CD59", "ABL2", "EGR4",
    "S100A10", "MUC20", "S100A4", "ARPC1B", "ITGB8", "STEAP4", "TPT1", "DNAJB1", "RDH10", "IDS", "HSPE1",
    "CSTB", "KYNU", "S100P", "CLDN1", "HSPA1B", "CDK2AP1", "KLF5", "C15orf48", "LAMP2", "DCBLD1", "EGR3",
    "UBC", "PRRG4", "ENC1", "FN1", "MAF", "SDC4", "LATS2", "ZFAND2A", "SCPEP1", "BRI3", "TKT", "CHI3L2",
    "BCAP31", "TNFRSF12A", "SNCG", "TGFBI", "ANXA2", "TGM2", "LTF", "SFN", "CAST", "TNNC2", "CLIC3",
    "CMTM6", "ATF3", "CAMK2N1", "TMSB10", "RCAN1", "DEFB1", "SERPINA1", "DSC2", "CD63", "MAFK", "RPL17",
    "SSPN", "SLC7A5", "CLTB", "P4HB", "SDCBP2", "NQO1", "ITM2B", "CARD16", "MGST1", "SLC25A37", "SSR4"
]

#Filter Fisher results for amplified in MET vs BREAST
amplified_met_genes = fisher_df[
    (fisher_df["comparison"] == "Amplified_MET_vs_BREAST") &
    (fisher_df["gene"].isin(gene_list))
]

# Save or display results
output_path = r"path\Amplified_in_MET_genes_filtered.xlsx"
amplified_met_genes.to_excel(output_path, index=False)

print("Filtered list of amplified genes in MET vs BREAST saved to:")
print(output_path)

##plot heatmap of the results
# === Load matrix ===
file_path = r"path\gene_state_matrix_with_biopsy_group.xlsx"
matrix = pd.read_excel(file_path, index_col=0)

# === Extract and remove BIOPSY_GROUP row ===
biopsy_group = matrix.loc["BIOPSY_GROUP"]
matrix = matrix.drop(index="BIOPSY_GROUP")

# === List of significantly amplified genes ===
amplified_genes = [
    "FN1", "STEAP4", "ANXA1", "PTN", "HSPD1", "HSPE1", "PRRG4", "CD59",
    "MGST1", "ZFAND2A", "MAFK", "TGFBI", "NME2", "PHLDA2", "DSC2", "CD81"
]

# === Subset to amplified genes and BREAST + MET samples ===
selected_columns = biopsy_group[biopsy_group.isin(["MET", "BREAST"])].index
subset = matrix.loc[amplified_genes, selected_columns].fillna(3).astype(float)

# === Generate compact sample names ===
new_labels = []
b_count = 1
m_count = 1
label_to_color = {}
for col in selected_columns:
    if biopsy_group[col] == "BREAST":
        label = f"B_G{b_count}"
        b_count += 1
        color = "#4575b4"
    else:
        label = f"M_G{m_count}"
        m_count += 1
        color = "#d73027"
    new_labels.append(label)
    label_to_color[label] = color

# Rename columns in matrix
subset.columns = new_labels

# === Find index where MET starts for black separator line ===
met_start_idx = next(i for i, label in enumerate(new_labels) if label.startswith("M"))

# === Plot ===
plt.figure(figsize=(14, 8))
ax = sns.heatmap(
    subset,
    cmap=sns.diverging_palette(240, 10, as_cmap=True),
    linewidths=0.3,
    linecolor="gray",
    cbar_kws={"label": "Gene State (1–6)", "orientation": "horizontal"},
    xticklabels=True,
    yticklabels=True,
    vmin=1,
    vmax=6
)

# === Draw black vertical separator between BREAST and MET ===
ax.axvline(x=met_start_idx, color='black', linewidth=3)

# === Add biopsy group color stripe ===
for idx, label in enumerate(new_labels):
    ax.add_patch(plt.Rectangle(
        (idx, -0.5), 1, 0.2,
        color=label_to_color[label],
        transform=ax.transData,
        clip_on=False
    ))

# === Format plot ===
plt.title("TP21 Amplified Genes Across MET vs.Breast Samples", fontsize=14, pad=20)
plt.xticks(rotation=90, fontsize=8)
plt.yticks(fontsize=10)
plt.xlabel("Samples", fontsize=11)
plt.ylabel("Genes", fontsize=11)

# === Legend for sample type ===
legend_elements = [
    Patch(facecolor="#4575b4", edgecolor='black', label='BREAST'),
    Patch(facecolor="#d73027", edgecolor='black', label='MET')
]
plt.legend(handles=legend_elements, title='BIOPSY_GROUP',
           bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)

plt.tight_layout()
plt.show()


### PD1 and LAG3 Correlation with TP ###

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# === Paths ===
excel_path_pd1lag3 = r"path to exal file contaning Mean_PD1_LAG3_combined_per_sample.xlsx"

# === Mean Program_4 per P_P (drop missing P_P) ===
obs = adata_tnbc.obs[['P_P', 'Program_4']].copy()
obs = obs.dropna(subset=['P_P'])

mean_prog4 = (
    obs.groupby('P_P', as_index=False)['Program_4']
    .mean()
    .rename(columns={'Program_4': 'Mean_Program_4'})
)

# === Load Excel and clean (drop missing P_P if exist) ===
pd1lag3_df = pd.read_excel(excel_path_pd1lag3)
pd1lag3_df.columns = pd1lag3_df.columns.str.strip()
pd1lag3_df = pd1lag3_df.dropna(subset=['P_P'])

# === Merge and drop NA numeric rows ===
merged_df = pd.merge(mean_prog4, pd1lag3_df[['P_P', 'Mean_PD1_LAG3']], on='P_P', how='inner')
merged_df = merged_df.dropna(subset=['Mean_Program_4', 'Mean_PD1_LAG3'])

# === Remove sample with the highest Mean_Program_4 ===
max_idx = merged_df['Mean_Program_4'].idxmax()
max_pp = merged_df.loc[max_idx, 'P_P']
merged_df = merged_df.drop(index=max_idx)
print(f"Removed highest Mean_Program_4 sample: {max_pp}")

# === Pearson correlation ===
pearson_r, pearson_p = pearsonr(merged_df['Mean_Program_4'], merged_df['Mean_PD1_LAG3'])
formatted_r = f"r = {pearson_r:.2f}"
formatted_p = f"P = {pearson_p:.2f}"  # scientific notation

# === Plot ===
figsize = (8.27 / 2, 11.69 / 4.2)
plt.figure(figsize=figsize)
ax = plt.gca()

# Red regression line with scatter and CI
sns.regplot(
    data=merged_df,
    x='Mean_Program_4',
    y='Mean_PD1_LAG3',
    scatter_kws={'s': 60, 'alpha': 0.7, 'edgecolor': 'k'},
    line_kws={'color': 'red', 'linewidth': 2},
    ci=95,
    ax=ax
)

# Labels and title
ax.set_xlabel('Mean TP4 Activity', fontsize=8)
ax.set_ylabel('Mean PD1 and LAG3 Expression', fontsize=8)
ax.set_title('TP4 Activity vs. PD1 and LAG3 Expression', fontsize=9)
ax.tick_params(axis='both', labelsize=8)
ax.grid(False)

# Correlation text
ax.text(0.98, 0.02, f"{formatted_r}, {formatted_p}",
        transform=ax.transAxes, ha='right', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig("TP4_vs_PD1_LAG3_scatter.pdf", format="pdf", bbox_inches="tight")
plt.show()

