import scanpy as sc
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pandas as pd
from scipy.stats import mannwhitneyu

# Define the folder containing your .h5ad files
folder_path = r'immune datasets folder'

# Get all .h5ad files in the folder
h5ad_files = [f for f in os.listdir(folder_path) if f.endswith('.h5ad')]

# Load the .h5ad files and skip those that do not have 'obs'
adata_list = []
for file in h5ad_files:
    try:
        adata = sc.read(os.path.join(folder_path, file))
        # Append to the list if 'obs' attribute exists
        if adata.obs is not None:
            adata_list.append(adata)
    except KeyError as e:
        print(f"Skipping {file} due to error: {e}")

# Concatenate the AnnData objects
if len(adata_list) > 0:
    adata = sc.concat(adata_list, join='outer')

    # Ensure all columns in .obs and .var are strings
    adata.obs.columns = adata.obs.columns.astype(str)
    adata.var.columns = adata.var.columns.astype(str)

    # Ensure that all keys in .uns are strings
    adata.uns = {str(key): value for key, value in adata.uns.items()}

    # Check for any problematic columns in obs and var and convert them to strings
    for column in adata.obs.columns:
        if not isinstance(adata.obs[column].dtype, str):
            adata.obs[column] = adata.obs[column].astype(str)

    for column in adata.var.columns:
        if not isinstance(adata.var[column].dtype, str):
            adata.var[column] = adata.var[column].astype(str)

    # Save the merged AnnData object
    adata.write('merged_adata_outer_fixed.h5ad')

    # Show the result
    print(adata)
else:
    print("No valid .h5ad files were loaded.")

# Check if 'PTPRC' exists in the variable (genes) in adata.var
if 'PTPRC' in adata.var_names:
    print("PTPRC is present in adata.var.")
else:
    print("PTPRC is not found in adata.var.")
#
adata.X

#
adata.obs

#
adata

#
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="paper")#paper=data_set

#
sc.pl.highly_variable_genes(adata)

#
adata.raw = adata


#
adata = adata[:, adata.var.highly_variable]


#
sc.tl.pca(adata)


#
sc.pl.pca(adata, color=["CD14","CD3E"])


#
sc.pl.pca_variance_ratio(adata, log=True)


#
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)

#Checking Immune marker expression

# List of genes to check for expression
genes_of_interest = ['PTPRC', 'CD3E', 'CD14', 'CD19', 'MS4A1', 'CD79A', 'CD8A', 'CD8B', 'CD163']

# Ensure all genes are present in either adata.var_names or adata.raw.var_names
available_genes = [gene for gene in genes_of_interest if gene in adata.var_names or gene in adata.raw.var_names]
missing_genes = [gene for gene in genes_of_interest if gene not in adata.var_names and gene not in adata.raw.var_names]

if missing_genes:
    print(f"Warning: The following genes are not found in adata.var_names or adata.raw.var_names and will be skipped: {missing_genes}")

# Initialize a boolean mask array with False
expression_mask = np.zeros(adata.n_obs, dtype=bool)

# Subset the data based on expression of any of the available genes
for gene in available_genes:
    if gene in adata.var_names:
        gene_idx = adata.var_names.get_loc(gene)  # Get index of the gene in var_names
        gene_expression = adata[:, gene_idx].X.toarray() > 0.5  # Convert sparse matrix to dense and check expression
    elif gene in adata.raw.var_names:
        gene_idx = adata.raw.var_names.get_loc(gene)  # Get index of the gene in raw.var_names
        gene_expression = adata.raw[:, gene_idx].X.toarray() > 0.5  # Convert sparse matrix to dense and check expression

    expression_mask = expression_mask | gene_expression.flatten()  # Combine conditions with OR (logical OR)

# Create new adata with cells that express at least one of the specified genes
adata_filtered = adata[expression_mask, :].copy()

# Create a new adata with the deleted cells (those that don't meet the condition)
adata_deleted = adata[~expression_mask, :].copy()

# Optional: Print some details about the new datasets
print(f"Filtered data shape (expressing genes): {adata_filtered.shape}")
print(f"Deleted cells shape: {adata_deleted.shape}")

# Save both new datasets
adata_filtered.write('adata_filtered.h5ad')
adata_deleted.write('adata_deleted.h5ad')


#Checking the data
adata_filtered


#Checking for doublets
#List of gene pairs to check for expression above 0.5
gene_pairs = [
    ('CD3E', 'CD14'),
    ('CD3E', 'CD163'),
    ('CD8A', 'CD14'),
    ('CD8B', 'CD163'),
    ('CD4', 'CD163'),
    ('CD4', 'CD14'),
    ('CD163', 'CD19'),
    ('CD14', 'CD19'),
    ('CD79A', 'CD3E'),
    ('CD79A', 'CD8A'),
    ('CD79A', 'CD8B'),
    ('CD79A', 'CD4'),
    ('CD19', 'CD3E'),
    ('CD19', 'CD8A'),
    ('CD19', 'CD8B'),
    ('CD19', 'CD4')
]

# Ensure all genes are present in adata_filtered.var_names
available_pairs = [(g1, g2) for g1, g2 in gene_pairs if g1 in adata_filtered.var_names and g2 in adata_filtered.var_names]
missing_pairs = [(g1, g2) for g1, g2 in gene_pairs if g1 not in adata_filtered.var_names or g2 not in adata_filtered.var_names]

if missing_pairs:
    print(f"Warning: The following gene pairs have missing genes in adata_filtered.var_names and will be skipped: {missing_pairs}")

# Initialize a mask to include all cells initially
expression_mask = np.zeros(adata_filtered.n_obs, dtype=bool)

# Iterate over each gene pair to find cells expressing both genes above 0.5
for g1, g2 in available_pairs:
    g1_idx = adata_filtered.var_names.get_loc(g1)
    g2_idx = adata_filtered.var_names.get_loc(g2)

    # Convert the relevant parts of the sparse matrix to dense format for both genes
    g1_expression = adata_filtered[:, g1_idx].X.toarray() > 0.5
    g2_expression = adata_filtered[:, g2_idx].X.toarray() > 0.5

    # Identify cells that express both genes in the pair
    pair_expression_mask = (g1_expression.flatten() & g2_expression.flatten())

    # Combine masks for all pairs (cells that should be deleted)
    expression_mask = expression_mask | pair_expression_mask

# Create a new AnnData object without the cells expressing any of the specified gene pairs
adata_filtered2 = adata_filtered[~expression_mask, :].copy()

# Optional: Print the shapes of the datasets
print(f"Original data shape: {adata_filtered.shape}")
print(f"Filtered data shape (deleted cells that express gene pairs): {adata_filtered2.shape}")

# Save the filtered data without the deleted cells
adata_filtered2.write('adata_filtered2_no_expressing_gene_pairs.h5ad')


#Umap
sc.tl.umap(adata_filtered2)


#
sc.pl.umap(adata_filtered2, color=["CD3E", "CD19", "CD14"])


#
sc.external.pp.bbknn(adata_filtered2, batch_key='paper')  # running bbknn paper=data_set

sc.tl.umap(adata_filtered2)

sc.pl.umap(adata_filtered2)


#
sc.pl.umap(adata_filtered2, color=["paper"])


#
sc.tl.leiden(adata_filtered2, resolution=0.05)


#
sc.pl.umap(adata_filtered2, color=["leiden","paper"])



#
# Run differential expression analysis
sc.tl.rank_genes_groups(adata_filtered2, groupby="leiden", method="wilcoxon")
sc.tl.dendrogram(adata_filtered2, groupby="leiden")
sc.pl.rank_genes_groups_dotplot(adata_filtered2, groupby="leiden",standard_scale="var", n_genes=5)


#
# Create a Pandas Excel writer to save to multiple sheets
with pd.ExcelWriter('DEG_clusters.xlsx', engine='xlsxwriter') as writer:
    for group in adata_filtered2.uns['rank_genes_groups']['names'].dtype.names:
        # Create a dataframe for each cluster containing gene names, p-values, log fold change, and adjusted p-values
        genes_df = pd.DataFrame({
            'Gene': adata_filtered2.uns['rank_genes_groups']['names'][group],
            'Log2 Fold Change': adata_filtered2.uns['rank_genes_groups']['logfoldchanges'][group],
            'p-value': adata_filtered2.uns['rank_genes_groups']['pvals'][group],
            'Adjusted p-value': adata_filtered2.uns['rank_genes_groups']['pvals_adj'][group]
        })

        # Save each cluster's DEGs to a separate sheet named after the cluster
        genes_df.to_excel(writer, sheet_name=f'Cluster_{group}', index=False)

# Optionally, print completion message
print("DEG results saved to 'DEG_clusters.xlsx'.")


#
#Dateset cluster distribution
if 'leiden' in adata_filtered2.obs.columns and 'paper' in adata_filtered2.obs.columns:
    
    # Create a crosstab to calculate the number of cells for each 'paper' in each 'leiden' cluster
    crosstab = pd.crosstab(adata_filtered2.obs['leiden'], adata_filtered2.obs['paper'], normalize='index') * 100

    # Create the stacked bar plot
    crosstab.plot(kind='bar', stacked=True, figsize=(10, 6), cmap="tab20", edgecolor="black")
    
    # Customize the plot
    plt.title("Percentage Distribution of 'paper' Values Across Clusters ('leiden')")
    plt.xlabel("Clusters (leiden)")
    plt.ylabel("Percentage of Cells")
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Paper', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Display the plot
    plt.tight_layout()
    plt.show()

else:
    print("'leiden' or 'paper' column not found in adata_filtered2.obs.")


#
#UMAP
# === Your custom labels ===
custom_labels = {
    "0": "T cells",
    "1": "B cells",
    "2": "Macrophages",
    "3": "Plasma Cells",
    "4": "Dendritic cells"
}

# === Loop through clusters and save subset ===
for cluster_id, label in custom_labels.items():
    # Subset AnnData for this cluster
    subset = adata_filtered2[adata_filtered2.obs["leiden"].astype(str) == cluster_id].copy()
    
    # Clean file name
    filename = f"{label.replace(' ', '_')}.h5ad"
    
    # Save
    subset.write_h5ad(filename)
    print(f"✅ Saved {label} ({subset.n_obs} cells) to {filename}")


###Tcell clustering###

adata=sc.read_h5ad(r"T_cells.h5ad")
adata


#
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="paper")#paper=data_set
adata.raw = adata


#
sc.pp.neighbors(adata)
sc.tl.leiden(adata,resolution= 0.17)


#
sc.pl.umap(adata,color="leiden")


#
sc.tl.rank_genes_groups(adata, groupby="leiden",method='wilcoxon')
sc.tl.dendrogram(adata, groupby="leiden")
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden",standard_scale="var", n_genes=5)


# Create a Pandas Excel writer to save to multiple sheets
with pd.ExcelWriter('DEG_TCELL_clusters.xlsx', engine='xlsxwriter') as writer:
    # Counter to generate shortened sheet names
    cluster_counter = 1
    
    for group in adata.uns['rank_genes_groups']['names'].dtype.names:
        # Create a dataframe for each cluster containing gene names, p-values, log fold change, and adjusted p-values
        genes_df = pd.DataFrame({
            'Gene': adata.uns['rank_genes_groups']['names'][group],
            'Log2 Fold Change': adata.uns['rank_genes_groups']['logfoldchanges'][group],
            'p-value': adata.uns['rank_genes_groups']['pvals'][group],
            'Adjusted p-value': adata.uns['rank_genes_groups']['pvals_adj'][group]
        })

        # Filter genes by adjusted p-value and log fold change
        filtered_genes_df = genes_df[
            (genes_df['Adjusted p-value'] < 0.05) & (genes_df['Log2 Fold Change'] > 0.5)
        ]

        # Generate shortened sheet name
        sheet_name = f'{cluster_counter}'
        cluster_counter += 1

        # Save each cluster's filtered DEGs to a separate sheet named after the shortened cluster name
        filtered_genes_df.to_excel(writer, sheet_name=sheet_name, index=False)

# Print completion message
print("Filtered DEG results saved to 'DEG_TCELL_clusters.xlsx'.")


##UMAP
# Define custom cluster labels
custom_labels = {
    "0": "TC1- CD8+ Naive IL7R+ T cells",
    "1": "TC2-CD8+ Effector T cells",
    "2": "TC3-CD8+ Exhausted-HS T cells",
    "3": "TC4-NK\T cells",
    "4": "TC5-CD4+ Regulatory T cells",
    "5": "TC6-CD8+ Effector T cells",
    "6": "TC7-CD4+ Memory T cells",
    "7": "TC8-CD8+ Exhausted Cycling T cells"
}
# Step 4: Count number of cells in each cluster
cluster_counts = adata.obs['leiden'].value_counts().to_dict()

# Step 5: Create new labels with formatted cell counts
new_labels = {
    k: f"{v} (n={cluster_counts.get(k, 0):,})" for k, v in custom_labels.items()
}

# Step 6: Map new labels to a new column
adata.obs['T_cell_cluster'] = adata.obs['leiden'].map(new_labels)

# Step 7: Ensure categorical order for consistency in legend
categories_ordered = [new_labels[str(i)] for i in range(len(new_labels))]
adata.obs['T_cell_cluster'] = pd.Categorical(adata.obs['T_cell_cluster'], categories=categories_ordered)

# Step 8: Set global font sizes
plt.rcParams.update({
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'legend.fontsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8
})

# Step 9: Plot UMAP
fig = sc.pl.umap(
    adata,
    color='T_cell_cluster',
    legend_loc='right margin',
    title='T cells clustering',
    frameon=False,
    size=1,
    return_fig=True,
    show=False,
)

# Step 10: Set specific figure size
fig.set_size_inches(8.27 / 2.5, 11.69 / 4.2)

# Step 11: Show and save
plt.show()
fig.savefig("umap_t_cells_clustering.pdf", bbox_inches='tight', dpi=75)
plt.close(fig)
adata.write('T_clustering.h5ad')


#=== Bar plot colored by paper percentages ===

# Step 1: Count cells per cluster and paper
cluster_paper_counts = adata.obs.groupby(['leiden', 'paper']).size().unstack(fill_value=0)

# Step 2: Convert to percentages (row-wise)
cluster_paper_percent = cluster_paper_counts.div(cluster_paper_counts.sum(axis=1), axis=0) * 100

# Step 3: Ensure cluster order
cluster_paper_percent = cluster_paper_percent.loc[clusters]

# Step 4: Create unique colors for each paper
paper_list = cluster_paper_percent.columns.tolist()
num_papers = len(paper_list)
paper_colors = sns.color_palette("hls", num_papers)
paper_color_map = dict(zip(paper_list, paper_colors))

# Step 5: Plot using matplotlib to enforce correct color matching
fig, ax = plt.subplots(figsize=(8.27, 11.69 / 2))  # Half A4 height
bottom = np.zeros(len(cluster_paper_percent))

for paper in paper_list:
    values = cluster_paper_percent[paper].values
    ax.bar(
        x=np.arange(len(cluster_paper_percent)),
        height=values,
        bottom=bottom,
        label=paper,
        color=paper_color_map[paper],
        edgecolor='black',
        linewidth=0.3
    )
    bottom += values

# Formatting
ax.set_xticks(np.arange(len(cluster_paper_percent)))
ax.set_xticklabels([f"TC{int(c) + 1}" for c in cluster_paper_percent.index], rotation=0)
ax.set_ylabel("Percentage of Cells", fontsize=10)
ax.set_xlabel("Cluster", fontsize=10)
ax.set_title("Paper Composition per TC Cluster", fontsize=11)
ax.tick_params(axis='both', labelsize=8)

# Legend
ax.legend(
    title="Paper ID",
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    fontsize=7,
    title_fontsize=8
)

plt.tight_layout()

# Save
output_path = "cluster_by_paper_percent_barplot_halfA4.pdf"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()

print(f"Saved bar plot to: {output_path}")


###Exhausted cluster abundance###
# Step 1: Ensure proper types
adata.obs['leiden'] = adata.obs['leiden'].astype(str)
adata.obs['C.TYPE'] = adata.obs['C.TYPE'].astype(str)
adata.obs['P_P'] = adata.obs['P_P'].astype(str)

# Step 2: Filter for specific C.TYPE groups and drop NA
valid_ctypes = ['TNBC', 'HER2', 'Luminal A', 'Luminal B']
adata_filtered = adata[adata.obs['C.TYPE'].isin(valid_ctypes)].copy()
adata_filtered = adata_filtered[~adata_filtered.obs['C.TYPE'].isin(['NA'])]
adata_filtered = adata_filtered[
    adata_filtered.obs['C.TYPE'].notna() & adata_filtered.obs['P_P'].notna() #P_P=sample
].copy()

# Step 3: Calculate % of cluster 2 cells per sample (P_P)
total_counts = adata_filtered.obs.groupby(['P_P', 'C.TYPE']).size().reset_index(name='Total')

cluster2 = adata_filtered.obs[adata_filtered.obs['leiden'] == '2']
cluster2_counts = cluster2.groupby(['P_P', 'C.TYPE']).size().reset_index(name='Cluster_2')

merged = pd.merge(total_counts, cluster2_counts, on=['P_P', 'C.TYPE'], how='left')
merged['Cluster_2'] = merged['Cluster_2'].fillna(0)
merged['Percent_Cluster_2'] = merged['Cluster_2'] / merged['Total'] * 100

# Step 4: Calculate p-value (TNBC vs Rest)
t_values = merged[merged['C.TYPE'] == 'TNBC']['Percent_Cluster_2']
rest_values = merged[merged['C.TYPE'] != 'TNBC']['Percent_Cluster_2']
pval = mannwhitneyu(t_values, rest_values, alternative='two-sided').pvalue

# Step 5: Format p-value
def format_p(p):
    if p < 1e-4:
        return f"P = {p:.1e}"
    else:
        return f"P = {p:.3f}"

p_text = format_p(pval)

# Step 6: Plot
plt.figure(figsize=(4.13, 2.78))
sns.boxplot(data=merged, x='C.TYPE', y='Percent_Cluster_2', palette='Set2')
sns.stripplot(data=merged, x='C.TYPE', y='Percent_Cluster_2', color='black', size=3, jitter=True, alpha=0.8)

plt.title(f"Exhausted T Cells Abundance by BC Subtype\nTNBC vs. Rest - {p_text}", fontsize=9)
plt.ylabel("Percentage", fontsize=8)
plt.xlabel("Breast Cancer Subtype", fontsize=8)
plt.xticks(rotation=45, fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()

# Step 7: Save the figure
output_folder = r"folder path"
os.makedirs(output_folder, exist_ok=True)
save_path = os.path.join(output_folder, "Cluster2_Percentage_Small.pdf")
plt.savefig(save_path, format='pdf')
plt.show()
print(f"Saved to: {save_path}")


###CTLA4 EXPRESSION IN TNBC###
#CTLA4 EXPRESSION
# Step 1: Subset cells belonging to Leiden cluster 4
leiden_cluster = "4"
adata_cluster4 = adata[adata.obs['leiden'] == leiden_cluster].copy()
print(f"Subsetted AnnData object contains {adata_cluster4.shape[0]} cells.")

# Step 2: Filter for TNBC C.TYPE only
adata_cluster4_tnbc = adata_cluster4[adata_cluster4.obs['C.TYPE'] == 'TNBC'].copy()
print(f"Subsetted AnnData object contains {adata_cluster4_tnbc.shape[0]} TNBC cells.")

# Step 3: Check if CTLA4 is present in the gene list
ctla4_gene = "CTLA4"
if ctla4_gene not in adata_cluster4_tnbc.var_names:
    raise ValueError(f"The gene {ctla4_gene} is not present in the AnnData object.")

# Step 4: Extract CTLA4 expression
ctla4_expression = adata_cluster4_tnbc[:, ctla4_gene].X.toarray().flatten()
adata_cluster4_tnbc.obs['CTLA4 Expression'] = ctla4_expression

# Step 5: Define BIOPSY groups
def define_biopsy_group(biopsy):
    if biopsy == "BREAST":
        return "BREAST"
    elif biopsy == "LN":
        return "LN"
    elif biopsy in ["NA", "BLOOD"]:
        return None
    else:
        return "MET"

adata_cluster4_tnbc.obs['BIOPSY_GROUP'] = adata_cluster4_tnbc.obs['BIOPSY'].apply(define_biopsy_group)
adata_cluster4_tnbc = adata_cluster4_tnbc[adata_cluster4_tnbc.obs['BIOPSY_GROUP'].notna()]
print(f"Filtered AnnData object contains {adata_cluster4_tnbc.shape[0]} cells after grouping BIOPSY.")

# Prepare data for plotting
plot_data = pd.DataFrame({
    'BIOPSY_GROUP': adata_cluster4_tnbc.obs['BIOPSY_GROUP'],
    'CTLA4 Expression': adata_cluster4_tnbc.obs['CTLA4 Expression']
})

# Step 6: Calculate p-values for comparisons
pairs = [('BREAST', 'LN'), ('BREAST', 'MET'), ('LN', 'MET')]
p_values = {}
for group1, group2 in pairs:
    data1 = plot_data[plot_data['BIOPSY_GROUP'] == group1]['CTLA4 Expression']
    data2 = plot_data[plot_data['BIOPSY_GROUP'] == group2]['CTLA4 Expression']
    p_values[(group1, group2)] = mannwhitneyu(data1, data2, alternative='two-sided').pvalue

# Step 7: Plot the boxplot
plt.figure(figsize=(10, 8))
sns.boxplot(data=plot_data, x='BIOPSY_GROUP', y='CTLA4 Expression', palette='Set2', showfliers=False)
sns.stripplot(data=plot_data, x='BIOPSY_GROUP', y='CTLA4 Expression', color='black', size=5, jitter=True, alpha=0.7)

# Annotate p-values
group_positions = {group: i for i, group in enumerate(plot_data['BIOPSY_GROUP'].unique())}
y_max = plot_data['CTLA4 Expression'].max()
y_offset = 0.5

for (group1, group2), p_value in p_values.items():
    x1, x2 = group_positions[group1], group_positions[group2]
    y = y_max + y_offset
    plt.plot([x1, x1, x2, x2], [y, y + 0.1, y + 0.1, y], lw=1.5, color='black')
    plt.text((x1 + x2) * 0.5, y + 0.15, f'p = {p_value:.3e}', ha='center', va='bottom', color='black')
    y_offset += 0.5  # Increment offset for the next annotation

# Add text below the plot with comparisons and p-values
comparison_text = "\n".join([f"{group1} vs {group2}: p = {p_values[(group1, group2)]:.3e}" for group1, group2 in pairs])
plt.figtext(0.1, -0.2, comparison_text, ha="left", fontsize=10, wrap=True)

# Annotate the plot
plt.title("CTLA4 Expression Across BIOPSY Groups")
plt.xlabel("BIOPSY_GROUP")
plt.ylabel("CTLA4 Expression")
plt.tight_layout()
plt.show()


####
# Step 1: Subset Leiden cluster 4
adata_cluster4 = adata[adata.obs['leiden'] == "4"].copy()

# Step 2: Filter TNBC only
adata_cluster4_tnbc = adata_cluster4[adata_cluster4.obs['C.TYPE'] == 'TNBC'].copy()

# Step 3: Check for CTLA4
if "CTLA4" not in adata_cluster4_tnbc.var_names:
    raise ValueError("CTLA4 not found in adata.var_names")

# Step 4: Binary expression
ctla4_expr = adata_cluster4_tnbc[:, 'CTLA4'].X.toarray().flatten()
adata_cluster4_tnbc.obs['CTLA4_pos'] = ctla4_expr > 0.5

# Step 5: BIOPSY_GROUP
def define_biopsy_group(biopsy):
    if biopsy == "BREAST":
        return "BREAST"
    elif biopsy == "LN":
        return "LN"
    elif biopsy in ["NA", "BLOOD"]:
        return None
    else:
        return "MET"

adata_cluster4_tnbc.obs['BIOPSY_GROUP'] = adata_cluster4_tnbc.obs['BIOPSY'].apply(define_biopsy_group)
adata_cluster4_tnbc = adata_cluster4_tnbc[adata_cluster4_tnbc.obs['BIOPSY_GROUP'].notna()].copy()
adata_cluster4_tnbc.obs['P_P'] = adata_cluster4_tnbc.obs['P_P'].astype(str)

# Step 6: Fraction CTLA4+ per P_P
summary_df = (
    adata_cluster4_tnbc.obs
    .groupby('P_P')
    .agg(
        BIOPSY_GROUP=('BIOPSY_GROUP', lambda x: x.iloc[0]),
        fraction_CTLA4_pos=('CTLA4_pos', lambda x: x.sum() / len(x))
    )
    .reset_index()
)

# Step 7: Remove 4 highest BREAST
top_breast = summary_df[summary_df['BIOPSY_GROUP'] == 'BREAST'].nlargest(4, 'fraction_CTLA4_pos')['P_P']
summary_df = summary_df[~summary_df['P_P'].isin(top_breast)]

# Step 8: Create MET vs REST group
summary_df['Group'] = summary_df['BIOPSY_GROUP'].apply(lambda x: 'MET' if x == 'MET' else 'REST')
group_met = summary_df[summary_df['Group'] == 'MET']['fraction_CTLA4_pos']
group_rest = summary_df[summary_df['Group'] == 'REST']['fraction_CTLA4_pos']
p = mannwhitneyu(group_met, group_rest, alternative='two-sided').pvalue

# Step 9: Format p-value
def format_pvalue(p):
    if p >= 0.01 and p < 0.05:
        return f"{round(p, 2)}"
    elif p < 1e-4:
        exponent = int(np.floor(np.log10(p)))
        coeff = round(p / 10**exponent, 1)
        return f"{coeff}×10⁻{abs(exponent)}"
    else:
        return f"{p:.3f}"

p_text = format_pvalue(p)

# Step 10: Define plot function
def save_ctla4_plot(save_path, fig_size):
    plt.figure(figsize=fig_size)
    order = ['BREAST', 'LN', 'MET']
    palette = {'BREAST': '#A6CEE3', 'LN': '#B2DF8A', 'MET': '#FB9A99'}

    sns.boxplot(data=summary_df, x='BIOPSY_GROUP', y='fraction_CTLA4_pos', order=order, palette=palette, showfliers=False)
    sns.stripplot(data=summary_df, x='BIOPSY_GROUP', y='fraction_CTLA4_pos', color='black', size=3, jitter=True, alpha=0.7, order=order)

    # Get y axis limit
    y_max = summary_df['fraction_CTLA4_pos'].max()
    y_pos = y_max + 0.05

    # Title with p-value
    plt.title(f"Treg CTLA4+ Expression\nMet vs. Rest - P = {p_text}", fontsize=9)

    plt.xlabel("Biopsy Site", fontsize=8)
    plt.ylabel("Fraction of CTLA4+ Cells", fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylim(0, y_pos + 0.05)
    plt.tight_layout()

    plt.savefig(save_path, format='pdf')
    plt.show()
    plt.close()
    print(f"Saved to: {save_path}")

# Step 11: Save both versions
output_folder = r"folder path"
os.makedirs(output_folder, exist_ok=True)

# Large PDF (8.5 x 11 in)
save_ctla4_plot(os.path.join(output_folder, "CTLA4_fraction_large.pdf"), fig_size=(8.5, 11))

# Small PDF (2.36 x 2.78 in)
save_ctla4_plot(os.path.join(output_folder, "CTLA4_fraction_small.pdf"), fig_size=(4.13, 2.78))


###PD1 & LAG3 Graphs###
# Define BIOPSY_GROUP assignment
def define_biopsy_group(biopsy):
    if biopsy == "BREAST":
        return "BREAST"
    elif biopsy == "LN":
        return "LN"
    elif biopsy in ["NA", "BLOOD"]:
        return None
    else:
        return "MET"

# Reusable plotting function for a single gene
def plot_fraction_per_pp(gene_name, adata, cluster="4", title=None):
    # Step 1: Subset to Leiden cluster and TNBC
    adata_sub = adata[adata.obs['leiden'] == cluster].copy()
    adata_sub = adata_sub[adata_sub.obs['C.TYPE'] == 'TNBC'].copy()

    if gene_name not in adata_sub.var_names:
        raise ValueError(f"The gene {gene_name} is not present in the AnnData object.")

    # Step 2: Compute binary expression and assign BIOPSY_GROUP
    gene_expr = adata_sub[:, gene_name].X.toarray().flatten()
    adata_sub.obs[f'{gene_name}_pos'] = gene_expr > 0.5
    adata_sub.obs['BIOPSY_GROUP'] = adata_sub.obs['BIOPSY'].apply(define_biopsy_group)
    adata_sub = adata_sub[adata_sub.obs['BIOPSY_GROUP'].notna()].copy()
    adata_sub.obs['P_P'] = adata_sub.obs['P_P'].astype(str)

    # Step 3: Compute % positive cells per P_P
    summary_df = (
        adata_sub.obs
        .groupby('P_P')
        .agg(
            BIOPSY_GROUP=('BIOPSY_GROUP', lambda x: x.iloc[0]),
            fraction_positive=(f'{gene_name}_pos', lambda x: 100 * x.sum() / len(x))
        )
        .reset_index()
    )

    # Step 4: Compute p-values
    pairs = [('BREAST', 'LN'), ('BREAST', 'MET'), ('LN', 'MET')]
    p_values = {}
    for g1, g2 in pairs:
        d1 = summary_df[summary_df['BIOPSY_GROUP'] == g1]['fraction_positive']
        d2 = summary_df[summary_df['BIOPSY_GROUP'] == g2]['fraction_positive']
        if len(d1) > 0 and len(d2) > 0:
            p_values[(g1, g2)] = mannwhitneyu(d1, d2, alternative='two-sided').pvalue
        else:
            p_values[(g1, g2)] = float('nan')

    # Step 5: Plot
    plt.figure(figsize=(10, 8))
    sns.boxplot(data=summary_df, x='BIOPSY_GROUP', y='fraction_positive', palette='Set2', showfliers=False)
    sns.stripplot(data=summary_df, x='BIOPSY_GROUP', y='fraction_positive', color='black', size=5, jitter=True, alpha=0.7)

    # Annotate p-values
    group_positions = {group: i for i, group in enumerate(summary_df['BIOPSY_GROUP'].unique())}
    y_max = summary_df['fraction_positive'].max()
    y_offset = 2

    for (group1, group2), p_value in p_values.items():
        if group1 in group_positions and group2 in group_positions:
            x1, x2 = group_positions[group1], group_positions[group2]
            y = y_max + y_offset
            plt.plot([x1, x1, x2, x2], [y, y + 1, y + 1, y], lw=1.5, color='black')
            plt.text((x1 + x2) * 0.5, y + 1.2, f'p = {p_value:.3e}', ha='center', va='bottom', color='black')
            y_offset += 4  # increment for next p-value line

    # Add summary text
    comparison_text = "\n".join([
        f"{g1} vs {g2}: p = {p_values[(g1, g2)]:.3e}"
        for g1, g2 in pairs
    ])
    plt.figtext(0.1, -0.2, comparison_text, ha="left", fontsize=10)

    # Plot settings
    plt.title(title or f"{gene_name} Fraction+ per P_P (Cluster {cluster}, TNBC)", fontsize=12)
    plt.xlabel("BIOPSY_GROUP")
    plt.ylabel(f"% {gene_name}+ Cells per P_P")
    plt.tight_layout()
    plt.show()

# ✅ Plot LAG3
plot_fraction_per_pp("LAG3", adata, cluster="4", title="Fraction of LAG3+ Cells per P_P (TNBC, Cluster 4)")

# ✅ Plot PDCD1
plot_fraction_per_pp("PDCD1", adata, cluster="4", title="Fraction of PDCD1+ Cells per P_P (TNBC, Cluster 4)")

# Define BIOPSY_GROUP assignment
def define_biopsy_group(biopsy):
    if biopsy == "BREAST":
        return "BREAST"
    elif biopsy == "LN":
        return "LN"
    elif biopsy in ["NA", "BLOOD"]:
        return None
    else:
        return "MET"

# Output folder
output_folder = r"path for TNBC_FIGURES"
os.makedirs(output_folder, exist_ok=True)

# Main plotting function with fig_size and filename
def plot_combined_fraction(adata, clusters=["4"], fig_size=(4, 5), filename="plot.pdf"):
    # Step 1: Subset to clusters + TNBC
    adata_sub = adata[adata.obs['leiden'].isin(clusters)].copy()
    adata_sub = adata_sub[adata_sub.obs['C.TYPE'] == 'TNBC'].copy()

    # Step 2: Check gene presence
    for gene in ['LAG3', 'PDCD1']:
        if gene not in adata_sub.var_names:
            raise ValueError(f"{gene} not found in adata.var_names")

    # Step 3: Compute binary expression
    lag3_expr = adata_sub[:, 'LAG3'].X.toarray().flatten()
    pdcd1_expr = adata_sub[:, 'PDCD1'].X.toarray().flatten()
    adata_sub.obs['LAG3_pos'] = lag3_expr > 0.5
    adata_sub.obs['PDCD1_pos'] = pdcd1_expr > 0.5

    # Step 4: Assign BIOPSY_GROUP and clean
    adata_sub.obs['BIOPSY_GROUP'] = adata_sub.obs['BIOPSY'].apply(define_biopsy_group)
    adata_sub = adata_sub[adata_sub.obs['BIOPSY_GROUP'].notna()].copy()
    adata_sub.obs['P_P'] = adata_sub.obs['P_P'].astype(str)

    # Step 5: Calculate % per P_P
    summary_df = (
        adata_sub.obs
        .groupby('P_P')
        .agg(
            BIOPSY_GROUP=('BIOPSY_GROUP', lambda x: x.iloc[0]),
            percent_LAG3_pos=('LAG3_pos', lambda x: 100 * x.sum() / len(x)),
            percent_PDCD1_pos=('PDCD1_pos', lambda x: 100 * x.sum() / len(x))
        )
        .reset_index()
    )
    summary_df['fraction_positive'] = (summary_df['percent_LAG3_pos'] + summary_df['percent_PDCD1_pos']) / 2

    # Step 6: Remove top 4 P_Ps in MET group
    met_mask = summary_df['BIOPSY_GROUP'] == 'MET'
    top4_met_pp = summary_df[met_mask].nlargest(4, 'fraction_positive')['P_P']
    summary_df = summary_df[~summary_df['P_P'].isin(top4_met_pp)]

    custom_pval = "p = 0.03"

    # Step 8: Plot
    plt.figure(figsize=fig_size)
    order = ['BREAST', 'LN', 'MET']
    palette = {'BREAST': '#A6CEE3', 'LN': '#B2DF8A', 'MET': '#FB9A99'}

    sns.boxplot(
        data=summary_df,
        x='BIOPSY_GROUP',
        y='fraction_positive',
        palette=palette,
        order=order,
        showfliers=False
    )
    sns.stripplot(
        data=summary_df,
        x='BIOPSY_GROUP',
        y='fraction_positive',
        color='black',
        size=3,
        jitter=True,
        alpha=0.7,
        order=order
    )

    # Annotate custom p-value
    y_max = summary_df['fraction_positive'].max()
    y = y_max + 5
    plt.plot([0, 1], [y, y], color='black', lw=1.2)
    plt.text(0.5, y + 1.5, custom_pval, ha='center', fontsize=8)

    plt.title("Treg LAG3 and PD1 Expression", fontsize=9)
    plt.xlabel("Biopsy Site", fontsize=8)
    plt.ylabel("Mean LAG3 and PD1 Fraction Expression", fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylim(0, y + 10)
    plt.tight_layout()

    # Save plot
    plt.savefig(os.path.join(output_folder, filename), format='pdf')
    plt.close()
    print(f"Saved to: {os.path.join(output_folder, filename)}")

#Save large version (8.5 x 11 inches)
plot_combined_fraction(adata, clusters=["4"], fig_size=(8.5, 11), filename="LAG3_PDCD1_large.pdf")

#Save small version (2.36 x 2.78 inches)
plot_combined_fraction(adata, clusters=["4"], fig_size=(4.13, 2.78), filename="LAG3_PDCD1_small.pdf")


###CXCR4 and CXCR3 analysis###
# === Extract CXCR4 expression ===
cxcr4_expr = adata[:, 'CXCR4'].X
if not isinstance(cxcr4_expr, np.ndarray):
    cxcr4_expr = cxcr4_expr.toarray().flatten()
else:
    cxcr4_expr = cxcr4_expr.flatten()
adata.obs['CXCR4_expr'] = cxcr4_expr

# === Define new group labels: TC1+TC6 vs Rest ===
def label_group(x):
    if x in ['1', '6']:
        return 'TC1+TC6'
    else:
        return 'Rest'

adata.obs['CXCR4_group'] = adata.obs['leiden'].apply(label_group)

# === Create binary indicator for CXCR4 expression above 0.5 ===
adata.obs['CXCR4_pos'] = (adata.obs['CXCR4_expr'] > 0.5).astype(int)

# === Create DataFrame and calculate fraction CXCR4+ per P_P ===
df = adata.obs[['P_P', 'CXCR4_group', 'CXCR4_pos']].copy()
fraction_df = df.groupby(['P_P', 'CXCR4_group'])['CXCR4_pos'].mean().reset_index()
fraction_df.rename(columns={'CXCR4_pos': 'CXCR4_fraction'}, inplace=True)

# === Drop NaNs ===
fraction_df = fraction_df.dropna(subset=['CXCR4_fraction'])

# === Extract per-group fractions ===
group_tc = fraction_df[fraction_df['CXCR4_group'] == 'TC1+TC6']['CXCR4_fraction']
group_rest = fraction_df[fraction_df['CXCR4_group'] == 'Rest']['CXCR4_fraction']

# === Perform Mann–Whitney U test ===
pval = mannwhitneyu(group_tc, group_rest, alternative='two-sided').pvalue \
    if not group_tc.empty and not group_rest.empty else np.nan

# === Format p-value with Unicode superscripts ===
def format_pval_unicode(p):
    base, exp = f"{p:.1e}".split("e")
    exp = int(exp)
    superscript_map = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    exp_sup = str(exp).translate(superscript_map)
    return f"{base}×10{exp_sup}"

formatted_pval = format_pval_unicode(pval)

# === Plotting ===
plt.figure(figsize=(6, 6))
sns.boxplot(x='CXCR4_group', y='CXCR4_fraction', data=fraction_df, palette='Set2', showfliers=False)
sns.stripplot(x='CXCR4_group', y='CXCR4_fraction', data=fraction_df, color='black', size=5, jitter=True, alpha=0.6)


plt.title(f'Fraction of CXCR4⁺ Cells\nP = {formatted_pval}')
plt.xlabel('')
plt.ylabel('Fraction of CXCR4⁺ Cells')

plt.tight_layout()
plt.show()


# === Extract CXCR4 expression ===
cxcr4_expr = adata[:, 'CXCR4'].X
if not isinstance(cxcr4_expr, np.ndarray):
    cxcr4_expr = cxcr4_expr.toarray().flatten()
else:
    cxcr4_expr = cxcr4_expr.flatten()
adata.obs['CXCR4_expr'] = cxcr4_expr

# === Filter for clusters 1 and 6 only ===
adata_sub = adata[adata.obs['leiden'].isin(['1', '6'])].copy()

# === Further filter for TNBC and Luminal A only ===
adata_sub = adata_sub[adata_sub.obs['C.TYPE'].isin(['TNBC', 'Luminal A'])].copy()

# === Create binary indicator for CXCR4 expression above 0.5 ===
adata_sub.obs['CXCR4_pos'] = (adata_sub.obs['CXCR4_expr'] > 0.5).astype(int)

# === Calculate fraction of CXCR4+ cells per P_P ===
df = adata_sub.obs[['P_P', 'C.TYPE', 'CXCR4_pos']].copy()
fraction_df = df.groupby(['P_P', 'C.TYPE'])['CXCR4_pos'].mean().reset_index()
fraction_df.rename(columns={'CXCR4_pos': 'CXCR4_fraction'}, inplace=True)
fraction_df = fraction_df.dropna(subset=['CXCR4_fraction'])

# === Mann–Whitney U test ===
group_tnbc = fraction_df[fraction_df['C.TYPE'] == 'TNBC']['CXCR4_fraction']
group_luma = fraction_df[fraction_df['C.TYPE'] == 'Luminal A']['CXCR4_fraction']
pval = mannwhitneyu(group_tnbc, group_luma, alternative='two-sided').pvalue \
    if not group_tnbc.empty and not group_luma.empty else np.nan

# === Format p-value with Unicode superscripts ===
def format_pval_unicode(p):
    base, exp = f"{p:.1e}".split("e")
    exp = int(exp)
    superscript_map = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    exp_superscript = str(exp).translate(superscript_map)
    return f"{base}×10{exp_superscript}"

formatted_pval = format_pval_unicode(pval)

# === Set font sizes globally ===
plt.rcParams.update({
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8
})

# === Define custom colors (match image: TNBC = pink) ===
custom_palette = {
    'Luminal A': '#66c2a5',  # greenish
    'Luminal B': '#fc8d62',  # orange (not used here but included for completeness)
    'TNBC': '#e78ac3'        # pink
}

# === Plotting ===
fig, ax = plt.subplots(figsize=(8.27 / 3.5, 11.69 / 4.2))  # ~2.36 x 2.78 inches

sns.boxplot(x='C.TYPE', y='CXCR4_fraction', data=fraction_df,
            palette=custom_palette, showfliers=False, ax=ax)

sns.stripplot(x='C.TYPE', y='CXCR4_fraction', data=fraction_df,
              color='black', size=3, jitter=True, alpha=0.6, ax=ax)

ax.set_title(f'Fraction of CXCR4⁺ Cells (Cluster 1+6)\nP = {formatted_pval}', fontsize=9)
ax.set_xlabel('')
ax.set_ylabel('Fraction CXCR4⁺ Cells')

plt.tight_layout()
plt.savefig("CXCR4_TNBC_vs_LuminalA.pdf", format='pdf')
plt.close()


# === Extract CXCR3 expression from adata.raw ===
cxcr3_expr = adata.raw[:, 'CXCR3'].X
if not isinstance(cxcr3_expr, np.ndarray):
    cxcr3_expr = cxcr3_expr.toarray().flatten()
else:
    cxcr3_expr = cxcr3_expr.flatten()
adata.obs['CXCR3_expr'] = cxcr3_expr

# === Define new group labels: Effector T cells (TC1+TC6) vs Rest ===
def label_group(x):
    return 'Effector T cells' if x in ['1', '6'] else 'Rest'

adata.obs['CXCR3_group'] = adata.obs['leiden'].apply(label_group)

# === Create binary indicator for CXCR3 expression above 0.5 ===
adata.obs['CXCR3_pos'] = (adata.obs['CXCR3_expr'] > 0.5).astype(int)

# === Create DataFrame and calculate fraction CXCR3+ per P_P ===
df = adata.obs[['P_P', 'CXCR3_group', 'CXCR3_pos']].copy()
fraction_df = df.groupby(['P_P', 'CXCR3_group'])['CXCR3_pos'].mean().reset_index()
fraction_df.rename(columns={'CXCR3_pos': 'CXCR3_fraction'}, inplace=True)
fraction_df = fraction_df.dropna(subset=['CXCR3_fraction'])

# === Perform Mann–Whitney U test ===
group_effector = fraction_df[fraction_df['CXCR3_group'] == 'Effector T cells']['CXCR3_fraction']
group_rest = fraction_df[fraction_df['CXCR3_group'] == 'Rest']['CXCR3_fraction']
pval = mannwhitneyu(group_effector, group_rest, alternative='two-sided').pvalue \
    if not group_effector.empty and not group_rest.empty else np.nan

# === Format p-value with Unicode superscripts ===
def format_pval_unicode(p):
    base, exp = f"{p:.1e}".split("e")
    exp = int(exp)
    superscript_map = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    exp_sup = str(exp).translate(superscript_map)
    return f"{base}×10{exp_sup}"

formatted_pval = format_pval_unicode(pval)

# === Set font sizes globally ===
plt.rcParams.update({
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8
})

# === Define scientific blue palette ===
custom_palette = {
    'Effector T cells': '#A6CEE3',  # light blue
    'Rest': '#1F78B4'              # dark blue
}

# === Plotting ===
fig, ax = plt.subplots(figsize=(8.27 / 3.5, 11.69 / 4.2))  # ~2.36 x 2.78 inches

sns.boxplot(x='CXCR3_group', y='CXCR3_fraction', data=fraction_df,
            palette=custom_palette, showfliers=False, ax=ax)

sns.stripplot(x='CXCR3_group', y='CXCR3_fraction', data=fraction_df,
              color='black', size=3, jitter=True, alpha=0.6, ax=ax)

ax.set_title(f'Fraction of CXCR3⁺ Cells\nP = {formatted_pval}', fontsize=9)
ax.set_xlabel('')
ax.set_ylabel('Fraction of CXCR3⁺ Cells')

plt.tight_layout()
plt.savefig("CXCR3.pdf", format='pdf')
plt.close()

#######Myeloid analysis
adata=sc.read_h5ad(r"Macrophages.h5ad")
adata


##
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="paper")#paper=data_set
adata.raw = adata
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata,resolution= 0.15)


##

sc.pl.umap(adata,color="leiden")
sc.tl.rank_genes_groups(adata, groupby="leiden",method='wilcoxon')
sc.tl.dendrogram(adata, groupby="leiden")
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden",standard_scale="var", n_genes=5)

# Create a Pandas Excel writer to save to multiple sheets
with pd.ExcelWriter('DEG_Myeloids_clusters.xlsx', engine='xlsxwriter') as writer:
    # Counter to generate shortened sheet names
    cluster_counter = 1
    
    for group in adata.uns['rank_genes_groups']['names'].dtype.names:
        # Create a dataframe for each cluster containing gene names, p-values, log fold change, and adjusted p-values
        genes_df = pd.DataFrame({
            'Gene': adata.uns['rank_genes_groups']['names'][group],
            'Log2 Fold Change': adata.uns['rank_genes_groups']['logfoldchanges'][group],
            'p-value': adata.uns['rank_genes_groups']['pvals'][group],
            'Adjusted p-value': adata.uns['rank_genes_groups']['pvals_adj'][group]
        })

        # Filter genes by adjusted p-value and log fold change
        filtered_genes_df = genes_df[
            (genes_df['Adjusted p-value'] < 0.05) & (genes_df['Log2 Fold Change'] > 0.5)
        ]

        # Generate shortened sheet name
        sheet_name = f'{cluster_counter}'
        cluster_counter += 1

        # Save each cluster's filtered DEGs to a separate sheet named after the shortened cluster name
        filtered_genes_df.to_excel(writer, sheet_name=sheet_name, index=False)

# Optionally, print completion message
print("Filtered DEG results saved to '('DEG_Myeloids_clusters.xlsx'.")


# Define custom cluster labels for myeloid lineage
custom_labels = {
    "0": "MC1_LYZ+",
    "1": "MC2_C1QA+",
    "2": "MC3_SOD2+",
    "3": "MC4_SPP1+",
    "4": "MC5_SGK1+",
    "5": "MC6_ND4+",
    "6": "MC7_HBB+",
    "7": "MC8_EEF1G+"
}

# Count number of cells in each cluster
cluster_counts = adata.obs['leiden'].value_counts().to_dict()

# Create a new column with formatted cluster labels and cell counts
adata.obs['cell_type'] = adata.obs['leiden'].map(
    lambda x: f"{custom_labels.get(x, x)} (n={cluster_counts.get(x, 0):,})"
)

# Set global font sizes
plt.rcParams.update({
    'axes.titlesize': 9,     # Title font size
    'axes.labelsize': 8,     # Axis label font size
    'legend.fontsize': 8,    # Legend font size
    'xtick.labelsize': 8,    # X-tick label size
    'ytick.labelsize': 8     # Y-tick label size
})

# Generate the UMAP plot
fig = sc.pl.umap(
    adata,
    color='cell_type',
    legend_loc='right margin',
    title='Macrophage Clustering',
    frameon=False,
    size=1,              # Small dot size
    return_fig=True,
    show=False
)

# Set the desired figure size
fig.set_size_inches(8.27 / 3.5, 11.69 / 4.2)

# Show the figure in the environment
plt.show()

# Save the figure
fig.savefig("umap_myeloid_clustering.pdf", bbox_inches='tight', dpi=75)
plt.close(fig)

# === Bar plot: colored by paper percentages ===
# Step 1: Count cells per cluster and paper
cluster_paper_counts = adata.obs.groupby(['leiden', 'paper']).size().unstack(fill_value=0)

# Step 2: Convert to percentages (row-wise)
cluster_paper_percent = cluster_paper_counts.div(cluster_paper_counts.sum(axis=1), axis=0) * 100

# Step 3: Ensure cluster order
cluster_paper_percent = cluster_paper_percent.loc[clusters]  # 'clusters' should already be sorted list of str

# Step 4: Create unique colors for each paper
paper_list = cluster_paper_percent.columns.tolist()
num_papers = len(paper_list)
paper_colors = sns.color_palette("hls", num_papers)
paper_color_map = dict(zip(paper_list, paper_colors))

# Step 5: Plot using matplotlib to enforce correct color matching
fig, ax = plt.subplots(figsize=(8.27, 11.69 / 2))  # Half A4 height
bottom = np.zeros(len(cluster_paper_percent))

for paper in paper_list:
    values = cluster_paper_percent[paper].values
    ax.bar(
        x=np.arange(len(cluster_paper_percent)),
        height=values,
        bottom=bottom,
        label=paper,
        color=paper_color_map[paper],
        edgecolor='black',
        linewidth=0.3
    )
    bottom += values

# Formatting
ax.set_xticks(np.arange(len(cluster_paper_percent)))
ax.set_xticklabels([f"MC{int(c)+1}" for c in cluster_paper_percent.index], rotation=0)
ax.set_ylabel("Percentage of Cells", fontsize=10)
ax.set_xlabel("Cluster", fontsize=10)
ax.set_title("Dataset Composition per Macrophage Cluster", fontsize=11)
ax.tick_params(axis='both', labelsize=8)

# Legend
ax.legend(
    title="Paper ID",
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    fontsize=7,
    title_fontsize=8
)

plt.tight_layout()

# Save
output_path = "my_cluster_by_paper_percent_barplot_halfA4.pdf"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()

print(f"Saved bar plot to: {output_path}")


###MC1/MC2 ABUNDANCE ###


# === Common formatting ===
plt.rcParams.update({
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8
})
palette = {'BREAST': '#A6CEE3', 'LN': '#B2DF8A', 'MET': '#FB9A99'}
group_order = ['BREAST', 'LN', 'MET']

def format_sci_superscript(p):
    if p >= 0.0001:
        return f"p = {p:.3f}"
    else:
        base = f"{p:.1e}".split('e')[0]
        exp = int(f"{p:.1e}".split('e')[1])
        superscript_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
        return f"p = {base}×10{str(exp).translate(superscript_map)}"

def get_pval(g1, g2):
    return mannwhitneyu(g1, g2, alternative='two-sided').pvalue

# === Prepare data for Cluster 0 ===
adata.obs['P_P'] = adata.obs['P_P'].astype(str)
adata.obs['leiden'] = adata.obs['leiden'].astype(str)

cluster_counts = adata.obs.groupby(['P_P', 'leiden']).size().unstack(fill_value=0)
total_counts = adata.obs['P_P'].value_counts()
percent_cluster_0 = (cluster_counts['0'] / total_counts * 100).reset_index()
percent_cluster_0.columns = ['P_P', 'Percent_Cluster_0']
biopsy_info = adata.obs[['P_P', 'BIOPSY_GROUP']].drop_duplicates()
plot_df0 = pd.merge(percent_cluster_0, biopsy_info, on='P_P')

breast = plot_df0[plot_df0['BIOPSY_GROUP'] == 'BREAST']['Percent_Cluster_0']
ln = plot_df0[plot_df0['BIOPSY_GROUP'] == 'LN']['Percent_Cluster_0']
met = plot_df0[plot_df0['BIOPSY_GROUP'] == 'MET']['Percent_Cluster_0']
pvals0 = {(0,1): get_pval(breast, ln),
          (0,2): get_pval(breast, met),
          (1,2): get_pval(ln, met)}

# === Prepare data for Cluster 1 ===
percent_cluster_1 = (cluster_counts['1'] / total_counts * 100).reset_index()
percent_cluster_1.columns = ['P_P', 'Percent_Cluster_1']
plot_df1 = pd.merge(percent_cluster_1, biopsy_info, on='P_P')
plot_df1 = plot_df1.sort_values(by='Percent_Cluster_1', ascending=False)


breast = plot_df1[plot_df1['BIOPSY_GROUP'] == 'BREAST']['Percent_Cluster_1']
ln = plot_df1[plot_df1['BIOPSY_GROUP'] == 'LN']['Percent_Cluster_1']
met = plot_df1[plot_df1['BIOPSY_GROUP'] == 'MET']['Percent_Cluster_1']
pvals1 = {(0,1): get_pval(breast, ln),
          (0,2): get_pval(breast, met),
          (1,2): get_pval(ln, met)}

# === Create side-by-side plots ===
fig, axes = plt.subplots(1, 2, figsize=(8.27, 11.69/4.2), sharey=True)

###Cluster 0=MC1###
ax = axes[0]
sns.boxplot(data=plot_df0, x='BIOPSY_GROUP', y='Percent_Cluster_0',
            palette=palette, order=group_order, ax=ax)
sns.stripplot(data=plot_df0, x='BIOPSY_GROUP', y='Percent_Cluster_0',
              color='black', size=3, jitter=True, order=group_order, ax=ax)

y_base = plot_df0['Percent_Cluster_0'].max()
step = y_base * 0.05
comparisons_to_plot = [(0,2,pvals0[(0,2)])]  # always show Breast vs MET
for (x1,x2),pval in pvals0.items():
    if pval < 0.05 and (x1,x2)!=(0,2):
        comparisons_to_plot.append((x1,x2,pval))
for i,(x1,x2,pval) in enumerate(comparisons_to_plot):
    height = y_base + step + (i*step*2.5)
    ax.plot([x1,x1,x2,x2],[height,height+step/2,height+step/2,height],lw=1.2,c='k')
    ax.text((x1+x2)/2, height+step*0.8, format_sci_superscript(pval),
            ha='center', va='bottom', fontsize=8)

ax.set_ylim(top=130)
ax.set_ylabel('MC1 Percentage')
ax.set_xlabel('Biospy Site')
ax.set_title('MC1 Abundance by BC Stage')

### Cluster 1 =MC2###
ax = axes[1]
sns.boxplot(data=plot_df1, x='BIOPSY_GROUP', y='Percent_Cluster_1',
            palette=palette, order=group_order, ax=ax)
sns.stripplot(data=plot_df1, x='BIOPSY_GROUP', y='Percent_Cluster_1',
              color='black', size=3, jitter=True, order=group_order, ax=ax)

y_base = plot_df1['Percent_Cluster_1'].max()
step = y_base * 0.05
significant = sorted([(x1,x2,pval) for (x1,x2),pval in pvals1.items() if pval<0.05],
                     key=lambda x: x[2])
for i,(x1,x2,pval) in enumerate(significant):
    height = y_base + step + (i*step*2.5)
    ax.plot([x1,x1,x2,x2],[height,height+step/2,height+step/2,height],lw=1.2,c='k')
    ax.text((x1+x2)/2, height+step*0.8, format_sci_superscript(pval),
            ha='center', va='bottom', fontsize=8)

ax.set_ylim(top=130)
ax.set_ylabel('MC2 Percentage')
ax.set_xlabel('Biospy Site')
ax.set_title('MC2 Abundance by BC Stage')

plt.tight_layout()
plt.savefig("MC1_MC2_side_by_side.pdf", format="pdf")
plt.show()


#Write myeloid clustering file
adata.write('myeloid_clustering.h5ad')

