#Stroma integration and analysis#

import scanpy as sc
import os
import numpy as np
import scipy.sparse as sp
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import liana as li
import matplotlib as mpl
from scipy.stats import mannwhitneyu
import math

#'C.TYPE' =BC Subtype
if 'C.TYPE' in adata_filtered.obs.columns:
    # Replace NaN with 'NA'
    adata_filtered.obs['C.TYPE'] = adata_filtered.obs['C.TYPE'].replace(pd.NA, 'NA')
    
    # Replace 'LUMINALA' with 'Luminal A' and 'LUMINALB' with 'Luminal B'
    adata_filtered.obs['C.TYPE'] = adata_filtered.obs['C.TYPE'].replace({
        'LUMINALA': 'Luminal A',
        'LUMINALB': 'Luminal B'
    })

    # Print the updated 'C.TYPE' values to confirm
    print(adata_filtered.obs['C.TYPE'].value_counts())
else:
    print("'C.TYPE' column not found in adata_filtered.obs.")
# Define the folder containing your .h5ad files
folder_path = r'Stroma folder'

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

    print(adata)
else:
    print("No valid .h5ad files were loaded.")


#Data check
adata.obs

#
adata

#
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="paper")

#
sc.pl.highly_variable_genes(adata)

#
adata.raw = adata

#
adata = adata[:, adata.var.highly_variable]

#
sc.tl.pca(adata)

#
sc.pl.pca(adata, color="FAP")

#
sc.pl.pca_variance_ratio(adata, log=True)


#
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)


#
sc.tl.umap(adata)


#
sc.pl.umap(adata, color=["paper"])


#Checking for immune cells markers
# List of genes to check for expression
genes_of_interest = ['PTPRC', 'CD3E', 'CD14', 'CD19', 'MS4A1', 'CD79A', 'CD8A', 'CD8B', 'CD163']

# Ensure the genes are in adata.raw.var_names
genes_found = [gene for gene in genes_of_interest if gene in adata.raw.var_names]

if len(genes_found) == 0:
    print("None of the genes were found in adata.raw.")
else:
    # Get the expression data for the selected genes from adata.raw
    expression_data = adata.raw[:, genes_found].X
    
    # Handle sparse matrices and convert them to dense if necessary
    if sp.issparse(expression_data):
        expression_data = expression_data.toarray()
    
    # Create a boolean mask for cells expressing any of the genes of interest
    expression_mask = np.any(expression_data > 0.5, axis=1)
    
    # Count the number of cells to be deleted
    num_cells_to_delete = expression_mask.sum()

    # Subset the AnnData object to keep only cells that do not express the genes
    adata_filtered = adata[~expression_mask].copy()

    # Calculate how many cells were deleted
    num_cells_deleted = adata.shape[0] - adata_filtered.shape[0]
    print(f"Number of cells deleted: {num_cells_deleted}")


#
sc.tl.umap(adata_filtered)

#
sc.pl.umap(adata_filtered, color=["FAP"])

#
sc.external.pp.bbknn(adata_filtered, batch_key="paper")

#
sc.tl.umap(adata_filtered)
sc.pl.umap(adata_filtered, color=["paper"])

#
sc.pl.umap(adata_filtered,color=['VWF','COL1A1','ACTA2','PDGFRB','MCAM','CD34','PECAM1','GATA3'])


##Clustering of the stromal cells to:CAF,Endothelial cells, Luminal epithelial cells, PVS.

sc.tl.leiden(adata_filtered,resolution=0.1)

#
sc.pl.umap(adata_filtered,color=["leiden"])

# Perform differential expression analysis
sc.tl.rank_genes_groups(adata_filtered, groupby="leiden", method="wilcoxon")
sc.tl.dendrogram(adata_filtered, groupby="leiden")
sc.pl.rank_genes_groups_dotplot(adata_filtered, groupby="leiden",standard_scale="var", n_genes=5)


# Create a Pandas Excel writer to save to multiple sheets
with pd.ExcelWriter('DEG_STROMA_CLUSTERS.xlsx', engine='xlsxwriter') as writer:
    # Counter to generate shortened sheet names
    cluster_counter = 1
    
    for group in adata_filtered.uns['rank_genes_groups']['names'].dtype.names:
        # Create a dataframe for each cluster containing gene names, p-values, log fold change, and adjusted p-values
        genes_df = pd.DataFrame({
            'Gene': adata_filtered.uns['rank_genes_groups']['names'][group],
            'Log2 Fold Change': adata_filtered.uns['rank_genes_groups']['logfoldchanges'][group],
            'p-value': adata_filtered.uns['rank_genes_groups']['pvals'][group],
            'Adjusted p-value': adata_filtered.uns['rank_genes_groups']['pvals_adj'][group]
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


print("Filtered DEG results saved to 'DEG_STROMA_CLUSTERS.xlsx'.")
adata_filtered.write('stroma&epithel_clustering.h5ad')


#UMAP

# Define custom cluster name mapping
custom_names = {
    '0': 'CAF_1',
    '1': 'Perivascula cells',
    '2': 'CAF_2',
    '3': 'Normal epithelial cluster',
    '4': 'Endothelial cells'
}

# Count number of cells per leiden cluster
cluster_counts = adata_filtered.obs['leiden'].value_counts().sort_index()

# Create new labels with custom names and cell counts, adding "n="
new_labels = {
    cluster: f"{custom_names.get(cluster, cluster)} (n={count})"
    for cluster, count in cluster_counts.items()
}

# Add the new labels to adata.obs
adata_filtered.obs['leiden_named'] = adata_filtered.obs['leiden'].map(new_labels)

# Plot UMAP with updated labels
sc.pl.umap(adata_filtered, color='leiden_named', legend_loc='right margin')

# Dateset cluster distribution
if 'leiden' in adata_filtered.obs.columns and 'paper' in adata_filtered.obs.columns:
    
    # Create a crosstab to calculate the number of cells for each 'paper' in each 'leiden' cluster
    crosstab = pd.crosstab(adata_filtered.obs['leiden'], adata_filtered.obs['paper'], normalize='index') * 100

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


###Further fine clustering of the stromal cells###
adata=sc.read_h5ad("stroma&epithel_clustering.h5ad")

sc.tl.leiden(adata,resolution=0.5)

# Plot UMAP with updated labels and title
sc.pl.umap(
    adata,
    color='leiden',
    legend_loc='right margin',
    title='Stroma & Epithelia clustering'
)

# Run rank_genes_groups using simplified cluster names
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
sc.tl.dendrogram(adata, groupby="leiden")
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden",standard_scale="var", n_genes=5)


# Create a Pandas Excel writer to save to multiple sheets
with pd.ExcelWriter('DEG_STROMA_FINE_CLUSTERING.xlsx', engine='xlsxwriter') as writer:
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
print("Filtered DEG results saved to 'DEG_STROMA_FINE_CLUSTERING.xlsx'.")

#UMAP
# Define custom cluster name mapping
custom_names = {
    '0': 'CAF1-ECM-Depositing ',
    '1': 'PVC1-CD36+ Immature',
    '2': 'CAF2-EMT-like',
    '3': 'Normal Luminal Epithel',
    '4': 'CAF3-Infiltrative',
    '5': 'CAF4-Myofibroblast-Like',
    '6': 'CAF5-Secretory ECM',
    '7': 'PVC2-THY1+ Intermediate',
    '9': 'PVC3-MYH11⁺ Mature ',
    '10': 'Endothelial Cells'
}

# Count cells in each cluster after filtering
cluster_counts = adata.obs['leiden'].value_counts().sort_index()

# Generate labels with formatted counts
new_labels = {
    cluster: f"{custom_names.get(cluster, cluster)} (n={count:,})"
    for cluster, count in cluster_counts.items()
}

# Assign new labels
adata.obs['leiden_named'] = adata.obs['leiden'].map(new_labels)

# Set global font sizes
plt.rcParams.update({
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'legend.fontsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8
})

# Plot the UMAP
fig = sc.pl.umap(
    adata,
    color='leiden_named',
    legend_loc='right margin',
    title='Stroma and Epithelia clustering',
    frameon=False,
    size=1,
    return_fig=True,
    show=False
)

# Set figure size and save
fig.set_size_inches(8.27 / 3.5, 11.69 / 4.2)
plt.show()
fig.savefig("umap_stroma_epithelia_named.pdf", bbox_inches='tight', dpi=75)
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
ax.set_xticklabels([f"SC{int(c)+1}" for c in cluster_paper_percent.index], rotation=0)
ax.set_ylabel("Percentage of Cells", fontsize=10)
ax.set_xlabel("Cluster", fontsize=10)
ax.set_title("Paper Composition per Stroma Cluster", fontsize=11)
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
output_path = "STROMA_cluster_by_paper_percent_barplot_halfA4.pdf"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()

print(f"Saved bar plot to: {output_path}")

####STROMA-IMMUNE ANALYSIS
#CALCULATION OF EACH STROMA CLUSTER SAMPLE PERCENTAGE
# Convert 'leiden' and 'P_P' to string format for consistency

adata.obs['leiden'] = adata.obs['leiden'].astype(str)
adata.obs['P_P'] = adata.obs['P_P'].astype(str)

# Count total number of cells per P_P
pp_counts = adata.obs['P_P'].value_counts()

# Count occurrences of each 'leiden' within each 'P_P'
pp_leiden_counts = adata.obs.groupby(['P_P', 'leiden']).size().unstack(fill_value=0)

# Convert counts to percentages
pp_leiden_percentage = pp_leiden_counts.div(pp_counts, axis=0) * 100

# Save the results to an Excel file
pp_leiden_percentage.to_excel("P_P_leiden_percentage.xlsx")

# Print the first few rows for quick verification
print("Percentage of cells in each Leiden cluster per P_P group:")
print(pp_leiden_percentage.head())
adata.write('stroma_fine_clustering.h5ad')


#Upload myeloid clustering
adata=sc.read_h5ad(r"path\Myeloid_clustering_LA.h5ad")

# Subset to Luminal A
adata_LumA = adata[adata.obs['C.TYPE'] == 'Luminal A'].copy()
adata_LumA.obs['P_P'] = adata_LumA.obs['P_P'].astype(str)


# Set font sizes
plt.rcParams.update({
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8
})

# Superscript formatter for p-value
def to_superscript(n):
    normal = "0123456789-"
    superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"
    return ''.join(superscript[normal.index(c)] for c in str(n))

# Prepare data (use adata_LumA)
adata_LumA.obs['P_P'] = adata_LumA.obs['P_P'].astype(str)
adata_LumA.obs['leiden'] = adata_LumA.obs['leiden'].astype(str)

stroma_file = r"path\P_P_leiden_percentage.xlsx"
stroma_df = pd.read_excel(stroma_file)
stroma_df['P_P'] = stroma_df['P_P'].astype(str)
stroma_cluster5 = stroma_df[['P_P', '5']].rename(columns={'5': 'Stroma_Cluster5'})

leiden_pct_df = (
    adata_LumA.obs.groupby(['P_P', 'leiden']).size().unstack(fill_value=0)
    .div(adata_LumA.obs.groupby('P_P').size(), axis=0) * 100
).reset_index()

merged = pd.merge(leiden_pct_df, stroma_cluster5, on='P_P', how='inner')

filtered = merged[merged['Stroma_Cluster5'] <= 50)]

# Plot and Save
def plot_and_save():
    x = filtered['1']  # 1=MC2
    y = filtered['Stroma_Cluster5']

    if x.nunique() > 1 and y.nunique() > 1:
        r, p = pearsonr(x, y)
    else:
        print("⚠️ Not enough variation to compute Pearson correlation.")
        return

    # Format p-value
    if p < 1e-1:
        base = round(p / (10 ** int(np.floor(np.log10(p)))), 1)
        base = int(base) if base == 1.0 else base
        exponent = int(np.floor(np.log10(p)))
        pval_text = f"P = {base} × 10{to_superscript(exponent)}"
    else:
        pval_text = f"P = {p:.1g}"

    # File sizes
    sizes = {
        "MC2_vs_Stroma5_small.pdf": (8.27 / 2.5, 11.69 / 4.2),
        "MC2_vs_Stroma5_large.pdf": (8.27, 11.69)
    }

    # Plot and save
    for filename, figsize in sizes.items():
        fig, ax = plt.subplots(figsize=figsize)
        sns.regplot(x=x, y=y,
                    scatter_kws={'color': 'royalblue', 's': 30},
                    line_kws={'color': 'red'},
                    ax=ax)

        ax.set_xlabel("MC2 Percentage")
        ax.set_ylabel("CAF_5 Percentage")

        full_title = f"MC2 Correlation with CAF_5 Abundance\n r = {r:.2f}, {pval_text}"
        ax.set_title(full_title, fontsize=9, pad=10)

        ax.grid(False)
        plt.tight_layout()
        plt.savefig(filename)
        if "large" in filename:
            plt.show()
        else:
            plt.close()

# Run
plot_and_save()

##ligand receptor stroma cluster 5 and MC2- luminal A##

# Load the stroma and T cell datasets
stroma = sc.read_h5ad(r"path\stroma_fine_clustering.h5ad")
# Subset to specific clusters
stroma = stroma[stroma.obs['leiden'] == '5'].copy()
# Make sure leiden is a string
adata_LumA.obs['leiden'] = adata_LumA.obs['leiden'].astype(str)

# Subset to leiden cluster 1 =MC2
adata_LumA_cluster2 = adata_LumA[adata_LumA.obs['leiden'] == '1'].copy()#leiden['1']=MC2

# Assign unified TYPE labels
stroma.obs['TYPE'] = 'stroma_cluster_5'
adata_LumA_cluster2.obs['TYPE'] = 'adata_LumA_cluster1'

# Merge the two AnnData objects
adata_combined = stroma.concatenate(adata_LumA_cluster2, batch_key=None)
adata_combined.obs['TYPE'] = adata_combined.obs['TYPE'].astype('category')

# Run LIANA's rank aggregation method
li.method.rank_aggregate(
    adata=adata_combined,
    groupby='TYPE',              # Group cells by this column
    resource_name='consensus',   # Use consensus across LR databases
    expr_prop=0.1,               # Minimum expression in fraction of cells
    min_cells=5,                 # Minimum cells per group
    key_added='liana_res',       # Results stored in adata.uns
    use_raw=False,               # Use .X (not .raw)
    n_perms=1000,                # Permutations for null distribution
    seed=1337,                   # For reproducibility
    verbose=True
)

# Step 6: Export results
df = adata_combined.uns['liana_res']
df.to_csv("liana_stroma5_vs_my_2.csv", index=False)


##all ligand receptor interactions were filtered to include only pvalues<0.05
##All the significant interaction were filtered to include biological relavent interactions
##summerized in this excel file:stroma_5_mac_2_l_r.xlsx

#Ploting results
#Load Excel ===
file_path = r"path\stroma_5_mac_2_l_r.xlsx"
df = pd.read_excel(file_path)

# === Add origin tags ===
def add_origin(row, molecule_type):
    if molecule_type == "ligand":
        if row["source"] == "stroma_cluster_5":
            return f"{row['ligand_complex']} (S)"
        elif row["source"] == "adata_LumA_cluster1":
            return f"{row['ligand_complex']} (M)"
    elif molecule_type == "receptor":
        if row["target"] == "stroma_cluster_5":
            return f"{row['receptor_complex']} (S)"
        elif row["target"] == "adata_LumA_cluster1":
            return f"{row['receptor_complex']} (M)"
    return None

df["ligand_tagged"] = df.apply(lambda row: add_origin(row, "ligand"), axis=1)
df["receptor_tagged"] = df.apply(lambda row: add_origin(row, "receptor"), axis=1)

# === Create interaction label ===
df["interaction"] = df["ligand_tagged"] + " → " + df["receptor_tagged"]

# === Order Y-axis ===
interaction_order = (
    df.groupby("interaction")["lr_logfc"]
    .mean()
    .sort_values()
    .index
)

# === Plot ===
fig, ax = plt.subplots(figsize=(8.27/2, 11.69 / 4.2))  # ≈ 2.36 x 2.78 inches

sns.scatterplot(
    data=df,
    x="lr_logfc",
    y=pd.Categorical(df["interaction"], categories=interaction_order, ordered=True),
    size="magnitude_rank",
    hue="lr_logfc",
    palette="coolwarm",
    sizes=(10, 80),
    edgecolor="black",
    ax=ax,
    legend=False
)

# === Style ===
ax.set_xlabel("Log Fold Change", fontsize=8)
ax.set_ylabel("", fontsize=8)
ax.set_title("Ligand–Receptor Interactions", fontsize=9)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

# === Custom colorbar (LogFC) ===
norm = mpl.colors.Normalize(vmin=df["lr_logfc"].min(), vmax=df["lr_logfc"].max())
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.01, shrink=0.6)
cbar.set_label("Log Fold Change", fontsize=8)
cbar.ax.tick_params(labelsize=8)

# === Move colorbar higher ===
cbar.ax.set_position([
    cbar.ax.get_position().x0,
    cbar.ax.get_position().y0 + 0.08,
    cbar.ax.get_position().width,
    cbar.ax.get_position().height * 0.9
])

# === Custom size legend (inverted size + custom label) ===
min_val = df["magnitude_rank"].min()
max_val = df["magnitude_rank"].max()

for size_val in [min_val, max_val]:
    size_scaled = (1 - (size_val - min_val) / (max_val - min_val))
    label = "<10\u207B\u2075" if size_val == min_val else f'{size_val:.3f}'
    plt.scatter([], [], 
                s=(10 + 70 * size_scaled),
                label=label,
                edgecolors='black',
                facecolors='gray')

size_legend = plt.legend(
    title="Magnitude Rank",
    bbox_to_anchor=(1.05, 0.01),
    loc='lower left',
    fontsize=8,
    title_fontsize=8,
    borderaxespad=0.
)
plt.gca().add_artist(size_legend)

# === Tight layout + Save ===
plt.tight_layout()
output_path = r"path\ligand_receptor_plot.pdf"
plt.savefig(output_path, format="pdf", bbox_inches='tight') 
# === Show ===
plt.show()


# In[ ]:


####stroma correlation with immune infiltration


# In[ ]:


# === Set font sizes ===
plt.rcParams.update({
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8
})

# === Helper function for superscript p-value ===
def to_superscript(n):
    normal = "0123456789-"
    superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"
    return ''.join(superscript[normal.index(c)] for c in str(n))

# === Subset Luminal A samples ===
adata_LumA = adata[adata.obs['C.TYPE'] == 'Luminal A'].copy()
adata_LumA.obs['P_P'] = adata_LumA.obs['P_P'].astype(str) #P_P=Sample

# === Load stroma cluster 5 percentages ===
stroma_file = r"path\P_P_leiden_percentage.xlsx"
stroma_df = pd.read_excel(stroma_file)
stroma_df['P_P'] = stroma_df['P_P'].astype(str)
stroma_cluster5 = stroma_df[['P_P', '5']].rename(columns={'5': 'Stroma_Cluster5'})

# Path to the Excel files contaning the high/low immune infiltration samples
high_file = r"path\high.xlsx"
low_file = r"path\low.xlsx"

# Load sample IDs from each Excel file
# Assuming the IDs are in the first column, no header
high_df = pd.read_excel(high_file, header=None)
low_df = pd.read_excel(low_file, header=None)

# Convert to Python lists of strings
high = high_df.iloc[:, 0].astype(str).tolist()
low = low_df.iloc[:, 0].astype(str).tolist()

print("High group sample count:", len(high))
print("Low group sample count:", len(low))
print("First few high samples:", high[:5])
print("First few low samples:", low[:5])

# === Drop missing values before checking valid P_Ps (P_P=Sample)===
stroma_cluster5 = stroma_cluster5.dropna(subset=['Stroma_Cluster5'])

# === Filter valid P_Ps and clean group lists ===
valid_pps = set(stroma_cluster5['P_P'].unique())
high_filtered = [pp for pp in high if pp in valid_pps]
low_filtered = [pp for pp in low if pp in valid_pps]

# === Assign immune groups ===
def assign_group(pp):
    if pp in high_filtered:
        return 'High'
    elif pp in low_filtered:
        return 'Low'
    else:
        return 'Unknown'

stroma_cluster5['Immune_Group'] = stroma_cluster5['P_P'].apply(assign_group)

# === Filter to just High and Low ===
plot_df = stroma_cluster5[stroma_cluster5['Immune_Group'].isin(['High', 'Low'])].copy()

# === Mann–Whitney U test ===
high_vals = plot_df[plot_df['Immune_Group'] == 'High']['Stroma_Cluster5']
low_vals = plot_df[plot_df['Immune_Group'] == 'Low']['Stroma_Cluster5']

if len(high_vals) == 0 or len(low_vals) == 0:
    raise ValueError(f"Empty group: High={len(high_vals)}, Low={len(low_vals)}")

u_stat, p_val = mannwhitneyu(high_vals, low_vals, alternative='two-sided')

# === Format p-value ===
if p_val < 1e-1:
    base = round(p_val / (10 ** int(np.floor(np.log10(p_val)))), 1)
    base = int(base) if base == 1.0 else base
    exponent = int(np.floor(np.log10(p_val)))
    pval_text = f"P = {base} × 10{to_superscript(exponent)}"
else:
    pval_text = f"P = {p_val:.1g}"

# === Define plot colors ===
colors = {'High': '#D8BFD8', 'Low': '#ADD8E6'}

# === Plot small version ===
fig, ax = plt.subplots(figsize=(8.27 / 2.5, 11.69 / 4.2))
sns.boxplot(data=plot_df, x='Immune_Group', y='Stroma_Cluster5', palette=colors, ax=ax)
sns.stripplot(data=plot_df, x='Immune_Group', y='Stroma_Cluster5',
              color='black', size=3, jitter=True, dodge=False, ax=ax)
ax.set_title(f"CAF_5 Abundance:\n High vs Low Immune Infiltration\n{pval_text}")
ax.set_xlabel("Immune Infiltration Group")
ax.set_ylabel("CAF_5 Percentage")
if ax.get_legend():
    ax.legend_.remove()
plt.tight_layout()
plt.savefig("cluster5_small.pdf")
plt.close()

# === Plot large version ===
fig, ax = plt.subplots(figsize=(8.27, 11.69))
sns.boxplot(data=plot_df, x='Immune_Group', y='Stroma_Cluster5', palette=colors, ax=ax)
sns.stripplot(data=plot_df, x='Immune_Group', y='Stroma_Cluster5',
              color='black', size=3, jitter=True, dodge=False, ax=ax)
ax.set_title(f"CAF_5 Abundance: High vs Low \n Immune Infiltration\n{pval_text}")
ax.set_xlabel("Immune Infiltration Group")
ax.set_ylabel("CAF_5 Percentage")
if ax.get_legend():
    ax.legend_.remove()
plt.tight_layout()
plt.savefig("cluster5_large.pdf")
plt.show()


