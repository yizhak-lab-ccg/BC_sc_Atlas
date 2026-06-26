import scanpy as sc
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

###Upload dataset
adata=sc.read_h5ad("path\Bassez_cohort.h5ad")
adata.obs['cluster_name'].unique()
adata

# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata,target_sum=1e4)

# Logarithmize the data
sc.pp.log1p(adata)

# Ensure 'cluster_name' column exists
if 'cluster_name' not in adata.obs.columns:
    raise ValueError("'cluster_name' column not found in adata.obs!")

# Subset the AnnData object to include only T cells
adata_t_cells = adata[adata.obs['cluster_name'] == 'T_cell'].copy()

# Print shape of the new AnnData object
print(f"Subset AnnData object created with shape: {adata_t_cells.shape}")

# Optionally, save the new object
# adata_t_cells.write("adata_T_cells.h5ad")

sc.pp.highly_variable_genes(adata_t_cells, n_top_genes=2000, batch_key="patient_id")

sc.pl.highly_variable_genes(adata_t_cells)

sc.tl.pca(adata_t_cells)

sc.pp.neighbors(adata_t_cells)

sc.tl.umap(adata_t_cells)

sc.pl.umap(adata_t_cells, color=['patient_id',"cellType"])

sc.tl.leiden(adata_t_cells,resolution=0.5)

sc.pl.umap(adata_t_cells, color=['patient_id',"leiden"])

sc.tl.rank_genes_groups(adata_t_cells, groupby="leiden", method="wilcoxon")

sc.pl.rank_genes_groups_dotplot(adata_t_cells, groupby="leiden",standard_scale="var", n_genes=5)import pandas as pd

# Create a Pandas Excel writer to save to multiple sheets
with pd.ExcelWriter('DEG_adata_t_cells.xlsx', engine='xlsxwriter') as writer:
    for group in adata_t_cells.uns['rank_genes_groups']['names'].dtype.names:
        # Create a dataframe for each group containing gene names, p-values, log fold change, and adjusted p-values
        genes_df = pd.DataFrame({
            'Gene': adata_t_cells.uns['rank_genes_groups']['names'][group],
            'Log2 Fold Change': adata_t_cells.uns['rank_genes_groups']['logfoldchanges'][group],
            'p-value': adata_t_cells.uns['rank_genes_groups']['pvals'][group],
            'Adjusted p-value': adata_t_cells.uns['rank_genes_groups']['pvals_adj'][group]
        })

        # Filter for adjusted p-value < 0.05 and log2 fold change > 0.5
        genes_df = genes_df[
            (genes_df['Adjusted p-value'] < 0.05) & (genes_df['Log2 Fold Change'] > 0.5)
        ]

        # Remove genes with specified prefixes
        genes_df = genes_df[~genes_df['Gene'].str.startswith(genes_to_remove)]

        # Save each group's DEGs to a separate sheet named after the group
        genes_df.to_excel(writer, sheet_name=f'Group_{group}', index=False)

# Print completion message
print("Filtered DEG results saved to 'DEG_adata_t_cells.xlsx'.")

adata_t_cells.obs['expansion'].unique()

# Subset the AnnData object to only include cells where BC_type is TNBC
adata_t_cells_tnbc = adata_t_cells[adata_t_cells.obs['BC_type'] == 'TNBC'].copy()

# Print the shape to confirm the subset
print(f"Subset AnnData object shape: {adata_t_cells_tnbc.shape}")

# Subset for Leiden clusters 2, 3, and 4
adata_t_cells_tnbc_leiden234 = adata_t_cells_tnbc[adata_t_cells_tnbc.obs['leiden'].isin(['2', '3', '4'])].copy()

genes = ['LAG3', 'PDCD1']
for gene in genes:
    if gene not in adata_t_cells_tnbc_leiden234.var_names:
        raise ValueError(f"The gene {gene} is not found in the dataset.")

# Convert to dense if sparse
X = adata_t_cells_tnbc_leiden234[:, genes].X
if not isinstance(X, np.ndarray):
    X = X.toarray()

# Create a DataFrame with expression > 0.5 binary matrix
expr_binary = (X > 0.5).astype(int)
sample_info = adata_t_cells_tnbc_leiden234.obs[['sample', 'response']].copy()
expr_df = pd.DataFrame(expr_binary, columns=genes, index=sample_info.index)
expr_df = pd.concat([expr_df, sample_info], axis=1)

# Group by sample, compute fraction per gene, then average
fraction_by_sample = (
    expr_df.groupby('sample')[genes]
    .mean()
    .assign(Response=expr_df.groupby('sample')['response'].first())
)
fraction_by_sample['Mean Fraction Expression'] = fraction_by_sample[genes].mean(axis=1)

# Prepare for plotting
plot_df = fraction_by_sample.reset_index()
plot_df['Response'] = pd.Categorical(plot_df['Response'], categories=['E', 'NE'], ordered=True)
plot_df = plot_df.dropna()

# Mann-Whitney U test
group_e = plot_df[plot_df['Response'] == 'E']['Mean Fraction Expression']
group_ne = plot_df[plot_df['Response'] == 'NE']['Mean Fraction Expression']
stat, p_value = mannwhitneyu(group_e, group_ne, alternative='two-sided')

# Plot
plt.figure(figsize=(8, 6))
sns.boxplot(data=plot_df, x='Response', y='Mean Fraction Expression', palette='Set2', showfliers=False)
sns.stripplot(data=plot_df, x='Response', y='Mean Fraction Expression', color='black', size=4, jitter=True, alpha=0.6)

max_y = plot_df['Mean Fraction Expression'].max()
plt.text(0.5, max_y + 0.01, f"p-value: {p_value:.4e}", ha='center', va='bottom', fontsize=12)

plt.title("LAG3 & PDCD1 Expression by Response Group")
plt.xlabel("Response Group")
plt.ylabel("Mean Fraction Expression")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()

