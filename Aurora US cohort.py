import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
import numpy as np
from lifelines.statistics import logrank_test

# Define the file path
file_path = r"path to raw data GSE209998.txt"

# Load the data
# Assumes genes are rows and cells are columns. If not, set transpose=True.
df = pd.read_csv(file_path, sep='\t', index_col=0)

# Create AnnData
adata = sc.AnnData(df.T)  # transpose to have cells as rows and genes as columns

# Check
print(adata)

#########Metadata1#############
# Step 1: Load metadata
metadata_path = r"path\metadata_1.xlsx"
metadata = pd.read_excel(metadata_path)

# Step 2: Remove rows where RNA Seq FreezeSet 123 == 'Failed'
metadata = metadata[metadata['RNA Seq FreezeSet 123'] != 'Failed']

# Step 3: Create 'sample_id' from adata.obs_names (first 4 parts)
adata.obs['sample_id'] = adata.obs_names.str.split('-').str[:4].str.join('-')

# Step 4: Find unmatched sample_ids
adata_sample_ids = set(adata.obs['sample_id'])
metadata_sample_ids = set(metadata['BCR Sample barcode'])
unmatched_sample_ids = adata_sample_ids - metadata_sample_ids

# Step 5: Filter out unmatched rows from adata
adata = adata[~adata.obs['sample_id'].isin(unmatched_sample_ids)].copy()

# Step 6: Drop duplicates in metadata and set index
metadata = metadata.drop_duplicates(subset='BCR Sample barcode')
metadata = metadata.set_index('BCR Sample barcode')

# Step 7: Map metadata by sample_id
meta_dict = metadata.to_dict(orient='index')
aligned_metadata = adata.obs['sample_id'].map(meta_dict)
aligned_metadata_df = pd.DataFrame(list(aligned_metadata), index=adata.obs.index)

# Step 8: Add metadata to adata.obs
adata.obs = pd.concat([adata.obs, aligned_metadata_df], axis=1)

# Final check
print("Final shape of adata:", adata.shape)
print("Metadata columns added:", aligned_metadata_df.columns.tolist())


########Metadata2############

# Step 1: Load the second metadata file
metadata2_path = r"path\metadata_2.xlsx"
metadata2 = pd.read_excel(metadata2_path)

# Step 2: Ensure types are strings for matching
adata.obs['Patient'] = adata.obs['Patient'].astype(str)
metadata2['Sample'] = metadata2['Sample'].astype(str)

# Step 3: Drop duplicates in metadata if needed
metadata2 = metadata2.drop_duplicates(subset='Sample')

# Step 4: Set index to 'Sample' in metadata2
metadata2 = metadata2.set_index('Sample')

# Step 5: Map metadata to adata.obs by matching 'Patient'
meta2_dict = metadata2.to_dict(orient='index')
aligned_meta2 = adata.obs['Patient'].map(meta2_dict)
aligned_meta2_df = pd.DataFrame(list(aligned_meta2), index=adata.obs.index)

# Step 6: Concatenate the aligned metadata
adata.obs = pd.concat([adata.obs, aligned_meta2_df], axis=1)

# Check result
print("Updated adata.obs with metadata_2 info:")
print(adata.obs.head())


##########Normalization############
# Step 1: Save raw counts to .layers['counts']
adata.layers["counts"] = adata.X.copy()

# Step 2: Compute CPM (Counts Per Million)
cell_sums = adata.X.sum(axis=1)
cpm = (adata.X / cell_sums[:, None]) * 1e6

# Step 3: Compute log2(CPM + 1)
log_cpm = np.log2(cpm + 1)

# Step 4: Overwrite .X with log-transformed CPM
adata.X = log_cpm

print("log2(CPM + 1) transformation complete and stored in .X")
print("Raw counts preserved in adata.layers['counts']")


# Step 1: Remove exact duplicate columns
adata.obs = adata.obs.loc[:, ~adata.obs.columns.duplicated()]
# Step 2: Subset TNBC cells
adata_tnbc = adata[adata.obs['Metastasis original receptors'] == 'TNBC'].copy()

# Check result
print("Subset complete. adata_tnbc shape: {adata_tnbc.shape}")


#####TNBC Samples#####


# Step 1: Extract data from adata_tnbc.obs
df = adata_tnbc.obs.copy()

# Step 2: Define survival time and event status
# Make sure they're numeric and drop NAs
survival_time = pd.to_numeric(df['Metastasis OS months'], errors='coerce')
event_status = pd.to_numeric(df['Metastasis OS Event'], errors='coerce')

km_data = pd.DataFrame({
    'time': survival_time,
    'event': event_status
}).dropna()

# Step 3: Fit and plot Kaplan–Meier curve
kmf = KaplanMeierFitter()
kmf.fit(durations=km_data['time'], event_observed=km_data['event'])

# Plot
plt.figure(figsize=(6, 4))
kmf.plot()
plt.title('Kaplan–Meier: Metastasis OS (TNBC)')
plt.xlabel('Time (months)')
plt.ylabel('Survival Probability')
plt.grid(True)
plt.tight_layout()
plt.show()


################

##km curve
# Step 1: Label cells as High if > 5.2(median), else Low
adata_tnbc.obs['expr_group'] = np.where(
    adata_tnbc.obs['mean_signature_log2CPM'] > 5.2, 'High', 'Low'
)

# Step 2: Group per patient
group_df = (
    adata_tnbc.obs.groupby('Patient')['expr_group']
    .apply(lambda x: 'High' if (x == 'High').mean() > 0.5 else 'Low')
    .reset_index()
    .rename(columns={'expr_group': 'Patient_group'})
)

# Step 3: Merge with survival data
surv_cols = ['Patient', 'Metastasis OS months', 'Metastasis OS Event']
surv_data = adata_tnbc.obs[surv_cols].drop_duplicates('Patient')
merged = pd.merge(group_df, surv_data, on='Patient')
merged = merged.dropna(subset=['Metastasis OS months', 'Metastasis OS Event', 'Patient_group'])

# Step 4: Split into High and Low groups
high = merged[merged['Patient_group'] == 'High']
low = merged[merged['Patient_group'] == 'Low']

# Step 5: Log-rank test
results = logrank_test(
    high['Metastasis OS months'], low['Metastasis OS months'],
    event_observed_A=high['Metastasis OS Event'],
    event_observed_B=low['Metastasis OS Event']
)
p_val = results.p_value

# Step 6: Plot with custom colors and font sizes
plt.figure(figsize=(6, 4))
kmf = KaplanMeierFitter()

# Define color map
color_map = {'High': '#5DA5DA', 'Low': '#FAA43A'}  # blue, yellow-orange

for label, group in merged.groupby("Patient_group"):
    kmf.fit(group["Metastasis OS months"], event_observed=group["Metastasis OS Event"], label=label)
    kmf.plot(ci_show=True, ci_alpha=0.2, linewidth=1.5, color=color_map[label])

# Formatting
plt.title(f"Kaplan–Meier: High vs Low Signature Expression\np = {p_val:.3f}", fontsize=9)
plt.xlabel("Time (months)", fontsize=8)
plt.ylabel("Survival Probability", fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.legend(title=None, fontsize=8)
plt.grid(False)
plt.tight_layout()
plt.show()

