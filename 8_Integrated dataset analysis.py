import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import scanpy as sc
import os
import liana as li
import networkx as nx
import numpy as np
from matplotlib.patches import PathPatch, Patch
from matplotlib.path import Path

#Integrating the IMMUNE/TUMOR/STROMA datasets to one integrated dataset
#P_P=sample
#C.TYPE=BC Subtype
import pandas as pd

file_path = r"path/Tumor.h5ad"
adata1 = sc.read(file_path)

file_path = r"path/Stroma.h5ad"
adata2 = sc.read(file_path)

file_path = r"path/Immune.h5ad"
adata3 = sc.read(file_path)

# Extract relevant information and add TYPE column
adata1_df = pd.DataFrame(adata1.obs['TYPE']).assign(TYPE='TUMOR')
adata2_df = pd.DataFrame(adata2.obs['TYPE']).assign(TYPE='STROMA')
adata3_df = pd.DataFrame(adata3.obs['TYPE']).assign(TYPE='IMMUNE')

# Concatenate all data
combined_df = pd.concat([adata1_df, adata2_df, adata3_df], axis=0)

# Save to CSV
combined_df.to_csv('combined_cells.csv', index_label='cell')

####Immune infiltration analysis####
#C.TYPE=BC subtype
#P_P=sample
#BIOPSY=Biopsy location
# Path to the .h5ad file
file_path = r"path\combined_df.h5ad"

# Load the AnnData object
adata = sc.read(file_path)
print("Data loaded successfully!")

# Modify the C.TYPE column to rename specific values
adata.obs['C.TYPE'] = adata.obs['C.TYPE'].replace({
    'LUMINALA': 'Luminal A',
    'LUMINALB': 'Luminal B',
    'TN': 'TNBC'
})

# Verify the renaming was successful
print("Updated C.TYPE values:")
print(adata.obs['C.TYPE'].unique())

# Ensure the 'C.TYPE' and 'TYPE' columns exist in adata.obs
if 'C.TYPE' in adata.obs.columns and 'TYPE' in adata.obs.columns:
    # Total number of cells for each C.TYPE
    total_cells_per_ctype = adata.obs['C.TYPE'].value_counts()

    # Total number of immune cells (where TYPE == 'IMMUNE') for each C.TYPE
    immune_cells_per_ctype = (
        adata.obs[adata.obs['TYPE'] == 'IMMUNE']['C.TYPE'].value_counts()
    )

    # Combine both counts into a DataFrame
    percentage_df = pd.DataFrame({
        'Total_Cells': total_cells_per_ctype,
        'Immune_Cells': immune_cells_per_ctype
    }).fillna(0)  # Replace NaN with 0 for missing immune counts

    # Calculate the percentage of immune cells
    percentage_df['Immune_Percentage'] = (
        (percentage_df['Immune_Cells'] / percentage_df['Total_Cells']) * 100
    )

    # === Exclude C.TYPE categories with only immune cells (100% immune) ===
    percentage_df = percentage_df[percentage_df['Immune_Percentage'] < 100]

    # Display the resulting DataFrame
    print("Immune cell percentages per C.TYPE (excluding all-immune samples):")
    print(percentage_df)

    # Plot the immune cell percentage for each C.TYPE
    plt.figure(figsize=(10, 6))
    plt.bar(percentage_df.index, percentage_df['Immune_Percentage'])
    plt.title("Percentage of Immune Cells in Each BC subtype")
    plt.xlabel("BC subtype")
    plt.ylabel("Percentage of Immune Cells (%)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()

else:
    print("Columns 'C.TYPE' or 'TYPE' not found in adata.obs. Please check your data.")


# Path to the .h5ad file
file_path = r"path\combined_df.h5ad"

# Load the AnnData object
adata = sc.read(file_path)
print("Data loaded successfully!")

# Modify the C.TYPE column to rename specific values
adata.obs['C.TYPE'] = adata.obs['C.TYPE'].replace({
    'LUMINALA': 'Luminal A',
    'LUMINALB': 'Luminal B',
    'IS': 'DCIS',
    'TN': 'TNBC',
    'NA': 'UNKNOWN',
    None: 'UNKNOWN'  # Replace NaN values
})

# Verify the renaming was successful
print("Updated C.TYPE values:")
print(adata.obs['C.TYPE'].unique())

# Ensure the 'C.TYPE' and 'TYPE' columns exist in adata.obs
if 'C.TYPE' in adata.obs.columns and 'TYPE' in adata.obs.columns:
    # Total number of cells for each C.TYPE
    total_cells_per_ctype = adata.obs['C.TYPE'].value_counts()

    # Total number of immune cells (where TYPE == 'IMMUNE') for each C.TYPE
    immune_cells_per_ctype = (
        adata.obs[adata.obs['TYPE'] == 'IMMUNE']['C.TYPE'].value_counts()
    )

    # Combine both counts into a DataFrame
    percentage_df = pd.DataFrame({
        'Total_Cells': total_cells_per_ctype,
        'Immune_Cells': immune_cells_per_ctype
    }).fillna(0)  # Replace NaN with 0 for missing immune counts

    # Calculate the percentage of immune cells
    percentage_df['Immune_Percentage'] = (
        (percentage_df['Immune_Cells'] / percentage_df['Total_Cells']) * 100
    )

    # Display the resulting DataFrame
    print("Immune cell percentages per C.TYPE:")
    print(percentage_df)

    # Plot the immune cell percentage for each C.TYPE
    plt.figure(figsize=(10, 6))
    plt.bar(percentage_df.index, percentage_df['Immune_Percentage'])
    plt.title("Percentage of Immune Cells in Each BC subtype")
    plt.xlabel("BC subtype")
    plt.ylabel("Percentage of Immune Cells (%)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()

else:
    print("Columns 'C.TYPE' or 'TYPE' not found in adata.obs. Please check your data.")


# Ensure 'C.TYPE' and 'P_P' columns are strings
adata.obs['C.TYPE'] = adata.obs['C.TYPE'].astype(str)
adata.obs['P_P'] = adata.obs['P_P'].astype(str)
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

adata.obs['BIOPSY_GROUP'] = adata.obs['BIOPSY'].apply(define_biopsy_group)
adata = adata[adata.obs['BIOPSY_GROUP'].notna()]
# Step 1: Remove 'unknown', 'DCIS' values from C.TYPE and rename to adata2
adata2 = adata[~adata.obs['C.TYPE'].isin(['unknown', 'DCIS', "UNKNOWN"])]

# Group data to calculate total cells and immune cells per P_P and C.TYPE
counts = (
    adata2.obs
    .groupby(['P_P', 'C.TYPE'], observed=True)
    .size()
    .reset_index(name='Total_Cells')
)

immune_counts = (
    adata2.obs[adata2.obs['TYPE'] == 'IMMUNE']
    .groupby(['P_P', 'C.TYPE'], observed=True)
    .size()
    .reset_index(name='Immune_Cells')
)

# Merge total and immune counts
merged_counts = pd.merge(counts, immune_counts, on=['P_P', 'C.TYPE'], how='left').fillna(0)

# Calculate immune cell percentages
merged_counts['Immune_Percentage'] = (
    merged_counts['Immune_Cells'] / merged_counts['Total_Cells']
) * 100

# Exclude P_P × C.TYPE groups with only immune cells (100%) ===
merged_counts = merged_counts[merged_counts['Immune_Percentage'] < 100]

# Step 2: Order C.TYPE by median Immune_Percentage
ctype_median_order = (
    merged_counts.groupby('C.TYPE')['Immune_Percentage']
    .median()
    .sort_values()
    .index
)

# Update C.TYPE as a categorical type with ordered categories
merged_counts['C.TYPE'] = pd.Categorical(
    merged_counts['C.TYPE'], categories=ctype_median_order, ordered=True
)

# Step 3: Calculate p-values for pairwise comparisons
def calculate_p_values(data, group_column, value_column):
    unique_groups = data[group_column].unique()
    p_values = {}

    for i, group1 in enumerate(unique_groups):
        for group2 in unique_groups[i + 1:]:
            # Extract data for the two groups
            values1 = data[data[group_column] == group1][value_column]
            values2 = data[data[group_column] == group2][value_column]
            
            # Perform Mann-Whitney U test
            stat, p_value = mannwhitneyu(values1, values2, alternative='two-sided')
            if p_value < 0.05:  # Only store significant p-values
                p_values[(group1, group2)] = p_value

    return p_values

# Calculate p-values for immune percentages between C.TYPE groups
p_values = calculate_p_values(merged_counts, 'C.TYPE', 'Immune_Percentage')

# Print significant p-values
print("Significant p-values (p < 0.05) between C.TYPE groups:")
for comparison, p_value in p_values.items():
    print(f"{comparison[0]} vs {comparison[1]}: p-value = {p_value:.4e}")

# Step 4: Plot the boxplot for immune cell percentages grouped by ordered C.TYPE
plt.figure(figsize=(12, 8))

# Use a color palette to differentiate C.TYPE groups
palette = sns.color_palette("Set2", n_colors=len(merged_counts['C.TYPE'].unique()))

# Boxplot with custom colors and ordered categories
sns.boxplot(
    data=merged_counts, 
    x='C.TYPE', 
    y='Immune_Percentage',
    palette=palette
)

# Stripplot for individual points
sns.stripplot(
    data=merged_counts, 
    x='C.TYPE', 
    y='Immune_Percentage', 
    color='black', 
    alpha=0.5, 
    jitter=True
)

# Annotate only significant p-values above the lines
y_offset = 3  # Adjust vertical spacing
line_spacing = 2  # Spacing between lines
max_y = merged_counts['Immune_Percentage'].max()

for (group1, group2), p_value in p_values.items():
    x1 = ctype_median_order.tolist().index(group1)
    x2 = ctype_median_order.tolist().index(group2)
    y = max_y + y_offset

    # Draw the connecting line
    plt.plot([x1, x1, x2, x2], [y, y + line_spacing, y + line_spacing, y], lw=1.5, c='black')
    
    # Place the p-value above the line with extra spacing
    plt.text((x1 + x2) / 2, y + line_spacing + 0.5, f"p={p_value:.4e}", 
             ha='center', va='bottom', fontsize=9)
    
    # Update vertical spacing for next line
    y_offset += line_spacing + 1

plt.title("Distribution of Immune Cell Percentages per BC subtype across samples")
plt.xlabel("BC subtype")
plt.ylabel("Immune Cell Percentage (%)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()


###Primary Breast cancer graph###

# Ensure 'C.TYPE' and 'P_P' columns are strings
adata.obs['C.TYPE'] = adata.obs['C.TYPE'].astype(str)
adata.obs['P_P'] = adata.obs['P_P'].astype(str)

# Step 1: Keep only cells with BIOPSY = BREAST or LN, and remove unwanted C.TYPE values
adata2 = adata[
    (adata.obs['BIOPSY_GROUP'].isin(['BREAST', 'LN'])) &
    (~adata.obs['C.TYPE'].isin(['unknown', 'DCIS', "UNKNOWN"]))
]

# Group data to calculate total cells and immune cells per P_P and C.TYPE
counts = (
    adata2.obs
    .groupby(['P_P', 'C.TYPE'], observed=True)
    .size()
    .reset_index(name='Total_Cells')
)

immune_counts = (
    adata2.obs[adata2.obs['TYPE'] == 'IMMUNE']
    .groupby(['P_P', 'C.TYPE'], observed=True)
    .size()
    .reset_index(name='Immune_Cells')
)

# Merge total and immune counts
merged_counts = pd.merge(counts, immune_counts, on=['P_P', 'C.TYPE'], how='left').fillna(0)

# Calculate immune cell percentages
merged_counts['Immune_Percentage'] = (
    merged_counts['Immune_Cells'] / merged_counts['Total_Cells']
) * 100

# === Exclude P_P × C.TYPE groups with only immune cells (100%) ===
merged_counts = merged_counts[merged_counts['Immune_Percentage'] < 100]

# Order C.TYPE by median Immune_Percentage
ctype_median_order = (
    merged_counts.groupby('C.TYPE')['Immune_Percentage']
    .median()
    .sort_values()
    .index
)

# Update C.TYPE as a categorical type with ordered categories
merged_counts['C.TYPE'] = pd.Categorical(
    merged_counts['C.TYPE'], categories=ctype_median_order, ordered=True
)

# Calculate p-values for pairwise comparisons
def calculate_p_values(data, group_column, value_column):
    unique_groups = data[group_column].unique()
    p_values = {}

    for i, group1 in enumerate(unique_groups):
        for group2 in unique_groups[i + 1:]:
            values1 = data[data[group_column] == group1][value_column]
            values2 = data[data[group_column] == group2][value_column]
            
            # Perform Mann-Whitney U test
            stat, p_value = mannwhitneyu(values1, values2, alternative='two-sided')
            if p_value < 0.05:  # Only store significant p-values
                p_values[(group1, group2)] = p_value

    return p_values

# Calculate p-values for immune percentages between C.TYPE groups
p_values = calculate_p_values(merged_counts, 'C.TYPE', 'Immune_Percentage')

# Print significant p-values
print("Significant p-values (p < 0.05) between C.TYPE groups (BIOPSY = BREAST or LN):")
for comparison, p_value in p_values.items():
    print(f"{comparison[0]} vs {comparison[1]}: p-value = {p_value:.4e}")

# Plot the boxplot for immune cell percentages grouped by ordered C.TYPE
plt.figure(figsize=(12, 8))

palette = sns.color_palette("Set2", n_colors=len(merged_counts['C.TYPE'].unique()))

sns.boxplot(
    data=merged_counts, 
    x='C.TYPE', 
    y='Immune_Percentage',
    palette=palette
)

sns.stripplot(
    data=merged_counts, 
    x='C.TYPE', 
    y='Immune_Percentage', 
    color='black', 
    alpha=0.5, 
    jitter=True
)

# Annotate only significant p-values above the lines
y_offset = 3
line_spacing = 2
max_y = merged_counts['Immune_Percentage'].max()

for (group1, group2), p_value in p_values.items():
    x1 = ctype_median_order.tolist().index(group1)
    x2 = ctype_median_order.tolist().index(group2)
    y = max_y + y_offset

    plt.plot([x1, x1, x2, x2], [y, y + line_spacing, y + line_spacing, y], lw=1.5, c='black')
    plt.text((x1 + x2) / 2, y + line_spacing + 0.5, f"p={p_value:.4e}", 
             ha='center', va='bottom', fontsize=9)
    
    y_offset += line_spacing + 1

plt.title("Distribution of Immune Cell Percentages per BC subtype (BIOPSY = BREAST or LN)")
plt.xlabel("BC subtype")
plt.ylabel("Immune Cell Percentage (%)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()

### Metastatic Breast cancer graph ###


# Ensure 'C.TYPE' and 'P_P' columns are strings
adata.obs['C.TYPE'] = adata.obs['C.TYPE'].astype(str)
adata.obs['P_P'] = adata.obs['P_P'].astype(str)

# Keep only cells with BIOPSY_GROUP = MET, and remove unwanted C.TYPE values
adata2 = adata[
    (adata.obs['BIOPSY_GROUP'] == 'MET') &
    (~adata.obs['C.TYPE'].isin(['unknown', 'DCIS',"UNKNOWN"]))
]

# Group data to calculate total cells and immune cells per P_P and C.TYPE
counts = (
    adata2.obs
    .groupby(['P_P', 'C.TYPE'], observed=True)
    .size()
    .reset_index(name='Total_Cells')
)

immune_counts = (
    adata2.obs[adata2.obs['TYPE'] == 'IMMUNE']
    .groupby(['P_P', 'C.TYPE'], observed=True)
    .size()
    .reset_index(name='Immune_Cells')
)

# Merge total and immune counts
merged_counts = pd.merge(counts, immune_counts, on=['P_P', 'C.TYPE'], how='left').fillna(0)

# Calculate immune cell percentages
merged_counts['Immune_Percentage'] = (
    merged_counts['Immune_Cells'] / merged_counts['Total_Cells']
) * 100

# === Exclude P_P × C.TYPE groups with only immune cells (100%) ===
merged_counts = merged_counts[merged_counts['Immune_Percentage'] < 100]

# Step 2: Order C.TYPE by median Immune_Percentage
ctype_median_order = (
    merged_counts.groupby('C.TYPE')['Immune_Percentage']
    .median()
    .sort_values()
    .index
)

# Update C.TYPE as a categorical type with ordered categories
merged_counts['C.TYPE'] = pd.Categorical(
    merged_counts['C.TYPE'], categories=ctype_median_order, ordered=True
)

# Calculate p-values for pairwise comparisons
def calculate_p_values(data, group_column, value_column):
    unique_groups = data[group_column].unique()
    p_values = {}

    for i, group1 in enumerate(unique_groups):
        for group2 in unique_groups[i + 1:]:
            values1 = data[data[group_column] == group1][value_column]
            values2 = data[data[group_column] == group2][value_column]
            
            stat, p_value = mannwhitneyu(values1, values2, alternative='two-sided')
            if p_value < 0.05:  # Only store significant p-values
                p_values[(group1, group2)] = p_value

    return p_values

# Calculate p-values for immune percentages between C.TYPE groups
p_values = calculate_p_values(merged_counts, 'C.TYPE', 'Immune_Percentage')

# Print significant p-values
print("Significant p-values (p < 0.05) between C.TYPE groups (BIOPSY_GROUP = MET):")
for comparison, p_value in p_values.items():
    print(f"{comparison[0]} vs {comparison[1]}: p-value = {p_value:.4e}")

# Plot the boxplot for immune cell percentages grouped by ordered C.TYPE
plt.figure(figsize=(12, 8))

palette = sns.color_palette("Set2", n_colors=len(merged_counts['C.TYPE'].unique()))

sns.boxplot(
    data=merged_counts, 
    x='C.TYPE', 
    y='Immune_Percentage',
    palette=palette
)

sns.stripplot(
    data=merged_counts, 
    x='C.TYPE', 
    y='Immune_Percentage', 
    color='black', 
    alpha=0.5, 
    jitter=True
)

# Annotate only significant p-values above the lines
y_offset = 3
line_spacing = 2
max_y = merged_counts['Immune_Percentage'].max()

for (group1, group2), p_value in p_values.items():
    x1 = ctype_median_order.tolist().index(group1)
    x2 = ctype_median_order.tolist().index(group2)
    y = max_y + y_offset

    plt.plot([x1, x1, x2, x2], [y, y + line_spacing, y + line_spacing, y], lw=1.5, c='black')
    plt.text((x1 + x2) / 2, y + line_spacing + 0.5, f"p={p_value:.4e}", 
             ha='center', va='bottom', fontsize=9)
    
    y_offset += line_spacing + 1

plt.title("Distribution of Immune Cell Percentages per BC subtype (BIOPSY_GROUP = MET)")
plt.xlabel("BC subtype")
plt.ylabel("Immune Cell Percentage (%)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()


###Ligand_receptor analysis###
#Luminal A ligand-receptor analysis

all=sc.read_h5ad(r'path/integrated dataset.h5ad')
all

# Replace values in 'C.TYPE'
if 'C.TYPE' in all.obs.columns:
    # Replace NaN with 'NA'
    all.obs['C.TYPE'] = all.obs['C.TYPE'].replace(pd.NA, 'NA')
    
    # Replace 'LUMINALA' with 'Luminal A' and 'LUMINALB' with 'Luminal B'
    all.obs['C.TYPE'] = all.obs['C.TYPE'].replace({
        'LUMINALA': 'Luminal A',
        'LUMINALB': 'Luminal B',
        'TN': 'TNBC'
    })

    # Print the updated 'C.TYPE' values to confirm
    print(all.obs['C.TYPE'].value_counts())
else:
    print("'C.TYPE' column not found in all.obs.")


#
if 'C.TYPE' in all.obs.columns:
    # Filter to keep only cells where C.TYPE is 'Luminal A'
    all_luminal_a = all[all.obs['C.TYPE'] == 'Luminal A']
    
    # Optionally, check the number of observations remaining
    print(f"Number of observations after filtering for Luminal A: {all_luminal_a.obs.shape[0]}")
else:
    print("'C.TYPE' column not found in the obs DataFrame.")

# Subset to only IMMUNE and STROMA cells
adata_subset = all_luminal_a[all_luminal_a.obs['TYPE'].isin(['IMMUNE', 'STROMA'])].copy()

adata_subset.obs['cell_type'] = adata_subset.obs['TYPE']

# Run LIANA's rank aggregation method
li.method.rank_aggregate(
    adata=adata_subset,
    groupby='cell_type',      # Column in adata.obs to group cells by
    resource_name='consensus',# Resource database for ligand-receptor pairs
    expr_prop=0.1,            # Minimum expression proportion threshold
    min_cells=5,              # Minimum number of cells per group
    key_added='liana_res',    # Key to store results in adata.uns
    use_raw=False,            # Use raw attribute of adata if present
    n_perms=1000,             # Number of permutations for statistical testing
    seed=1337,                # Random seed for reproducibility
    verbose=True              # Verbosity flag
)


# Extract results
results_df = adata_subset.uns['liana_res']

# Save to Excel
output_path = "liana_results_STROMA_IMMUNE.xlsx"
results_df.to_excel(output_path, index=False)

print(f"LIANA results saved to {output_path}")


#Subset to IMMUNE and TUMOR cells#

adata_subset = all_luminal_a[all_luminal_a.obs['TYPE'].isin(['IMMUNE', 'TUMOR'])].copy()

# Define the column that contains the cell type label
adata_subset.obs['cell_type'] = adata_subset.obs['TYPE']

# Run LIANA's rank aggregation method
li.method.rank_aggregate(
    adata=adata_subset,
    groupby='cell_type',           # Column in adata.obs to group cells by
    resource_name='consensus',# Resource database for ligand-receptor pairs
    expr_prop=0.1,            # Minimum expression proportion threshold
    min_cells=5,              # Minimum number of cells per group
    key_added='liana_res',    # Key to store results in adata.uns
    use_raw=False,             # Use raw attribute of adata if present
    n_perms=1000,             # Number of permutations for statistical testing
    seed=1337,                # Random seed for reproducibility
    verbose=True              # Verbosity flag
)

# Save the results to Excel
results_df = adata_subset.uns['liana_res']
results_df.to_excel("liana_results_TUMOR_IMMUNE.xlsx", index=False)

print("LIANA results saved as 'liana_results_TUMOR_IMMUNE.xlsx'")


##Subset to STROMA and TUMOR cells##
adata_subset = all_luminal_a[all_luminal_a.obs['TYPE'].isin(['STROMA', 'TUMOR'])].copy()

#Assign 'cell_type' from TYPE
adata_subset.obs['cell_type'] = adata_subset.obs['TYPE']


# Run LIANA's rank aggregation method
li.method.rank_aggregate(
    adata=adata_subset,
    groupby='cell_type',           # Column in adata.obs to group cells by
    resource_name='consensus',# Resource database for ligand-receptor pairs
    expr_prop=0.1,            # Minimum expression proportion threshold
    min_cells=5,              # Minimum number of cells per group
    key_added='liana_res',    # Key to store results in adata.uns
    use_raw=False,             # Use raw attribute of adata if present
    n_perms=1000,             # Number of permutations for statistical testing
    seed=1337,                # Random seed for reproducibility
    verbose=True              # Verbosity flag
)
#Save results to Excel
results_df = adata_subset.uns['liana_res']
results_df.to_excel("liana_results_STROMA_TUMOR.xlsx", index=False)

print("LIANA results saved as 'liana_results_STROMA_TUMOR.xlsx'")


##SUMARY OF COMMON RESULTS##
# File paths
file_immune_stroma = r'path\liana_results_STROMA_IMMUNE.xlsx'
file_stroma_tumor = r'path\liana_results_STROMA_TUMOR.xlsx'
file_tumor_immune = r'path\liana_results_TUMOR_IMMUNE.xlsx'

# Load the files
df_immune_stroma = pd.read_excel(file_immune_stroma)
df_stroma_tumor = pd.read_excel(file_stroma_tumor)
df_tumor_immune = pd.read_excel(file_tumor_immune)

# Function to filter out self-interactions (same source and target cell type)
def filter_self_interactions(df):
    return df[df['source'] != df['target']]

# Apply filtering
df_immune_stroma = filter_self_interactions(df_immune_stroma)
df_stroma_tumor = filter_self_interactions(df_stroma_tumor)
df_tumor_immune = filter_self_interactions(df_tumor_immune)

# Create unique interaction ID that includes ligand → receptor direction
def create_interaction_id(df):
    return df['ligand_complex'] + '->' + df['receptor_complex']

df_immune_stroma['interaction'] = create_interaction_id(df_immune_stroma)
df_stroma_tumor['interaction'] = create_interaction_id(df_stroma_tumor)
df_tumor_immune['interaction'] = create_interaction_id(df_tumor_immune)

# Find common directional interactions across all 3
common_interactions = set(df_immune_stroma['interaction']) & \
                      set(df_stroma_tumor['interaction']) & \
                      set(df_tumor_immune['interaction'])

# Extract matching rows
common_df_immune_stroma = df_immune_stroma[df_immune_stroma['interaction'].isin(common_interactions)]
common_df_stroma_tumor = df_stroma_tumor[df_stroma_tumor['interaction'].isin(common_interactions)]
common_df_tumor_immune = df_tumor_immune[df_tumor_immune['interaction'].isin(common_interactions)]

# Add group tag and combine
common_df_immune_stroma['pair'] = 'IMMUNE-STROMA'
common_df_stroma_tumor['pair'] = 'STROMA-TUMOR'
common_df_tumor_immune['pair'] = 'TUMOR-IMMUNE'

# Combine all into one summary DataFrame
summary_df = pd.concat([
    common_df_immune_stroma,
    common_df_stroma_tumor,
    common_df_tumor_immune
], ignore_index=True)

# Save the result
output_path = r'path\liana_common_interactions.xlsx'
summary_df.to_excel(output_path, index=False)

print(f"Summary saved to: {output_path}")

# File paths
files = [
   r'path\liana_results_STROMA_IMMUNE.xlsx'
   r'path\liana_results_STROMA_TUMOR.xlsx'
   r'path\liana_results_TUMOR_IMMUNE.xlsx'

]

# Chemokine-related prefixes
chemokine_prefixes = ['CCL', 'CXCL', 'CX3CL', 'XCL']
receptor_prefixes = ['CCR', 'CXCR', 'CX3CR', 'XCR']

# List to collect filtered data
filtered_dfs = []

# Loop through each file
for file in files:
    df = pd.read_excel(file)

    # Drop NA and filter for p-value
    df = df[df['cellphone_pvals'].notna()]
    df = df[df['cellphone_pvals'] < 0.05]

    # Chemokine filtering
    is_chemokine_ligand = df['ligand_complex'].str.startswith(tuple(chemokine_prefixes))
    is_chemokine_receptor = df['receptor_complex'].str.startswith(tuple(receptor_prefixes))
    df = df[is_chemokine_ligand | is_chemokine_receptor]

    filtered_dfs.append(df)

# Combine all filtered interactions
combined_df = pd.concat(filtered_dfs, ignore_index=True)

# Remove self-interactions: STROMA-STROMA, IMMUNE-IMMUNE, TUMOR-TUMOR
self_interactions = (combined_df['source'] == combined_df['target']) & \
                    combined_df['source'].isin(['STROMA', 'IMMUNE', 'TUMOR'])
filtered_combined_df = combined_df[~self_interactions].copy()

# Save to Excel
filtered_combined_df.to_excel("Significant_Chemokine_Interactions_LA.xlsx", index=False)


###The same prossess was done on TNBC samples

all_TNBC = all[all.obs['C.TYPE'] == 'TNBC']
# Subset to only IMMUNE and STROMA cells
adata_subset = all_TNBC[all_TNBC.obs['TYPE'].isin(['IMMUNE', 'STROMA'])].copy()
adata_subset.obs['cell_type'] = adata_subset.obs['TYPE']

# Run LIANA's rank aggregation method
li.method.rank_aggregate(
    adata=adata_subset,
    groupby='cell_type',      # Column in adata.obs to group cells by
    resource_name='consensus',# Resource database for ligand-receptor pairs
    expr_prop=0.1,            # Minimum expression proportion threshold
    min_cells=5,              # Minimum number of cells per group
    key_added='liana_res',    # Key to store results in adata.uns
    use_raw=False,             # Use raw attribute of adata if present
    n_perms=1000,             # Number of permutations for statistical testing
    seed=1337,                # Random seed for reproducibility
    verbose=True              # Verbosity flag
)

# Extract results
results_df = adata_subset.uns['liana_res']
# Save to Excel
output_path = "TNBC_IMMUNE-STROMA_liana_results.xlsx"
results_df.to_excel(output_path, index=False)

print(f"LIANA results saved to {output_path}")

#Subset to IMMUNE and TUMOR cells
adata_subset = all_TNBC[all_TNBC.obs['TYPE'].isin(['IMMUNE', 'TUMOR'])].copy()

#Define the column that contains the cell type label
adata_subset.obs['cell_type'] = adata_subset.obs['TYPE']

# Run LIANA's rank aggregation method
li.method.rank_aggregate(
    adata=adata_subset,
    groupby='cell_type',      # Column in adata.obs to group cells by
    resource_name='consensus',# Resource database for ligand-receptor pairs
    expr_prop=0.1,            # Minimum expression proportion threshold
    min_cells=5,              # Minimum number of cells per group
    key_added='liana_res',    # Key to store results in adata.uns
    use_raw=False,             # Use raw attribute of adata if present
    n_perms=1000,             # Number of permutations for statistical testing
    seed=1337,                # Random seed for reproducibility
    verbose=True              # Verbosity flag
)

#Save the results to Excel
results_df = adata_subset.uns['liana_res']
results_df.to_excel("TNBC_liana_results_TUMOR_IMMUNE.xlsx", index=False)

print("LIANA results saved as 'liana_results_TUMOR_IMMUNE.xlsx'")

##Subset to STROMA and TUMOR cells##
adata_subset = all_TNBC[all_TNBC.obs['TYPE'].isin(['STROMA', 'TUMOR'])].copy()

#Assign 'cell_type' from TYPE
adata_subset.obs['cell_type'] = adata_subset.obs['TYPE']


# Run LIANA's rank aggregation method
li.method.rank_aggregate(
    adata=adata_subset,
    groupby='cell_type',      # Column in adata.obs to group cells by
    resource_name='consensus',# Resource database for ligand-receptor pairs
    expr_prop=0.1,            # Minimum expression proportion threshold
    min_cells=5,              # Minimum number of cells per group
    key_added='liana_res',    # Key to store results in adata.uns
    use_raw=False,            # Use raw attribute of adata if present
    n_perms=1000,             # Number of permutations for statistical testing
    seed=1337,                # Random seed for reproducibility
    verbose=True              # Verbosity flag
)
# Save results to Excel
results_df = adata_subset.uns['liana_res']
results_df.to_excel("TNBC_liana_results_STROMA_TUMOR.xlsx", index=False)

print("LIANA results saved as 'liana_results_STROMA_TUMOR.xlsx'")

# File paths
files = [
   r'path\TNBC_liana_results_STROMA_IMMUNE.xlsx'
   r'path\TNBC_liana_results_STROMA_TUMOR.xlsx'
   r'path\TNBC_liana_results_TUMOR_IMMUNE.xlsx'

]

# Chemokine-related prefixes
chemokine_prefixes = ['CCL', 'CXCL', 'CX3CL', 'XCL']
receptor_prefixes = ['CCR', 'CXCR', 'CX3CR', 'XCR']

# List to collect filtered data
filtered_dfs = []

# Loop through each file
for file in files:
    df = pd.read_excel(file)

    # Drop NA and filter for p-value
    df = df[df['cellphone_pvals'].notna()]
    df = df[df['cellphone_pvals'] < 0.05]

    # Chemokine filtering
    is_chemokine_ligand = df['ligand_complex'].str.startswith(tuple(chemokine_prefixes))
    is_chemokine_receptor = df['receptor_complex'].str.startswith(tuple(receptor_prefixes))
    df = df[is_chemokine_ligand | is_chemokine_receptor]

    filtered_dfs.append(df)

# Combine all filtered interactions
combined_df = pd.concat(filtered_dfs, ignore_index=True)

# Remove self-interactions: STROMA-STROMA, IMMUNE-IMMUNE, TUMOR-TUMOR
self_interactions = (combined_df['source'] == combined_df['target']) & \
                    combined_df['source'].isin(['STROMA', 'IMMUNE', 'TUMOR'])
filtered_combined_df = combined_df[~self_interactions].copy()

# Save to Excel
filtered_combined_df.to_excel("Significant_Chemokine_Interactions_TNBC.xlsx", index=False)


##The SignificantChemokine Interactions were manually annotated based on known litrature functions as:
##"Tumor Promoting, Tumor suppresing or Unknown/dual role.
##The two excel files were merged to one file with a column called C.TYPE contaning information about the ##Luminal A or TNBC BC subtype
##The annotated excel file was named:Annotated_CXC_Interactions

##The graph
# === Load and preprocess ===
def load_and_process(file_path, cancer_type):
    df = pd.read_excel(file_path)  # <-- FIXED: use read_excel for .xlsx files
    df = df[df['C.TYPE'] == cancer_type].copy()
    df['ligand_complex'] = df['ligand_complex'].str.upper()
    df['receptor_complex'] = df['receptor_complex'].str.upper()
    df['interaction_type'] = df['interaction_type'].fillna('Unknown')

    def get_cell_type(x):
        x = str(x).upper()
        if 'IMMUNE' in x:
            return 'IMMUNE'
        elif 'TUMOR' in x:
            return 'TUMOR'
        elif 'STROMA' in x:
            return 'STROMA'
        return 'OTHER'

    df['source_type'] = df['source'].apply(get_cell_type)
    df['target_type'] = df['target'].apply(get_cell_type)
    return df

def build_graph(df):
    G = nx.DiGraph()
    for _, row in df.iterrows():
        ligand = row['ligand_complex']
        receptor = row['receptor_complex']
        G.add_node(ligand, cell_type=row['source_type'])
        G.add_node(receptor, cell_type=row['target_type'])
        G.add_edge(ligand, receptor, interaction=row['interaction_type'])
    return G

# === Color maps ===
node_color_map = {
    'IMMUNE': '#aec7e8',
    'TUMOR': '#ff9896',
    'STROMA': '#98df8a',
    'OTHER': '#d3d3d3'
}
edge_color_map = {
    'Tumor Promoting': '#d62728',
    'Tumor suppresing': '#2a9df4',
    'Unknown': '#d3d3d3'
}
cell_type_priority = {'STROMA': 0, 'TUMOR': 1, 'IMMUNE': 2, 'OTHER': 3}

# === File path ===
file_path = r"path\Annotated_CXC_Interactions.xlsx"

# === Load and build graphs ===
df_luminal = load_and_process(file_path, "Luminal A")
df_tnbc = load_and_process(file_path, "TNBC")
G_luminal = build_graph(df_luminal)
G_tnbc = build_graph(df_tnbc)

# === Create figure with 2 plots ===
fig, axes = plt.subplots(1, 2, figsize=(8.27, 11.69 / 3.7), gridspec_kw={"wspace": 0.02})
titles = ["Luminal A L-R Interactions", "TNBC L-R Interactions"]
graphs = [G_luminal, G_tnbc]

for ax, G, title in zip(axes, graphs, titles):
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=9)

    # Extract and order nodes by cell type
    nodes_by_type = {k: [] for k in cell_type_priority.keys()}
    for node in G.nodes:
        cell_type = G.nodes[node].get("cell_type", "OTHER")
        nodes_by_type[cell_type].append(node)

    ordered_nodes = []
    for group in ['STROMA', 'TUMOR', 'IMMUNE', 'OTHER']:
        ordered_nodes.extend(sorted(nodes_by_type[group]))

    radius = 7
    theta = np.linspace(0, 2 * np.pi, len(ordered_nodes), endpoint=False)
    pos = {node: (radius * np.cos(t), radius * np.sin(t)) for node, t in zip(ordered_nodes, theta)}

    for node in G.nodes():
        x, y = pos[node]
        color = node_color_map.get(G.nodes[node]['cell_type'], '#d3d3d3')
        ax.scatter(x, y, s=500, color=color, zorder=2)

        angle_rad = np.arctan2(y, x)
        angle_deg = np.degrees(angle_rad)
        rotation = angle_deg if -90 <= angle_deg <= 90 else angle_deg + 180

        ax.text(
            x, y, node, fontsize=6.5,
            ha='center', va='center',
            rotation=rotation,
            rotation_mode='anchor',
            zorder=3
        )

    for u, v in G.edges():
        if u not in pos or v not in pos:
            continue
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        verts = [(x1, y1), (0, 0), (x2, y2)]
        codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
        path = Path(verts, codes)
        edge_type = G[u][v]['interaction']
        edge_color = edge_color_map.get(edge_type, '#d3d3d3')
        patch = PathPatch(path, edgecolor=edge_color, lw=1.2, facecolor='none', zorder=1)
        ax.add_patch(patch)

# === Add legend to the right side ===
legend_elements = [
    Patch(facecolor=node_color_map[k], edgecolor='black', label=k)
    for k in ['IMMUNE', 'TUMOR', 'STROMA']
] + [
    Patch(facecolor='none', edgecolor=edge_color_map['Pro Tumor'], lw=1.5, label='Pro Tumor'),
    Patch(facecolor='none', edgecolor=edge_color_map['Anti Tumor'], lw=1.5, label='Anti Tumor'),
    Patch(facecolor='none', edgecolor=edge_color_map['Unknown'], lw=1.5, label='Unknown / Dual')
]

fig.legend(
    handles=legend_elements,
    loc='center left',
    bbox_to_anchor=(1.01, 0.5),
    ncol=1,
    fontsize=7.5,
    frameon=False
)

# === Save ===
output_path = r"path\combined_luminalA_TNBC_chord_LEGEND_RIGHT.pdf"
plt.savefig(output_path, format='pdf', bbox_inches='tight')
plt.show()

print(f"Plot saved with legend on the right:\n{output_path}")

####Luminal A Stage specific interactions
#the same prosess was used for ligand receptor analysis in local vs.Metastatic samples:
#3 files were generated for local and for Metastatic samples:

#Local Breast Luminal A samples:
# File paths
files = [
    "liana_results_STROMA_TUMOR_LA_BREAST.xlsx",
    "liana_results_TUMOR_IMMUNE-LA_BREAST.xlsx",
    "IMMUNE-STROMA_liana_results_LA_BREAST.xlsx"
]

# Chemokine-related prefixes
chemokine_prefixes = ['CCL', 'CXCL', 'CX3CL', 'XCL']
receptor_prefixes = ['CCR', 'CXCR', 'CX3CR', 'XCR']

# List to collect filtered data
filtered_dfs = []

# Loop through each file
for file in files:
    df = pd.read_excel(file)

    # Drop NA and filter for p-value
    df = df[df['cellphone_pvals'].notna()]
    df = df[df['cellphone_pvals'] < 0.05]

    # Chemokine filtering
    is_chemokine_ligand = df['ligand_complex'].str.startswith(tuple(chemokine_prefixes))
    is_chemokine_receptor = df['receptor_complex'].str.startswith(tuple(receptor_prefixes))
    df = df[is_chemokine_ligand | is_chemokine_receptor]

    filtered_dfs.append(df)

# Combine all filtered interactions
combined_df = pd.concat(filtered_dfs, ignore_index=True)

# Remove self-interactions: STROMA-STROMA, IMMUNE-IMMUNE, TUMOR-TUMOR
self_interactions = (combined_df['source'] == combined_df['target']) & \
                    combined_df['source'].isin(['STROMA', 'IMMUNE', 'TUMOR'])
filtered_combined_df = combined_df[~self_interactions].copy()

# Save to Excel
filtered_combined_df.to_excel("Significant_Chemokine_Interactions_LA_BREAST.xlsx", index=False)

main_file = "Significant_Chemokine_Interactions_LA_BREAST.xlsx"
main_df = pd.read_excel(main_file)

# === Load the annotation file with interaction_type ===
annotation_file = r'C:\Users\kerenylab\Documents\breast cancer\filtered_CXC_interactions.xlsx'
annotation_df = pd.read_excel(annotation_file)

# === Merge based on all relevant columns ===
merged_df = pd.merge(
    main_df,
    annotation_df[['source', 'target', 'ligand_complex', 'receptor_complex', 'interaction_type']],
    how='left',
    on=['source', 'target', 'ligand_complex', 'receptor_complex']
)

# === Fill unmatched interaction_type as NA ===
merged_df['interaction_type'] = merged_df['interaction_type'].fillna(pd.NA)

# === Save the final annotated file ===
output_file = "Significant_Chemokine_Interactions_LA_BREAST_annotated.xlsx"
merged_df.to_excel(output_file, index=False)

#Metastatic Luminal A samples
# File paths
files = [
    "liana_results_STROMA_TUMOR_LA_MET.xlsx",
    "liana_results_TUMOR_IMMUNE-LA_MET.xlsx",
    "IMMUNE-STROMA_liana_results_LA_MET.xlsx"
]

# Chemokine-related prefixes
chemokine_prefixes = ['CCL', 'CXCL', 'CX3CL', 'XCL']
receptor_prefixes = ['CCR', 'CXCR', 'CX3CR', 'XCR']

# List to collect filtered data
filtered_dfs = []

# Loop through each file
for file in files:
    df = pd.read_excel(file)

    # Drop NA and filter for p-value
    df = df[df['cellphone_pvals'].notna()]
    df = df[df['cellphone_pvals'] < 0.05]

    # Chemokine filtering
    is_chemokine_ligand = df['ligand_complex'].str.startswith(tuple(chemokine_prefixes))
    is_chemokine_receptor = df['receptor_complex'].str.startswith(tuple(receptor_prefixes))
    df = df[is_chemokine_ligand | is_chemokine_receptor]

    filtered_dfs.append(df)

# Combine all filtered interactions
combined_df = pd.concat(filtered_dfs, ignore_index=True)

# Remove self-interactions: STROMA-STROMA, IMMUNE-IMMUNE, TUMOR-TUMOR
self_interactions = (combined_df['source'] == combined_df['target']) & \
                    combined_df['source'].isin(['STROMA', 'IMMUNE', 'TUMOR'])
filtered_combined_df = combined_df[~self_interactions].copy()

# Save to Excel
filtered_combined_df.to_excel("Significant_Chemokine_Interactions_LA_MET.xlsx", index=False)

main_file = "Significant_Chemokine_Interactions_LA_MET.xlsx"
main_df = pd.read_excel(main_file)

# === Load the annotation file with interaction_type ===
annotation_file = r'C:\Users\kerenylab\Documents\breast cancer\Annotated_CXC_Interactions.csv'
annotation_df = pd.read_csv(annotation_file)

# === Merge based on all relevant columns ===
merged_df = pd.merge(
    main_df,
    annotation_df[['source', 'target', 'ligand_complex', 'receptor_complex', 'interaction_type']],
    how='left',
    on=['source', 'target', 'ligand_complex', 'receptor_complex']
)

# === Fill unmatched interaction_type as NA ===
merged_df['interaction_type'] = merged_df['interaction_type'].fillna(pd.NA)

# === Save the final annotated file ===
output_file = "Significant_Chemokine_Interactions_LA_MET_annotated.xlsx"
merged_df.to_excel(output_file, index=False)


##Graphs:
#Metastatic samples interactions-graph

# === Load data ===
file_path = r"path\Significant_Chemokine_Interactions_LA_MET_annotated.xlsx"
df = pd.read_excel(file_path)
#df = df[df['C.TYPE'] == 'TNBC'].copy()

# Clean and classify
df['ligand_complex'] = df['ligand_complex'].str.upper()
df['receptor_complex'] = df['receptor_complex'].str.upper()
df['interaction_type'] = df['interaction_type'].fillna('Unknown')

def get_cell_type(x):
    x = str(x).upper()
    if 'IMMUNE' in x:
        return 'IMMUNE'
    elif 'TUMOR' in x:
        return 'TUMOR'
    elif 'STROMA' in x:
        return 'STROMA'
    return 'OTHER'

df['source_type'] = df['source'].apply(get_cell_type)
df['target_type'] = df['target'].apply(get_cell_type)

# === Create graph ===
G = nx.DiGraph()
for _, row in df.iterrows():
    ligand = row['ligand_complex']
    receptor = row['receptor_complex']
    G.add_node(ligand, cell_type=row['source_type'])
    G.add_node(receptor, cell_type=row['target_type'])
    G.add_edge(ligand, receptor, interaction=row['interaction_type'])

# === Color maps ===
node_color_map = {
    'IMMUNE': '#aec7e8',
    'TUMOR': '#ff9896',
    'STROMA': '#98df8a',
    'OTHER': '#d3d3d3'
}
edge_color_map = {
    'Tumor Promoting': '#d62728',
    'Tumor Suppresing': '#2a9df4',
    'Unknown': '#d3d3d3'
}

# === Group and sort nodes ===
grouped_nodes = []
for group in ['IMMUNE', 'TUMOR', 'STROMA', 'OTHER']:
    group_nodes = sorted([n for n in G.nodes if G.nodes[n]['cell_type'] == group])
    grouped_nodes.extend(group_nodes)

# === Circular layout ===
radius = 2
theta = np.linspace(0, 2 * np.pi, len(grouped_nodes), endpoint=False)
pos = {node: (radius * np.cos(t), radius * np.sin(t)) for node, t in zip(grouped_nodes, theta)}

# === Plot ===
fig, ax = plt.subplots(figsize=(8.27, 11.69 / 4.2))
ax.set_aspect('equal')
ax.axis('off')

# Draw nodes and labels with rotation
for node in G.nodes():
    x, y = pos[node]
    color = node_color_map.get(G.nodes[node]['cell_type'], '#d3d3d3')
    ax.scatter(x, y, s=500, color=color, zorder=2)
    
    # Compute angle for rotation
    angle_rad = theta[grouped_nodes.index(node)]
    angle_deg = np.degrees(angle_rad)
    rotation = angle_deg if angle_deg <= 180 else angle_deg - 180
    
    ax.text(
        x, y, node, fontsize=7,
        ha='center', va='center',
        rotation=rotation,
        rotation_mode='anchor',
        zorder=3
    )

# Draw curved edges
for u, v in G.edges():
    x1, y1 = pos[u]
    x2, y2 = pos[v]
    verts = [(x1, y1), (0, 0), (x2, y2)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    path = Path(verts, codes)
    edge_color = edge_color_map.get(G[u][v]['interaction'], '#d3d3d3')
    patch = PathPatch(path, edgecolor=edge_color, lw=1.2, facecolor='none', zorder=1)
    ax.add_patch(patch)

# === Title and save ===
plt.title("LA Metastatic Samples\n L-R Interactions", fontsize=9)
plt.tight_layout()
output_path = r"path\MET_grouped_tight_chord_FINAL.pdf"
plt.savefig(output_path, format='pdf', bbox_inches='tight')
plt.show()

print(f"✅ Plot saved as PDF to:\n{output_path}")

#Primary breast samples interactions-graph
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import PathPatch
from matplotlib.path import Path

# === Load data ===
file_path = r"path\Significant_Chemokine_Interactions_LA_BREAST_annotated.xlsx"
df = pd.read_excel(file_path)
#df = df[df['C.TYPE'] == 'TNBC'].copy()

# Clean and classify
df['ligand_complex'] = df['ligand_complex'].str.upper()
df['receptor_complex'] = df['receptor_complex'].str.upper()
df['interaction_type'] = df['interaction_type'].fillna('Unknown')

def get_cell_type(x):
    x = str(x).upper()
    if 'IMMUNE' in x:
        return 'IMMUNE'
    elif 'TUMOR' in x:
        return 'TUMOR'
    elif 'STROMA' in x:
        return 'STROMA'
    return 'OTHER'

df['source_type'] = df['source'].apply(get_cell_type)
df['target_type'] = df['target'].apply(get_cell_type)

# === Create graph ===
G = nx.DiGraph()
for _, row in df.iterrows():
    ligand = row['ligand_complex']
    receptor = row['receptor_complex']
    G.add_node(ligand, cell_type=row['source_type'])
    G.add_node(receptor, cell_type=row['target_type'])
    G.add_edge(ligand, receptor, interaction=row['interaction_type'])

# === Color maps ===
node_color_map = {
    'IMMUNE': '#aec7e8',
    'TUMOR': '#ff9896',
    'STROMA': '#98df8a',
    'OTHER': '#d3d3d3'
}
edge_color_map = {
    'Pro Tumor': '#d62728',
    'Anti Tumor': '#2a9df4',
    'Unknown': '#d3d3d3'
}

# === Group and sort nodes ===
grouped_nodes = []
for group in ['IMMUNE', 'TUMOR', 'STROMA', 'OTHER']:
    group_nodes = sorted([n for n in G.nodes if G.nodes[n]['cell_type'] == group])
    grouped_nodes.extend(group_nodes)

# === Circular layout ===
radius = 2
theta = np.linspace(0, 2 * np.pi, len(grouped_nodes), endpoint=False)
pos = {node: (radius * np.cos(t), radius * np.sin(t)) for node, t in zip(grouped_nodes, theta)}

# === Plot ===
fig, ax = plt.subplots(figsize=(8.27, 11.69 / 4.2))
ax.set_aspect('equal')
ax.axis('off')

# Draw nodes and labels with rotation
for node in G.nodes():
    x, y = pos[node]
    color = node_color_map.get(G.nodes[node]['cell_type'], '#d3d3d3')
    ax.scatter(x, y, s=500, color=color, zorder=2)
    
    # Compute angle for rotation
    angle_rad = theta[grouped_nodes.index(node)]
    angle_deg = np.degrees(angle_rad)
    rotation = angle_deg if angle_deg <= 180 else angle_deg - 180
    
    ax.text(
        x, y, node, fontsize=7,
        ha='center', va='center',
        rotation=rotation,
        rotation_mode='anchor',
        zorder=3
    )

# Draw curved edges
for u, v in G.edges():
    x1, y1 = pos[u]
    x2, y2 = pos[v]
    verts = [(x1, y1), (0, 0), (x2, y2)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    path = Path(verts, codes)
    edge_color = edge_color_map.get(G[u][v]['interaction'], '#d3d3d3')
    patch = PathPatch(path, edgecolor=edge_color, lw=1.2, facecolor='none', zorder=1)
    ax.add_patch(patch)

# === Title and save ===
plt.title("LA Primary Samples\n L-R Interactions", fontsize=9)
plt.tight_layout()
output_path = r"path\BREAST_grouped_tight_chord_FINAL.pdf"
plt.savefig(output_path, format='pdf', bbox_inches='tight')
plt.show()

print(f"✅ Plot saved as PDF to:\n{output_path}")
