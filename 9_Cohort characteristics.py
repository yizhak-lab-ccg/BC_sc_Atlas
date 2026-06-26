import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import pandas as pd
import numpy as np
from matplotlib import rcParams

#load integrated dataset
all=sc.read_h5ad(r'integrated dataset.h5ad')
all


#Pie chart-Geographic distribution of cohort samples by country
#C.TYPE=BC Subtype
#P_P=sample
#BIOPSY=Biopsy location
# Custom pink shades palette
pink_shades = [
    '#fde0dd', '#fcc5c0', '#fa9fb5', '#f768a1',
    '#dd3497', '#ae017e', '#7a0177', '#49006a'
]

# Filter out 'nan'
filtered = all.obs[~all.obs['ORIGIN'].isin(['nan'])].copy()

# Remove unused categories if ORIGIN is categorical
if filtered['ORIGIN'].dtype.name == 'category':
    filtered['ORIGIN'] = filtered['ORIGIN'].cat.remove_unused_categories()

# Capitalize first letter, keep 'USA' unchanged
filtered['ORIGIN'] = filtered['ORIGIN'].apply(
    lambda x: 'USA' if str(x).upper() == 'USA' else str(x).capitalize()
)

# Calculate percentages
origin_counts = filtered['ORIGIN'].value_counts(normalize=True) * 100

# Plot pie chart
fig, ax = plt.subplots(figsize=(4.135, 2.783))  # 8.27 / 2 by 11.69 / 4.2

# Fonts
plt.rcParams.update({
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8
})

# Pie chart
wedges, texts, autotexts = ax.pie(
    origin_counts, 
    labels=origin_counts.index, 
    autopct='%1.1f%%',
    startangle=140,
    colors=pink_shades[:len(origin_counts)],
    textprops={'fontsize': 8}
)

# Title
ax.set_title('Samples Origin by Country', fontsize=9, pad=15)
ax.axis('equal')

# Save and show
plt.tight_layout()
plt.savefig("origin_piechart.pdf", format='pdf', bbox_inches='tight')
plt.show()

#BC subtype pie chart disstribution
#Font settings
rcParams.update({
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8
})

# Create BIOPSY_GROUP column and filter
all.obs['BIOPSY_GROUP'] = all.obs['BIOPSY'].apply(
    lambda x: 'BREAST' if x == 'BREAST' else ('LN' if x == 'LN' else ('MET' if x not in ['BLOOD', 'NA'] else None))
)
all = all[all.obs['BIOPSY_GROUP'].notnull()]
adata_filtered = all[all.obs['C.TYPE'] != 'UNKNOWN'].copy()

# Calculate C.TYPE percentages
ctype_counts = adata_filtered.obs['C.TYPE'].value_counts(normalize=True) * 100
ctype_labels = ctype_counts.index.tolist()
ctype_values = ctype_counts.values

# Define color palette
main_colors = ['#fee8c8', '#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f', '#b30000', '#7f0000']

# Create figure
fig, ax = plt.subplots(figsize=(4.135, 2.783))

# Plot pie chart with external labels and fixed radius
wedges, texts, autotexts = ax.pie(
    ctype_values,
    labels=ctype_labels,
    autopct='%1.1f%%',
    startangle=140,
    colors=main_colors[:len(ctype_values)],
    textprops={'fontsize': 8},
    radius=1.0,
    labeldistance=1.2
)

# Title and layout
ax.set_title('BC Subtype Distribution', fontsize=9, pad=15)
ax.axis('equal')

# Manually expand pie plot to use more space in axis
ax.set_position([0.05, 0.05, 0.9, 0.9])  # [left, bottom, width, height]

plt.savefig("C_TYPE_piechart_matched_with_labels.pdf", format='pdf', bbox_inches='tight')
plt.show()


#Bar Plot-Breakdown of sample collection sites and corresponding cell types by breast cancer subtype.
# Filter only major subtypes
major_ctypes = ['TNBC', 'Luminal A', 'Luminal B', 'HER2']
adata_bar = adata_filtered[adata_filtered.obs['C.TYPE'].isin(major_ctypes)].copy()

# Calculate percentages of TYPE per C.TYPE
type_props = (
    adata_bar.obs.groupby(['C.TYPE', 'TYPE'])
    .size()
    .groupby(level=0)
    .apply(lambda x: x / x.sum() * 100)
    .unstack(fill_value=0)
    .loc[major_ctypes]
)

# Calculate percentages of BIOPSY_GROUP per C.TYPE
biopsy_props = (
    adata_bar.obs.groupby(['C.TYPE', 'BIOPSY_GROUP'])
    .size()
    .groupby(level=0)
    .apply(lambda x: x / x.sum() * 100)
    .unstack(fill_value=0)
    .loc[major_ctypes]
)

# Set figure size
fig, ax = plt.subplots(figsize=(5, 2.783))

# Font size settings
plt.rcParams.update({
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8
})

width = 0.35
x = np.arange(len(major_ctypes))

# === Left bars (TYPE) ===
bottom = np.zeros(len(major_ctypes))
for group, color in sub_colors.items():
    values = type_props[group] if group in type_props else np.zeros(len(major_ctypes))
    bars = ax.bar(x - width/2, values, width, bottom=bottom, color=color)
    for i, val in enumerate(values):
        if val > 1:
            ax.text(x[i] - width/2, bottom[i] + val / 2, f"{val:.0f}%", ha='center', va='center', fontsize=8)
    bottom += values

# === Right bars (BIOPSY_GROUP) ===
bottom = np.zeros(len(major_ctypes))
for group, color in biopsy_colors.items():
    values = biopsy_props[group] if group in biopsy_props else np.zeros(len(major_ctypes))
    bars = ax.bar(x + width/2, values, width, bottom=bottom, color=color)
    for i, val in enumerate(values):
        if val > 1:
            ax.text(x[i] + width/2, bottom[i] + val / 2, f"{val:.0f}%", ha='center', va='center', fontsize=8)
    bottom += values

# Axis settings
ax.set_ylabel('Percentage', fontsize=8)
ax.set_title('Cell Type and Biopsy Group Distribution by BC Subtype', fontsize=9)
ax.set_xticks(x)
ax.set_xticklabels(major_ctypes, fontsize=8)
ax.tick_params(axis='y', labelsize=8)

# === Custom grouped legend on the right with renamed labels ===
legend_elements = []

# Group 1: TYPE (Cell Annotation)
legend_elements.append(Patch(color='white', label='Cell Annotation'))
type_labels = {'IMMUNE': 'Immune', 'STROMA': 'Stroma', 'TUMOR': 'Tumor'}
for key, color in sub_colors.items():
    legend_elements.append(Patch(facecolor=color, label=f'  {type_labels.get(key, key)}'))

# Group 2: BIOPSY_GROUP (Biopsy Site)
legend_elements.append(Patch(color='white', label='Biopsy Site'))
biopsy_labels = {'BREAST': 'Breast', 'LN': 'L/N', 'MET': 'Metastasis'}
for key, color in biopsy_colors.items():
    legend_elements.append(Patch(facecolor=color, label=f'  {biopsy_labels.get(key, key)}'))

# Add legend to the right
ax.legend(
    handles=legend_elements,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    fontsize=8,
    title=None
)

plt.tight_layout()
plt.savefig("C_TYPE_barplot_with_renamed_legend.pdf", format='pdf', bbox_inches='tight')
plt.show()

