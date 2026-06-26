import scanpy as sc
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import numpy as np
import pandas as pd

#load integrated tumor file
adata=sc.read_h5ad('path\tumor.h5ad')
adata

###Fraction of ERBB2 positive cells
#C.TYPE=BC Subtype
#P_P=sample
#BIOPSY=Biopsy location

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import numpy as np
import itertools

# Define valid C.TYPE categories
valid_types = ['TNBC', 'HER2', 'Luminal A', 'Luminal B']

# Subset the AnnData
adata_filtered = adata[adata.obs['C.TYPE'].isin(valid_types)].copy()

# Get ERBB2 expression
erbb2_expr = adata_filtered[:, 'ERBB2'].X
if not isinstance(erbb2_expr, np.ndarray):
    erbb2_expr = erbb2_expr.toarray().flatten()
else:
    erbb2_expr = erbb2_expr.flatten()

# Define ERBB2-positive
adata_filtered.obs['ERBB2_positive'] = erbb2_expr > 0.5

# Create DataFrame
df = adata_filtered.obs[['C.TYPE', 'P_P', 'ERBB2_positive']].copy()

# Calculate fraction of ERBB2-positive cells per P_P and C.TYPE
fraction_df = (
    df.groupby(['P_P', 'C.TYPE'])
    .agg(Fraction=('ERBB2_positive', 'mean'))
    .reset_index()
)

# Order C.TYPE by mean fraction
ordered_types = (
    fraction_df.groupby('C.TYPE')['Fraction']
    .mean()
    .sort_values()
    .index
)
fraction_df['C.TYPE'] = pd.Categorical(fraction_df['C.TYPE'], categories=ordered_types, ordered=True)

# Plot
plt.figure(figsize=(10, 6))
sns.boxplot(data=fraction_df, x='C.TYPE', y='Fraction', order=ordered_types, palette='Set2')
sns.stripplot(data=fraction_df, x='C.TYPE', y='Fraction', color='black', jitter=True, alpha=0.6, order=ordered_types)

# Calculate p-values for pairwise comparisons
p_values = {}
pairs = list(itertools.combinations(ordered_types, 2))
max_y = fraction_df['Fraction'].max()
offset = 0.05

for i, (g1, g2) in enumerate(pairs):
    vals1 = fraction_df[fraction_df['C.TYPE'] == g1]['Fraction'].dropna()
    vals2 = fraction_df[fraction_df['C.TYPE'] == g2]['Fraction'].dropna()
    
    if len(vals1) > 0 and len(vals2) > 0:
        stat, p = mannwhitneyu(vals1, vals2, alternative='two-sided')
        p_values[(g1, g2)] = p

        # Annotate plot
        x1, x2 = ordered_types.tolist().index(g1), ordered_types.tolist().index(g2)
        y = max_y + (i + 1) * offset
        plt.plot([x1, x1, x2, x2], [y, y + 0.01, y + 0.01, y], lw=1.5, c='black')
        plt.text((x1 + x2) / 2, y + 0.015, f"p={p:.3g}", ha='center', va='bottom', fontsize=8)

# Final plot settings
plt.title("Fraction of ERBB2-Positive Cells by Subtype", fontsize=14)
plt.ylabel("Fraction of ERBB2+ Cells", fontsize=12)
plt.xlabel("BC Subtype", fontsize=12)
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.savefig("ERBB2_Fraction_Boxplot_with_pvalues.pdf", format='pdf', bbox_inches='tight')
plt.show()


###ER positive cells###
# Define pastel scientific color palette
custom_palette = {
    'BREAST': '#A6CEE3',  # Light Blue
    'LN': '#B2DF8A',      # Light Green
    'MET': '#FB9A99'      # Light Red/Pink
}

# Superscript translation map
superscript_map = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

# Ordered groups
ordered_groups = ['BREAST', 'LN', 'MET']

# Save function
def save_plot(fig_size, filename, show_plot=False):
    fig, ax = plt.subplots(figsize=fig_size)

    # Boxplot
    sns.boxplot(
        data=counts_df,
        x='BIOPSY_GROUP',
        y='ESR1_positive_percentage',
        order=ordered_groups,
        showfliers=False,
        palette=custom_palette,
        ax=ax
    )

    # Stripplot
    sns.stripplot(
        data=counts_df,
        x='BIOPSY_GROUP',
        y='ESR1_positive_percentage',
        color='black',
        size=3,
        jitter=True,
        alpha=0.6,
        order=ordered_groups,
        ax=ax
    )

    # Calculate p-values
    pairs = [('BREAST', 'LN'), ('BREAST', 'MET'), ('LN', 'MET')]
    p_values = {}
    for pair in pairs:
        group1 = counts_df[counts_df['BIOPSY_GROUP'] == pair[0]]['ESR1_positive_percentage']
        group2 = counts_df[counts_df['BIOPSY_GROUP'] == pair[1]]['ESR1_positive_percentage']
        stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
        p_values[pair] = p

    # Determine Y-limits before adding annotations
    data_max = counts_df['ESR1_positive_percentage'].max()
    annotation_space = 12 * len(p_values)
    ymax = data_max + annotation_space
    ax.set_ylim(0, ymax)

    # Add annotation lines + p-values
    line_height = 2
    text_offset = 2.5
    base_y = data_max + 5

    for i, (pair, p) in enumerate(p_values.items()):
        if p < 0.05:
            x1, x2 = ordered_groups.index(pair[0]), ordered_groups.index(pair[1])
            y_line = base_y + i * 12
            y_text = y_line + line_height + text_offset

            ax.plot([x1, x1, x2, x2], [y_line, y_line + line_height, y_line + line_height, y_line],
                    lw=1.2, c='black')

            if p < 1e-16:
                p_label = "p < 1×10" + "⁻¹⁶"
            elif p >= 0.0001:
                p_label = f"p = {round(p, 3)}"
            else:
                exponent = int(np.floor(np.log10(p)))
                coeff = p / (10**exponent)
                exponent_str = str(exponent).translate(superscript_map)
                p_label = f"p = {coeff:.1f}×10{exponent_str}"

            ax.text((x1 + x2) / 2, y_text, p_label, ha='center', va='bottom', fontsize=8)

    # Labels
    ax.set_title("Percentage of ER+ Cells by Stage", fontsize=9)
    ax.set_xlabel("Biopsy Site", fontsize=8)
    ax.set_ylabel("Percentage of ER+ Cells", fontsize=8)
    ax.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig(filename, format='pdf', bbox_inches='tight')

    if show_plot:
        plt.show()
    else:
        plt.close()

# Output folder
save_dir = r"path for LUMINAL A figures"
os.makedirs(save_dir, exist_ok=True)

# File paths
small_path = os.path.join(save_dir, "ESR1_Positive_Percentage_FIXED_small.pdf")
full_path = os.path.join(save_dir, "ESR1_Positive_Percentage_FIXED_full.pdf")

# Save both figures
save_plot((8.27 / 3.5, 11.69 / 4.2), small_path)
save_plot((8.27, 11.69), full_path, show_plot=True)


###PR positive cells###

#Define pastel scientific color palette
custom_palette = {
    'BREAST': '#A6CEE3',  # Light Blue
    'LN': '#B2DF8A',      # Light Green
    'MET': '#FB9A99'      # Light Red/Pink
}

# Superscript translation map
superscript_map = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

# Ordered groups
ordered_groups = ['BREAST', 'LN', 'MET']

# === Extract PGR expression from adata_luminal_a ===
adata_luminal_a.obs['PGR_expr'] = adata_luminal_a[:, 'PGR'].X.toarray().flatten() \
    if not isinstance(adata_luminal_a[:, 'PGR'].X, np.ndarray) else adata_luminal_a[:, 'PGR'].X.flatten()

# Prepare dataframe
counts_df = adata_luminal_a.obs[['BIOPSY_GROUP', 'PGR_expr']].copy()
counts_df = counts_df[counts_df['BIOPSY_GROUP'].isin(ordered_groups)].copy()

# === Plotting Function ===
def save_plot(fig_size, filename, show_plot=False):
    fig, ax = plt.subplots(figsize=fig_size)

    # Boxplot
    sns.boxplot(
        data=counts_df,
        x='BIOPSY_GROUP',
        y='PGR_expr',
        order=ordered_groups,
        showfliers=False,
        palette=custom_palette,
        ax=ax
    )

    # Stripplot
    sns.stripplot(
        data=counts_df,
        x='BIOPSY_GROUP',
        y='PGR_expr',
        color='black',
        size=3,
        jitter=True,
        alpha=0.6,
        order=ordered_groups,
        ax=ax
    )

    # Calculate p-values
    pairs = [('BREAST', 'LN'), ('BREAST', 'MET'), ('LN', 'MET')]
    p_values = {}
    for pair in pairs:
        group1 = counts_df[counts_df['BIOPSY_GROUP'] == pair[0]]['PGR_expr']
        group2 = counts_df[counts_df['BIOPSY_GROUP'] == pair[1]]['PGR_expr']
        stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
        p_values[pair] = p

    # === Set fixed Y-limit and annotation range ===
    ax.set_ylim(0, 10)
    line_height = 0.3
    text_offset = 0.4
    base_y = 7.0

    # Add annotation lines + p-values
    for i, (pair, p) in enumerate(p_values.items()):
        if p < 0.05:
            x1, x2 = ordered_groups.index(pair[0]), ordered_groups.index(pair[1])
            y_line = base_y + i * 1.0
            y_text = y_line + line_height + text_offset

            ax.plot([x1, x1, x2, x2], [y_line, y_line + line_height, y_line + line_height, y_line],
                    lw=1.2, c='black')

            if p < 1e-16:
                p_label = "p < 1×10" + "⁻¹⁶"
            elif p >= 0.0001:
                p_label = f"p = {round(p, 3)}"
            else:
                exponent = int(np.floor(np.log10(p)))
                coeff = p / (10**exponent)
                exponent_str = str(exponent).translate(superscript_map)
                p_label = f"p = {coeff:.1f}×10{exponent_str}"

            ax.text((x1 + x2) / 2, y_text, p_label, ha='center', va='bottom', fontsize=8)

    # Labels
    ax.set_title("PGR Expression by Biopsy Site", fontsize=9)
    ax.set_xlabel("Biopsy Site", fontsize=8)
    ax.set_ylabel("PGR Expression", fontsize=8)
    ax.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig(filename, format='pdf', bbox_inches='tight')

    if show_plot:
        plt.show()
    else:
        plt.close()

# === Save plots ===
save_dir = r"path"
os.makedirs(save_dir, exist_ok=True)

# File paths
small_path = os.path.join(save_dir, "PGR_Expression_FIXED_small.pdf")
full_path = os.path.join(save_dir, "PGR_Expression_FIXED_full.pdf")

# Save and display
save_plot((8.27 / 3.5, 11.69 / 4.2), small_path, show_plot=True)
save_plot((8.27, 11.69), full_path, show_plot=False)


###ER positive vs negative DEG###
# Define the gene prefixes to remove
genes_to_remove_prefixes = ('RP', 'LIN', 'AC', 'MT-')

# Identify genes to remove based on their prefixes
genes_to_remove = [gene for gene in adata_luminal_a.var_names if gene.startswith(genes_to_remove_prefixes)]

# Subset the AnnData object to exclude the identified genes
adata_luminal_a_filtered = adata_luminal_a[:, ~adata_luminal_a.var_names.isin(genes_to_remove)].copy()

# Output the result
print(f"Filtered AnnData object shape: {adata_luminal_a_filtered.shape}")
adata_luminal_a=adata_luminal_a_filtered

# Check if 'ESR1' exists in adata.var_names
if 'ESR1' not in adata.var_names:
    raise ValueError("The gene 'ESR1' is not found in adata.var_names.")

# Extract ESR1 expression from adata.X
if isinstance(adata_luminal_a.X, np.ndarray):
    esr1_expression = adata_luminal_a[:, adata_luminal_a.var_names == 'ESR1'].X.flatten()
else:
    esr1_expression = adata_luminal_a[:, adata_luminal_a.var_names == 'ESR1'].X.toarray().flatten()

# Add ESR1 expression as a column in adata.obs
adata_luminal_a.obs['ESR1_expression'] = esr1_expression

# Define cells expressing ESR1 (expression > 0.5) vs non-expressing
adata_luminal_a.obs['ESR1_positive'] = adata_luminal_a.obs['ESR1_expression'].apply(
    lambda x: 'POSITIVE' if x > 0.5 else 'NEGATIVE'
).astype('category')

# Perform rank_genes_groups for DEG analysis
sc.tl.rank_genes_groups(
    adata_luminal_a,
    groupby='ESR1_positive',
    use_raw=False,  # Explicitly state not to use raw data
    method='wilcoxon'
)


###EXCEL FILE FOR DEG###
# Create a Pandas Excel writer to save to multiple sheets
with pd.ExcelWriter('DEG_ESR1_POSITIVE_vs_NEGATIVE.xlsx', engine='xlsxwriter') as writer:
    for group in adata_luminal_a_filtered.uns['rank_genes_groups']['names'].dtype.names:
        # Create a dataframe for each group containing gene names, p-values, log fold change, and adjusted p-values
        genes_df = pd.DataFrame({
            'Gene': adata_luminal_a_filtered.uns['rank_genes_groups']['names'][group],
            'Log2 Fold Change': adata_luminal_a_filtered.uns['rank_genes_groups']['logfoldchanges'][group],
            'p-value': adata_luminal_a_filtered.uns['rank_genes_groups']['pvals'][group],
            'Adjusted p-value': adata_luminal_a_filtered.uns['rank_genes_groups']['pvals_adj'][group]
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
print("Filtered DEG results saved to 'DEG_ESR1_POSITIVE_vs_NEGATIVE.xlsx'.")


###HEATMAP OF THE TOP 20 DEG###
# === Extract top 20 genes per group ===
n_top_genes = 20
result = adata_luminal_a.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = []

for group in groups:
    top_genes.extend(result['names'][group][:n_top_genes])

top_genes = pd.unique(top_genes)
print("Top genes:")
print(top_genes)

# === Subset data ===
adata_top_genes = adata_luminal_a[:, top_genes]

# === Set font size globally ===
sc.set_figure_params(scanpy=True, fontsize=8, dpi=150)

# === Plot heatmap without showing ===
sc.pl.heatmap(
    adata_top_genes,
    var_names=top_genes,
    groupby="ESR1_positive",
    use_raw=False,
    cmap="viridis",
    swap_axes=True,
    dendrogram=False,
    figsize=(8.27, 11.69 / 2),
    show=False
)

# === Fix axis labels ===
ax = plt.gcf().axes[0]
ax.set_title("Top 20 DEGs by ESR1 Status", fontsize=9)

# Replace default labels
xticklabels = ax.get_xticklabels()
new_labels = []

for label in xticklabels:
    val = label.get_text().strip()
    if val == '0':
        new_labels.append("ER NEGATIVE")
    elif val == '1':
        new_labels.append("ER POSITIVE")
    else:
        new_labels.append(val)

# Apply horizontal labels
ax.set_xticklabels(new_labels, rotation=0, ha='center', fontsize=8)

# === Save before showing ===
plt.savefig("top20_deg_heatmap.pdf", bbox_inches="tight")

# === Show the plot ===
plt.show()

# === Final message ===
print("Saved heatmap as 'top20_deg_heatmap.pdf'")


###ER positive cells- Immune Infiltration###
# Subset to Luminal A cells
adata_luminal_a = adata[adata.obs['C.TYPE'] == 'Luminal A'].copy()

# Full Luminal A gene signature list
gene_list = [
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

# Path to the Excel files
high_file = r"path\high.xlsx"
low_file = r"path\low.xlsx"

# Load sample IDs from each Excel file
# Assuming the IDs are in the first column, no header
high_df = pd.read_excel(high_file, header=None)
low_df = pd.read_excel(low_file, header=None)

# Convert to Python lists of strings
High = high_df.iloc[:, 0].astype(str).tolist()
Low = low_df.iloc[:, 0].astype(str).tolist()

print("High group sample count:", len(high))
print("Low group sample count:", len(low))
print("First few high samples:", high[:5])
print("First few low samples:", low[:5])

selected_pp = High + Low
adata_sub = adata_luminal_a[adata_luminal_a.obs['P_P'].isin(selected_pp)].copy()
existing_genes = [g for g in gene_list if g in adata_sub.var_names]

# Expression matrix
X = adata_sub[:, existing_genes].X
if not isinstance(X, np.ndarray):
    X = X.toarray()

expr_df = pd.DataFrame(X, index=adata_sub.obs_names, columns=existing_genes)
expr_df["P_P"] = adata_sub.obs["P_P"].values

# Fraction > 0.5 per sample
fraction_df = expr_df.groupby("P_P")[existing_genes].apply(lambda df: (df > 0.5).sum() / len(df))
signature_fraction = fraction_df.mean(axis=1).reset_index()
signature_fraction.columns = ["P_P", "LuminalA_signature_fraction"]
signature_fraction["Immune_Group"] = signature_fraction["P_P"].apply(lambda x: "High" if x in high else "Low")

# Add simplified BIOPSY_GROUP: 'B' for BREAST, 'M' for MET
biopsy_map = adata_sub.obs.set_index("P_P")["BIOPSY_GROUP"].to_dict()
simple_biopsy = {k: ('B' if str(v).strip().upper() == 'BREAST' else 'M') for k, v in biopsy_map.items()}
signature_fraction["BIOPSY_GROUP"] = signature_fraction["P_P"].map(simple_biopsy)

# Mann–Whitney test
g1 = signature_fraction[signature_fraction["Immune_Group"] == "High"]["LuminalA_signature_fraction"]
g2 = signature_fraction[signature_fraction["Immune_Group"] == "Low"]["LuminalA_signature_fraction"]
stat, p_value = mannwhitneyu(g1, g2)
p_label = "p < 2.2×10⁻¹⁶" if p_value < 2.2e-16 else f"p = {p_value:.1e}"

# Plotting function
def show_and_save_plot(size, filename=None):
    plt.figure(figsize=size)
    order = ["High", "Low"]
    palette = ["#CAB2D6", "#A6CEE3"]

    sns.boxplot(data=signature_fraction, x="Immune_Group", y="LuminalA_signature_fraction",
                order=order, palette=palette, fliersize=0)

    biopsy_palette = sns.color_palette("Set1", n_colors=signature_fraction["BIOPSY_GROUP"].nunique())
    biopsy_lut = dict(zip(signature_fraction["BIOPSY_GROUP"].unique(), biopsy_palette))

    sns.stripplot(data=signature_fraction, x="Immune_Group", y="LuminalA_signature_fraction",
                  hue="BIOPSY_GROUP", dodge=False, jitter=True, alpha=0.6, size=3,
                  order=order, palette=biopsy_lut)

    plt.legend(title=None, fontsize=6, bbox_to_anchor=(1.05, 1), loc='upper left')
    title_text = "ER+ Signature Fraction\nLow vs. High Immune Infiltration Samples\n" + p_label
    plt.title(title_text, fontsize=9)
    plt.xlabel("Infiltration Group", fontsize=8)
    plt.ylabel("Mean Fraction Expression", fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()

    if filename:
        plt.savefig(filename, format='pdf', bbox_inches='tight')

    plt.show()

# Show and save plots
show_and_save_plot((8.27 / 3.5, 11.69 / 4.2), "LuminalA_infiltration_small.pdf")
show_and_save_plot((8.27, 11.69), "LuminalA_infiltration_full.pdf")


# === ER-negative signature plot ===

# Define ER-negative gene list
er_negative_genes = [
    "PTMA", "FTH1", "TMSB10", "FAU", "TMSB4X", "EEF1A1", "UBA52", "FTL", "NACA",
    "H3F3A", "SERF2", "EIF1", "TPT1", "GAPDH", "MYL6", "EEF1D", "MIF", "B2M",
    "OAZ1", "ATP5F1E", "COX7C", "H3F3B", "HSPB1", "FOS", "COX6C", "CALM2",
    "PFDN5", "RACK1", "SRP14", "PPDPF", "TMA7", "BTF3", "JUN", "CHCHD2", "CFL1",
    "S100A11", "HSP90AB1", "UBB", "NDUFA4", "COX6A1", "TOMM7", "UQCR11", "PFN1"]

# Filter existing genes
existing_genes_erneg = [g for g in er_negative_genes if g in adata_sub.var_names]

# Expression matrix
X_erneg = adata_sub[:, existing_genes_erneg].X
if not isinstance(X_erneg, np.ndarray):
    X_erneg = X_erneg.toarray()

expr_df_erneg = pd.DataFrame(X_erneg, index=adata_sub.obs_names, columns=existing_genes_erneg)
expr_df_erneg["P_P"] = adata_sub.obs["P_P"].values

# Calculate fraction > 0.5
fraction_df_erneg = expr_df_erneg.groupby("P_P")[existing_genes_erneg].apply(lambda df: (df > 0.5).sum() / len(df))
signature_fraction_erneg = fraction_df_erneg.mean(axis=1).reset_index()
signature_fraction_erneg.columns = ["P_P", "ERneg_signature_fraction"]
signature_fraction_erneg["Immune_Group"] = signature_fraction_erneg["P_P"].apply(lambda x: "High" if x in high else "Low")

# Map BIOPSY_GROUP
signature_fraction_erneg["BIOPSY_GROUP"] = signature_fraction_erneg["P_P"].map(simple_biopsy)

# Mann–Whitney U test
g1_erneg = signature_fraction_erneg[signature_fraction_erneg["Immune_Group"] == "High"]["ERneg_signature_fraction"]
g2_erneg = signature_fraction_erneg[signature_fraction_erneg["Immune_Group"] == "Low"]["ERneg_signature_fraction"]
stat_erneg, pval_erneg = mannwhitneyu(g1_erneg, g2_erneg)
p_label_erneg = "p < 2.2×10⁻¹⁶" if pval_erneg < 2.2e-16 else f"p = {pval_erneg:.1e}"

# Plotting function
def plot_erneg(size, filename=None):
    plt.figure(figsize=size)
    order = ["High", "Low"]
    palette = ["#CAB2D6", "#A6CEE3"]

    sns.boxplot(data=signature_fraction_erneg, x="Immune_Group", y="ERneg_signature_fraction",
                order=order, palette=palette, fliersize=0)

    biopsy_palette = sns.color_palette("Set1", n_colors=signature_fraction_erneg["BIOPSY_GROUP"].nunique())
    biopsy_lut = dict(zip(signature_fraction_erneg["BIOPSY_GROUP"].unique(), biopsy_palette))

    sns.stripplot(data=signature_fraction_erneg, x="Immune_Group", y="ERneg_signature_fraction",
                  hue="BIOPSY_GROUP", dodge=False, jitter=True, alpha=0.6, size=3,
                  order=order, palette=biopsy_lut)

    plt.legend(title=None, fontsize=6, bbox_to_anchor=(1.05, 1), loc='upper left')
    title = "ER– Signature Fraction\nLow vs. High Immune Infiltration Samples\n" + p_label_erneg
    plt.title(title, fontsize=9)
    plt.xlabel("Infiltration Group", fontsize=8)
    plt.ylabel("Mean Fraction Expression", fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()

    if filename:
        plt.savefig(filename, format='pdf', bbox_inches='tight')
    plt.show()

# Save and show ER-negative plots
plot_erneg((8.27 / 3.5, 11.69 / 4.2), "ERneg_signature_small.pdf")
plot_erneg((8.27, 11.69), "ERneg_signature_fullpdf")


###ER negative signiture graphs###
# Function to calculate per-gene fraction
def calculate_fraction_df(adata, gene, threshold=0.5):
    gene_expr = adata[:, gene].X
    if not isinstance(gene_expr, np.ndarray):
        gene_expr = gene_expr.toarray().flatten()
    else:
        gene_expr = gene_expr.flatten()

    adata.obs[f"{gene}_expression"] = gene_expr
    adata.obs[f"{gene}_positive"] = adata.obs[f"{gene}_expression"] > threshold

    df = (
        adata.obs
        .groupby("P_P")
        .agg(
            total_cells=(f"{gene}_positive", 'size'),
            gene_positive=(f"{gene}_positive", 'sum'),
            ESR1_positive=('ESR1_positive', 'first')
        )
        .reset_index()
    )
    df[f"{gene}_fraction"] = (df["gene_positive"] / df["total_cells"]) * 100
    df["Gene"] = gene
    df["ESR1_label"] = df["ESR1_positive"].map({"POSITIVE": "+", "NEGATIVE": "–"})
    return df[["P_P", "ESR1_label", f"{gene}_fraction", "Gene", "ESR1_positive"]].rename(columns={f"{gene}_fraction": "Fraction"})

# Genes
genes = ["MYC", "GAPDH", "LDHA"]
dfs = [calculate_fraction_df(adata_luminal_a.copy(), gene) for gene in genes]
combined_df = pd.concat(dfs, ignore_index=True)

# Add Group for x-axis
combined_df["Group"] = combined_df["Gene"] + "_" + combined_df["ESR1_label"]

# Order and palette
group_order = [
    "MYC_+", "MYC_–",
    "GAPDH_+", "GAPDH_–",
    "LDHA_+", "LDHA_–"
]
palette = {group: "#A6CEE3" if "+" in group else "#FDBF6F" for group in group_order}

# P-values
p_values = {}
for gene in genes:
    pos = combined_df[(combined_df["Gene"] == gene) & (combined_df["ESR1_positive"] == "POSITIVE")]["Fraction"]
    neg = combined_df[(combined_df["Gene"] == gene) & (combined_df["ESR1_positive"] == "NEGATIVE")]["Fraction"]
    stat, p = mannwhitneyu(pos, neg)
    if p < 2.2e-16:
        p_values[gene] = "p < 2.2×10⁻¹⁶"
    else:
        base, exp = f"{p:.1e}".split("e")
        superscript = exp.translate(str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))
        p_values[gene] = f"p = {float(base):.0f}×10{superscript}"

# Plot
plt.figure(figsize=(8.27 / 2.5, 11.69 / 4.2))

sns.boxplot(
    data=combined_df,
    x="Group",
    y="Fraction",
    order=group_order,
    palette=palette,
    showfliers=False
)
sns.stripplot(
    data=combined_df,
    x="Group",
    y="Fraction",
    order=group_order,
    color="black",
    size=3,
    jitter=True,
    alpha=0.6
)

# Y-axis limit
plt.ylim(0, 120)

# Add p-values above
x_positions = [0.5, 2.5, 4.5]
y_pval = 110
for xpos, gene in zip(x_positions, genes):
    plt.text(xpos, y_pval, p_values[gene], ha="center", fontsize=8, color="black")

# Add gene names below plot area
y_gene_label = -12
for i, gene in enumerate(genes):
    x_center = i * 2 + 0.5
    plt.text(x_center, y_gene_label, gene, ha="center", va="top", fontsize=8)

# Adjust ticks
plt.xticks(
    ticks=[0, 1, 2, 3, 4, 5],
    labels=["+", "–", "+", "–", "+", "–"],
    fontsize=8
)

plt.title("MYC, GAPDH, and LDHA\nFraction Expression by ER Status", fontsize=9)
plt.xlabel("")
plt.ylabel("Fraction Expression (%)", fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()

# Save
outdir = r"path"
os.makedirs(outdir, exist_ok=True)
plt.savefig(os.path.join(outdir, "ER negative genes.pdf"), format='pdf', bbox_inches='tight')
plt.show()


###ER positive signiture graphs###
# === Function to calculate % of pathway genes expressed > threshold ===
def calculate_pathway_fraction_df(adata, gene_list, pathway_name, threshold=0.5):
    adata = adata.copy()
    existing_genes = [g for g in gene_list if g in adata.var_names]

    if len(existing_genes) == 0:
        raise ValueError(f"No genes from {pathway_name} found in adata.")

    expr = adata[:, existing_genes].X
    if not isinstance(expr, np.ndarray):
        expr = expr.toarray()

    positive_matrix = (expr > threshold).astype(int)
    mean_positive_fraction = positive_matrix.mean(axis=1) * 100  # Percent

    adata.obs[f"{pathway_name}_fraction"] = mean_positive_fraction

    df = (
        adata.obs
        .groupby("P_P")
        .agg(
            Fraction=(f"{pathway_name}_fraction", 'mean'),
            ESR1_positive=('ESR1_positive', 'first')
        )
        .reset_index()
    )
    df["Gene"] = pathway_name
    df["ESR1_label"] = df["ESR1_positive"].map({"POSITIVE": "+", "NEGATIVE": "–"})
    return df[["P_P", "ESR1_label", "Fraction", "Gene", "ESR1_positive"]]

# === Pathway gene sets ===
her_genes = ["ERBB2", "ERBB3", "ERBB4"]
akt_genes = ["AKT1", "AKT2", "AKT3", "MTOR", "PIK3CA"]
mapk_genes = ["KRAS", "BRAF", "CERK"]

# === Calculate fractions ===
df_her = calculate_pathway_fraction_df(adata_luminal_a, her_genes, "HER")
df_akt = calculate_pathway_fraction_df(adata_luminal_a, akt_genes, "AKT/MTOR")
df_mapk = calculate_pathway_fraction_df(adata_luminal_a, mapk_genes, "MAPK")

# === Combine ===
combined_df = pd.concat([df_her, df_akt, df_mapk], ignore_index=True)
combined_df["Group"] = combined_df["Gene"] + "_" + combined_df["ESR1_label"]

# === Order and colors ===
group_order = [
    "HER_+", "HER_–",
    "AKT/MTOR_+", "AKT/MTOR_–",
    "MAPK_+", "MAPK_–"
]
palette = {group: "#A6CEE3" if "+" in group else "#FDBF6F" for group in group_order}

# === P-values ===
p_values = {}
for gene in ["HER", "AKT/MTOR", "MAPK"]:
    pos = combined_df[(combined_df["Gene"] == gene) & (combined_df["ESR1_positive"] == "POSITIVE")]["Fraction"]
    neg = combined_df[(combined_df["Gene"] == gene) & (combined_df["ESR1_positive"] == "NEGATIVE")]["Fraction"]
    stat, p = mannwhitneyu(pos, neg)
    if p < 2.2e-16:
        p_values[gene] = "p < 2.2×10⁻¹⁶"
    else:
        base, exp = f"{p:.1e}".split("e")
        superscript = exp.translate(str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))
        p_values[gene] = f"p = {float(base):.0f}×10{superscript}"

# === Plot ===
plt.figure(figsize=(8.27 / 2.5, 11.69 / 4.2))

sns.boxplot(
    data=combined_df,
    x="Group",
    y="Fraction",
    order=group_order,
    palette=palette,
    showfliers=False
)
sns.stripplot(
    data=combined_df,
    x="Group",
    y="Fraction",
    order=group_order,
    color="black",
    size=3,
    jitter=True,
    alpha=0.6
)

# === Y-axis ===
plt.ylim(0, 100)

# Add gene labels and p-values
x_positions = [0.5, 2.5, 4.5]
gene_labels = ["HER", "AKT/MTOR", "MAPK"]
y_pval = 90
y_label = -10

for xpos, gene in zip(x_positions, gene_labels):
    plt.text(xpos, y_pval, p_values[gene], ha="center", fontsize=8)
    plt.text(xpos, y_label, gene, ha="center", va="top", fontsize=8)

# X-axis +/-
plt.xticks(
    ticks=[0, 1, 2, 3, 4, 5],
    labels=["+", "–", "+", "–", "+", "–"],
    fontsize=8
)

# Labels and layout
plt.title("HER, AKT/MTOR, and MAPK Pathways\n Mean Fraction Expression by ER Status", fontsize=9)
plt.xlabel("")
plt.ylabel("Mean Fraction Expression", fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()

# === Save ===
outdir = r"path"
os.makedirs(outdir, exist_ok=True)
plt.savefig(os.path.join(outdir, "ER positive genes.pdf"), format='pdf', bbox_inches='tight')
plt.show()

