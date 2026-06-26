import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import scanpy as sc

# Define the file path
file_path = r"path to downloaded raw_counts txt file"

# Load the file into a DataFrame
df = pd.read_csv(file_path, sep="\t", index_col=0)

# Convert to an AnnData object (flipping obs and var)
adata = ad.AnnData(df)  # Genes should be in .var, cells in .obs

# Create a dictionary mapping cell prefixes to responses
response_map = {
    "001": "R", "002": "IR", "003": "AR", "004": "AR", "005": "IR", "006": "R",
    "007": "AR", "008": "AR", "009": "R", "010": "R", "011": "R", "012": "R",
    "013": "AR", "014": "IR", "015": "AR", "016": "AR", "017": "AR", "018": "AR",
    "019": "AR", "020": "AR", "021": "AR", "022": "AR", "023": "AR", "024": "AR",
    "025": "AR", "026": "AR", "029": "R", "030": "R", "031": "R", "032": "R",
    "033": "R", "034": "R", "035": "R"
}

adata
# Assign responses based on cell names
adata.obs["response"] = adata.obs.index.str[:3].map(response_map)

# Display the first few rows of adata.obs
print(adata.obs.head())

#Save the AnnData object
adata.write_h5ad(r"path\RNAseq_adata.h5ad")

##QC
# Filter genes expressed in fewer than 5 samples
min_cells = 5
filtered_genes = (adata.X > 0).sum(axis=0) >= min_cells
adata = adata[:, filtered_genes]


plt.hist(adata.X.sum(axis=1), bins=50)
plt.xlabel("Total Counts per Sample")
plt.ylabel("Number of Samples")
plt.show()

adata.layers["raw_counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

##ER POSITIVE DEG##

# List of genes
genes = [
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

# Ensure only genes present in adata are selected
genes_present = [gene for gene in genes if gene in adata.var_names]

# Extract mean expression of selected genes
mean_expression = pd.DataFrame(
    adata[:, genes_present].X.mean(axis=1),
    index=adata.obs.index,
    columns=["Mean Expression"]
)
mean_expression["Response"] = adata.obs["response"]

# Create a boxplot
plt.figure(figsize=(12, 6))
sns.boxplot(data=mean_expression, x="Response", y="Mean Expression", palette="Set2")

# Customize plot
plt.xlabel("Response Group")
plt.ylabel("Mean Gene Expression")
plt.title("Mean Expression of Selected Genes per Response Group")
plt.xticks(rotation=45)
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show plot
plt.show()

# List of genes
genes = [
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

# === Filter out AR group ===
adata = adata[~adata.obs["response"].isin(["AR"])].copy()

# Ensure only genes present in adata are selected
genes_present = [gene for gene in genes if gene in adata.var_names]

# Extract mean expression of selected genes
mean_expression = pd.DataFrame(
    adata[:, genes_present].X.mean(axis=1),
    index=adata.obs.index,
    columns=["Mean Expression"]
)
mean_expression["Response"] = adata.obs["response"]

# Perform Wilcoxon tests for each pair of response groups
response_groups = mean_expression["Response"].unique()
p_values = {}

for group1, group2 in combinations(response_groups, 2):
    group1_values = mean_expression[mean_expression["Response"] == group1]["Mean Expression"]
    group2_values = mean_expression[mean_expression["Response"] == group2]["Mean Expression"]
    
    if len(group1_values) > 1 and len(group2_values) > 1:
        stat, p_val = stats.mannwhitneyu(group1_values, group2_values, alternative="two-sided")
        p_values[(group1, group2)] = p_val

# Create a boxplot with p-value annotation
plt.figure(figsize=(12, 6))
ax = sns.boxplot(data=mean_expression, x="Response", y="Mean Expression", palette="Set2")

# Add p-values to the plot
y_max = mean_expression["Mean Expression"].max()
y_min = mean_expression["Mean Expression"].min()
y_range = y_max - y_min

for i, ((group1, group2), p_val) in enumerate(p_values.items()):
    x1, x2 = response_groups.tolist().index(group1), response_groups.tolist().index(group2)
    y = y_max + (i * 0.05 * y_range)
    plt.plot([x1, x1, x2, x2], [y, y + 0.02 * y_range, y + 0.02 * y_range, y], lw=1.5, color='black')
    plt.text((x1 + x2) / 2, y + 0.02 * y_range, f"P={p_val:.2e}", ha='center', va='bottom', fontsize=10, color='black')

# Customize plot
plt.xlabel("Response Group")
plt.ylabel("Mean Gene Expression")
plt.title("Mean Expression of Selected Genes per Response Group (Wilcoxon Test)")
plt.xticks(rotation=45)
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show plot
plt.show()

# Print p-values
print("\nPairwise Wilcoxon Test P-values:")
for (group1, group2), p_val in p_values.items():
    print(f"{group1} vs {group2}: P-value = {p_val:.4e}")


#ER NEGATIVE DEG

def format_pval(p):
    if p < 0.001:
        base = f"{p:.1e}".split("e")[0]
        exponent = int(f"{p:.1e}".split("e")[1])
        superscript_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
        return f"P = {base}×10{str(exponent).translate(superscript_map)}"
    else:
        return f"P = {p:.2f}".rstrip("0").rstrip(".")

# Filter adata to only include R and IR cells
adata_filtered = adata[adata.obs["response"].isin(["R", "AR"])].copy()

# List of genes for analysis
genes = [
    "PTMA", "FTH1", "TMSB10", "FAU", "TMSB4X", "EEF1A1", "UBA52", "FTL", "NACA", "H3F3A",
    "SERF2", "EIF1", "TPT1", "GAPDH", "MYL6", "EEF1D", "MIF", "B2M", "OAZ1", "ATP5F1E",
    "COX7C", "H3F3B", "HSPB1", "FOS", "COX6C", "CALM2", "PFDN5", "RACK1", "SRP14", "PPDPF",
    "TMA7", "BTF3", "JUN", "CHCHD2", "CFL1", "S100A11", "HSP90AB1", "UBB", "NDUFA4", "COX6A1",
    "TOMM7", "UQCR11", "PFN1", "PPIA", "COX6B1", "UBL5", "UQCRQ", "COX4I1", "UBC", "UQCRB"]

# Keep only genes present
genes_present = [gene for gene in genes if gene in adata_filtered.var_names]

# Calculate mean expression
mean_expression = pd.DataFrame(
    adata_filtered[:, genes_present].X.mean(axis=1),
    index=adata_filtered.obs.index,
    columns=["Mean Expression"]
)
mean_expression["Response"] = adata_filtered.obs["response"]

# Wilcoxon test
group_r = mean_expression[mean_expression["Response"] == "R"]["Mean Expression"]
group_ir = mean_expression[mean_expression["Response"] == "AR"]["Mean Expression"]
stat, p_val = stats.mannwhitneyu(group_r, group_ir, alternative="two-sided")
formatted_pval = format_pval(p_val)

# Function to plot and save
def plot_and_save(fig_width, fig_height, output_path):
    plt.figure(figsize=(fig_width, fig_height))
    ax = sns.boxplot(data=mean_expression, x="Response", y="Mean Expression", palette="Set2", fliersize=0)
    sns.stripplot(data=mean_expression, x="Response", y="Mean Expression",
                  color='black', size=3, jitter=True, dodge=True, alpha=0.7)

    # Y-axis scale
    y_max = mean_expression["Mean Expression"].max()
    y_min = mean_expression["Mean Expression"].min()
    y_range = y_max - y_min
    plt.ylim(y_min, y_max + 0.15 * y_range)

    # Annotate p-value
    x1, x2 = 0, 1
    y = y_max + 0.05 * y_range
    plt.plot([x1, x1, x2, x2], [y, y + 0.02 * y_range, y + 0.02 * y_range, y], lw=1.5, color='black')
    plt.text((x1 + x2) / 2, y + 0.025 * y_range, formatted_pval,
             ha='center', va='bottom', fontsize=8)

    # Aesthetics
    ax.set_xlabel("Response Group", fontsize=8)
    ax.set_ylabel("Mean Genes Expression", fontsize=8)
    ax.set_title("Mean ER− Signature: R vs AR", fontsize=9)
    ax.tick_params(axis='both', labelsize=8)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(output_path, format="pdf")
    plt.close()

# Save full-page A4
plot_and_save(8.27, 11.69, "ER_negative_signature_R_vs_AR_fullpage.pdf")

# Save smaller version
plot_and_save(8.27 / 2.5, 11.69 / 4.2, "ER_negative_signature_R_vs_AR_small.pdf")

# Print p-value
print(f"Wilcoxon P-value (R vs AR): {formatted_pval}")

