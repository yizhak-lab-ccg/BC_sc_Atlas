import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


# Define the file path
file_path = r"path to normalized log transformed tsv file"

# Load the TSV file into a Pandas DataFrame
data = pd.read_csv(file_path, sep='\t', index_col=0)

# Transpose the DataFrame so that rows are sample names and columns are genes
data_transposed = data.T

# Display basic info about the transposed dataset
print(data_transposed.info())
print(data_transposed.head())

# Filter samples where row names begin with specific prefixes
sample_prefixes = ('A1', 'A20', 'A11', 'A23', 'A30', 'A7', 'A15', 'A17', 'A26', 'A5')
data_filtered = data_transposed[data_transposed.index.str.startswith(sample_prefixes)]

# Display basic info about the transposed dataset
print(data_filtered.info())
print(data_filtered.head())

# Check if there is a column name that begins with 'PTPRC'
ptprc_exists = any(col.startswith('PTPRC') for col in data_transposed.columns)
print(f"Column name starting with 'PTPRC' exists: {ptprc_exists}")


##TP21 TP Genes##
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# Remove samples A1, A5, A26
data_filtered = data_filtered[~data_filtered.index.str.startswith(('A1', 'A5', 'A26'))]

# Extract relevant genes based on matching character names before '|'
gene_list = ["SOX9", "ITGB8", "FN1", "ANXA1", "SDC4", "CLDN1", "ZFAND2A", "BRI3",  
  "ARPC1B", "TNFRSF12A", "NQO1", "TGFBI", "PRRG4", "NME2", "LCN2", "LATS2",  
  "CSTB", "HSP90AA1", "ABL2", "CDK2AP1", "ITM2B", "MUC20", "ENC1", "PTN",  
  "DCBLD1", "TMSB10", "BCAP31", "BAG3", "ATF3", "CRYAB", "C15orf48", "TGM2",  
  "SERPINA1", "CLTB", "KYNU", "DEFB1"]

# Extract gene columns matching before '|'
gene_columns = [col for col in data_filtered.columns if any(col.split('|')[0] == gene for gene in gene_list)]
data_selected = data_filtered[gene_columns]

# Separate samples containing 'MET', 'Met', or 'met' in row names vs those that do not
met_samples = data_selected[data_selected.index.str.contains("MET|Met|met", regex=True)]
non_met_samples = data_selected[~data_selected.index.str.contains("MET|Met|met", regex=True)]

# Calculate mean expression for each sample
met_means = met_samples.mean(axis=1)
non_met_means = non_met_samples.mean(axis=1)

# Assign colors to sample prefixes using full A* category
unique_prefixes = sorted(set(sample.split('-')[0] for sample in met_means.index) | set(sample.split('-')[0] for sample in non_met_means.index))
colors = sns.color_palette("husl", len(unique_prefixes))
prefix_color_map = {prefix: colors[i] for i, prefix in enumerate(unique_prefixes)}

# Scatter plot with positions swapped
plt.figure(figsize=(12, 6))

# Plot Non-MET samples on the left
x_positions_non_met = list(range(len(non_met_means)))
for i, sample in enumerate(non_met_means.index):
    prefix = sample.split('-')[0]
    plt.scatter(x_positions_non_met[i], non_met_means.iloc[i], color=prefix_color_map[prefix],
                label=prefix if prefix not in plt.gca().get_legend_handles_labels()[1] else "")

# Plot MET samples on the right
x_positions_met = [i + len(non_met_means) + 5 for i in range(len(met_means))]
for i, sample in enumerate(met_means.index):
    prefix = sample.split('-')[0]
    plt.scatter(x_positions_met[i], met_means.iloc[i], color=prefix_color_map[prefix],
                label=prefix if prefix not in plt.gca().get_legend_handles_labels()[1] else "")

# Connect matching samples with lines
for prefix in unique_prefixes:
    non_met_match = [i for i, sample in enumerate(non_met_means.index) if sample.startswith(prefix)]
    met_match = [i + len(non_met_means) + 5 for i, sample in enumerate(met_means.index) if sample.startswith(prefix)]
    if met_match and non_met_match:
        for non_met_idx in non_met_match:
            for met_idx in met_match:
                plt.plot([x_positions_non_met[non_met_idx], x_positions_met[met_idx - len(non_met_means) - 5]],
                         [non_met_means.iloc[non_met_idx], met_means.iloc[met_idx - len(non_met_means) - 5]],
                         color='gray', linestyle='dotted', alpha=0.6)

# Calculate and add Mann-Whitney U test p-value
stat, p_value = mannwhitneyu(non_met_means, met_means)
# Format p-value with superscript notation
exponent = int(f"{p_value:.1e}".split('e')[1])
base = float(f"{p_value:.1e}".split('e')[0])
plt.text(len(non_met_means), max(max(non_met_means), max(met_means)), f'p-value: {base}×10$^{{{exponent}}}$',
         fontsize=10, ha='center', va='bottom')

plt.axvline(x=len(non_met_means) + 2.5, color='black', linestyle='--')  # Divider between Non-MET and MET
plt.xlabel("Samples")
plt.ylabel("Mean Gene Expression")
plt.title("Mean TP21 Gene Expression in Matched Metastasis vs Primary Samples")
plt.legend(title="Sample", fontsize='small', frameon=False, loc='best')
plt.xticks([])
plt.show()


#LAG3 & PDCD1
# Remove samples that are A1, A5, A26
data_filtered = data_filtered[~data_filtered.index.str.startswith(('A1', 'A5', 'A26'))]

# Extract relevant genes based on matching character names before '|'
gene_list = ["PDCD1","LAG3"]

# Extract gene columns matching before '|'
gene_columns = [col for col in data_filtered.columns if any(col.split('|')[0] == gene for gene in gene_list)]
data_selected = data_filtered[gene_columns]

# Separate samples containing 'MET', 'Met', or 'met' in row names vs those that do not
met_samples = data_selected[data_selected.index.str.contains("MET|Met|met", regex=True)]
non_met_samples = data_selected[~data_selected.index.str.contains("MET|Met|met", regex=True)]

# Calculate mean expression for each sample
met_means = met_samples.mean(axis=1)
non_met_means = non_met_samples.mean(axis=1)

# Assign colors to sample prefixes using full A* category
unique_prefixes = sorted(set(sample.split('-')[0] for sample in met_means.index) | set(sample.split('-')[0] for sample in non_met_means.index))
colors = sns.color_palette("husl", len(unique_prefixes))
prefix_color_map = {prefix: colors[i] for i, prefix in enumerate(unique_prefixes)}

# Scatter plot with positions swapped
plt.figure(figsize=(12, 6))

# Plot Non-MET samples on the left
x_positions_non_met = list(range(len(non_met_means)))
for i, sample in enumerate(non_met_means.index):
    prefix = sample.split('-')[0]
    plt.scatter(x_positions_non_met[i], non_met_means.iloc[i], color=prefix_color_map[prefix],
                label=prefix if prefix not in plt.gca().get_legend_handles_labels()[1] else "")

# Plot MET samples on the right
x_positions_met = [i + len(non_met_means) + 5 for i in range(len(met_means))]
for i, sample in enumerate(met_means.index):
    prefix = sample.split('-')[0]
    plt.scatter(x_positions_met[i], met_means.iloc[i], color=prefix_color_map[prefix],
                label=prefix if prefix not in plt.gca().get_legend_handles_labels()[1] else "")

# Connect matching samples with lines
for prefix in unique_prefixes:
    non_met_match = [i for i, sample in enumerate(non_met_means.index) if sample.startswith(prefix)]
    met_match = [i + len(non_met_means) + 5 for i, sample in enumerate(met_means.index) if sample.startswith(prefix)]
    if met_match and non_met_match:
        for non_met_idx in non_met_match:
            for met_idx in met_match:
                plt.plot([x_positions_non_met[non_met_idx], x_positions_met[met_idx - len(non_met_means) - 5]],
                         [non_met_means.iloc[non_met_idx], met_means.iloc[met_idx - len(non_met_means) - 5]],
                         color='gray', linestyle='dotted', alpha=0.6)

# Calculate and add Mann-Whitney U test p-value
stat, p_value = mannwhitneyu(non_met_means, met_means)
# Format p-value with superscript notation
exponent = int(f"{p_value:.1e}".split('e')[1])
base = float(f"{p_value:.1e}".split('e')[0])
plt.text(len(non_met_means), max(max(non_met_means), max(met_means)), f'p-value: {base}×10$^{{{exponent}}}$',
         fontsize=10, ha='center', va='bottom')

plt.axvline(x=len(non_met_means) + 2.5, color='black', linestyle='--')  # Divider between Non-MET and MET
plt.xlabel("Samples")
plt.ylabel("Mean Gene Expression")
plt.title("Mean TP21 Gene Expression in Matched Metastasis vs Primary Samples")
plt.legend(title="Sample", fontsize='small', frameon=False, loc='best')
plt.xticks([])
plt.show()

