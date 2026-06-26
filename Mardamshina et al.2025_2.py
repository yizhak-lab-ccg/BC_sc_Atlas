import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr


# ============================================================
# 0) Paths / settings
# ============================================================
out_dir = r"folder location\Proteomics"
os.makedirs(out_dir, exist_ok=True)

# --- Load the AnnData you already created before ---
adata_path = os.path.join(out_dir, "proteomics_S9B_with_S1A_metadata.h5ad")  # <-- change if your name differs
adata = sc.read_h5ad(adata_path)
print("Loaded:", adata_path)
print(adata)

# --- ESR1 extraction source sheet ---
expr_file = r"folder location\Proteomics\41467_2025_65477_MOESM11_ESM.xlsx"
esr_sheet = "S9A"   # S9A for ESR1 NaNs

# --- Stroma / aSMA sheet ---
stromal_file = r"folder location\Proteomics\41467_2025_65477_MOESM8_ESM.xlsx"
stromal_sheet = "S6D"
stromal_region_col = "Region"
stromal_value_col = "aSMA_percent_PROI"

# --- Luminal definition in metadata ---
path_ds_col = "Pathologist DS"
luminal_values = ["ERPR", "ER"]

# --- Plot sizing (your journal style)
plt.rcParams["font.size"] = 8
FIG_W = 8.27 / 3.5
FIG_H = 11.69 / 4.2


# ============================================================
# 1) Add ESR1 protein + missing indicator to adata.obs
# ============================================================
df = pd.read_excel(expr_file, sheet_name=esr_sheet)
gene_col = df.columns[-1]  # last column has gene names

df_genes = df.set_index(gene_col)

if "ESR1" not in df_genes.index:
    raise ValueError("ESR1 not found in the sheet used for ESR1 extraction (S9A).")

esr1_series = df_genes.loc["ESR1"]

# align to adata regions (obs_names should match Excel columns)
esr1_aligned = esr1_series.reindex(adata.obs_names)

adata.obs["ESR1_protein"] = pd.to_numeric(esr1_aligned, errors="coerce").values
adata.obs["ESR1_protein_is_na"] = adata.obs["ESR1_protein"].isna().astype(int)

# report
nan_cols = esr1_series.index[esr1_series.isna()].tolist()
print("\nColumns where ESR1 is NaN (from Excel):", len(nan_cols))
print("Example NaN columns:", nan_cols[:20])
print("\nESR1 status counts in adata:")
print(adata.obs["ESR1_protein_is_na"].value_counts(dropna=False))


# ============================================================
# 2) Subset Luminal regions (adata_la)
# ============================================================
if path_ds_col not in adata.obs.columns:
    raise ValueError(f"'{path_ds_col}' not found in adata.obs. Available columns: {list(adata.obs.columns)[:30]} ...")

adata_la = adata[adata.obs[path_ds_col].isin(luminal_values)].copy()
print("\nadata_la (Luminal) shape:", adata_la.shape)

# sanity
print("Luminal Pathologist DS counts:")
print(adata_la.obs[path_ds_col].value_counts(dropna=False))


# ============================================================
# 3) ER+ signature mean (Luminal) by ESR1 missing/present
# ============================================================
status_col = "ESR1_protein_is_na"

er_sig_genes = [
    "ESR1","TRPS1","THSD4","ERBB4","ANKRD30A","AFF3","PATJ","PBX1","LRBA","TTC6","MSI2",
    "FOXP1","RBM47","RERE","MED13L","ANKS1B","FTX","AGAP1","RORA","DACH1","NEBL","ASH1L",
    "APBB2","KIAA1217","ZFAND3","MAST4","MAGI1","SIPA1L3","MACF1","NFIA","IGF1R","NAALADL2",
    "LPP","AUTS2","TRIO","SHANK2","LDLRAD4","ZBTB20","MAGI3","BCL2","XIST","VPS13B","BMPR1B",
    "PVT1","TBC1D9","EXOC4","WWOX","SPIDR","TRAPPC9","WWP1","TBC1D5","MLPH","SRPK2","KDM4B",
    "MRTFB","PTK2","KAZN","HIVEP3","SMYD3","GALNT7","CERS6","ARID1B","NEDD4L","GSE1","PRLR",
    "RABEP1","NCOA2","ITPR2","GLI3","EVL","EIF4G3","NSMCE2","FAF1","IMMP2L","ADK","DOCK1",
    "JMJD1C","TCF12","KMT2C","SYTL2","USP34","GRHL2","SH3RF1","GOLGA4","SIK3","PTPRT","KAT6B",
    "FARP1","C1orf21","BIRC6","CASZ1","MYO6","PDXDC1","DNAH14","RAD51B","LIMCH1","SHROOM3",
    "ATXN1","FUT8","TULP4","PARD3","STAG1","ALCAM","MBD5","CACNA1D","ZNF609","ETV6","BRWD1",
    "TMCC1","RERG","FHIT","DNAJC1","TNRC6A","MIRLET7BHG","PPP6R3","RBM6","NF1","GATA3",
    "LDLRAD3","EHMT1","KIAA1324","FLNB","FBXL17","CELSR1","FAM172A","CUX1","UBR5","BCAS3",
    "ZMYND8","INPP4B","ARHGAP32","TTC3","KCNMA1","SAMD12","DENND1B","ANKRD17","PKP4","CTBP2",
    "ZCCHC7","KDM2A","DLEU2","TNRC18","ZNF91","SETD5","FBXL7","NIPBL","EXOC6B","DLG1",
    "BTBD9","GMDS-DT","ZNF638","FNDC3B","TANC2","KANSL1","EXT1","RALGPS2","BAZ2B","STRBP",
    "DTNB","ANK3","GTF2I","ERC1","FARS2","BTRC","MKLN1","SSBP2","VAV3","PLEKHA7","ARHGEF3",
    "GPHN","CSNK1A1","NCKAP5","ERBIN","AR","CBFA2T2","TBL1XR1","ARID4B","MAML3","POGZ","VWA8",
    "CLASP2","PAN3","ENAH","FAM214A","MAP7","MIR29B2CHG","IKBKB","MUC20-OT1","COP1","CADPS2",
    "C1QTNF3-AMACR","ATXN2","BPTF","NCOA1","PUM2","NSD1","TENT2","PTPRK","RALGAPA1","ZNF148"
]

genes_in_adata = [g for g in er_sig_genes if g in adata_la.var_names]
missing = [g for g in er_sig_genes if g not in adata_la.var_names]

print(f"\nER signature genes total: {len(er_sig_genes)}")
print(f"Found in proteomics: {len(genes_in_adata)}")
print(f"Missing: {len(missing)} (showing first 20): {missing[:20]}")

# Save present genes
pd.DataFrame({"Gene": genes_in_adata}).to_excel(
    os.path.join(out_dir, "ER_signature_genes_present_in_adata_la.xlsx"),
    index=False
)

if len(genes_in_adata) == 0:
    raise ValueError("None of the ER signature genes were found in adata_la.var_names.")

X = adata_la[:, genes_in_adata].X
if not isinstance(X, np.ndarray):
    X = X.toarray().astype(float)

adata_la.obs["ER_signature_mean"] = np.nanmean(X, axis=1)

plot_df = adata_la.obs[[status_col, "ER_signature_mean"]].copy()
plot_df = plot_df.dropna(subset=[status_col])  # keep defined ESR1 groups only
plot_df["ESR1_label"] = plot_df[status_col].map({0: "ER+ (ESR1 present)", 1: "ER− (ESR1 missing)"})

g0 = plot_df.loc[plot_df[status_col] == 0, "ER_signature_mean"].dropna()
g1 = plot_df.loc[plot_df[status_col] == 1, "ER_signature_mean"].dropna()

p_mwu = np.nan
if len(g0) > 0 and len(g1) > 0:
    _, p_mwu = mannwhitneyu(g0, g1, alternative="two-sided")

print("\nMWU p-value (ER signature mean):", p_mwu)

# Plot + save
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
sns.boxplot(data=plot_df, x="ESR1_label", y="ER_signature_mean", ax=ax)
sns.stripplot(data=plot_df, x="ESR1_label", y="ER_signature_mean", color="black", alpha=0.5, ax=ax)

ax.set_xlabel("ESR1 protein status")
ax.set_ylabel("Mean intensity of ER-related signature")
ax.set_title(f"ER signature intensity by ESR1 status (Luminal)\nMWU p = {p_mwu:.2e}" if not np.isnan(p_mwu) else
             "ER signature intensity by ESR1 status (Luminal)\nMWU p = NA")

plt.tight_layout()
out_pdf = os.path.join(out_dir, "ER_signature_boxplot_ESR1status_Luminal.pdf")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()
print("Saved:", out_pdf)


# ============================================================
# 4) Robustness analysis-
# ============================================================
X_all = X  # already dense matrix for genes_in_adata
status_vals = adata_la.obs[status_col].values

results = []
for k in range(10, len(genes_in_adata) + 1):
    sig_k = np.nanmean(X_all[:, :k], axis=1)
    df_k = pd.DataFrame({"sig": sig_k, "status": status_vals}).dropna(subset=["status"])
    group0 = df_k.loc[df_k["status"] == 0, "sig"].dropna()
    group1 = df_k.loc[df_k["status"] == 1, "sig"].dropna()
    if len(group0) > 0 and len(group1) > 0:
        _, p = mannwhitneyu(group0, group1, alternative="two-sided")
    else:
        p = np.nan
    results.append({"n_genes": k, "p_value": p})

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(out_dir, "ER_signature_bestK_MWU_results.csv"), index=False)

valid = results_df.dropna(subset=["p_value"])
if not valid.empty:
    best = valid.loc[valid["p_value"].idxmin()]
    best_n = int(best["n_genes"])
    best_p = float(best["p_value"])
    print(f"\nBest K by MWU: n_genes={best_n}, p={best_p:.3e}")

    sig_best = np.nanmean(X_all[:, :best_n], axis=1)
    df_best = pd.DataFrame({"sig": sig_best, "status": status_vals}).dropna(subset=["status"])
    df_best["ESR1_label"] = df_best["status"].map({0: "ER+ (ESR1 present)", 1: "ER− (ESR1 missing)"})

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    sns.boxplot(data=df_best, x="ESR1_label", y="sig", ax=ax)
    sns.stripplot(data=df_best, x="ESR1_label", y="sig", color="black", alpha=0.5, ax=ax)

    ax.set_xlabel("ESR1 protein status")
    ax.set_ylabel(f"Mean intensity of ER signature (first {best_n} genes)")
    ax.set_title(f"Best ER signature separation (Luminal)\nMWU p = {best_p:.2e}")

    plt.tight_layout()
    out_pdf = os.path.join(out_dir, "ER_signature_boxplot_bestK_Luminal.pdf")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print("Saved:", out_pdf)
else:
    print("No valid p-values computed for bestK search.")


# ============================================================
# 5) Add aSMA_percent_PROI from MOESM8 S6D to adata and adata_la
# ============================================================
df_stroma = pd.read_excel(stromal_file, sheet_name=stromal_sheet)

if stromal_region_col not in df_stroma.columns or stromal_value_col not in df_stroma.columns:
    raise ValueError(
        f"Missing required columns in {stromal_sheet}. "
        f"Need: {stromal_region_col}, {stromal_value_col}. Found: {list(df_stroma.columns)}"
    )

df_stroma_sub = df_stroma[[stromal_region_col, stromal_value_col]].copy()
df_stroma_sub = df_stroma_sub.set_index(stromal_region_col)

# add to full adata
adata.obs[stromal_value_col] = df_stroma_sub.reindex(adata.obs_names)[stromal_value_col].values

# add to luminal subset
adata_la.obs[stromal_value_col] = df_stroma_sub.reindex(adata_la.obs_names)[stromal_value_col].values

print("\naSMA matched (Luminal):",
      np.isfinite(pd.to_numeric(adata_la.obs[stromal_value_col], errors="coerce")).sum(),
      "/", adata_la.n_obs)


# ============================================================
# 6) Boxplot: aSMA_percent_PROI by ESR1 status (Luminal)
# ============================================================
plot_df2 = adata_la.obs[[status_col, stromal_value_col]].copy()
plot_df2[stromal_value_col] = pd.to_numeric(plot_df2[stromal_value_col], errors="coerce")
plot_df2 = plot_df2.dropna(subset=[status_col, stromal_value_col])
plot_df2["ESR1_label"] = plot_df2[status_col].map({0: "ER+", 1: "ER−"})

pos = plot_df2.loc[plot_df2["ESR1_label"] == "ER+", stromal_value_col]
neg = plot_df2.loc[plot_df2["ESR1_label"] == "ER−", stromal_value_col]

p_asma = np.nan
if len(pos) > 0 and len(neg) > 0:
    _, p_asma = mannwhitneyu(pos, neg, alternative="two-sided")

print("\nMWU p-value (aSMA by ESR1 status):", p_asma)

fig, ax = plt.subplots(figsize=(8.27 / 2.5, FIG_H))
palette = {"ER+": "#a6cee3", "ER−": "#fdbf6f"}

sns.boxplot(data=plot_df2, x="ESR1_label", y=stromal_value_col, order=["ER+", "ER−"], palette=palette, ax=ax)
sns.stripplot(data=plot_df2, x="ESR1_label", y=stromal_value_col, order=["ER+", "ER−"],
              color="black", alpha=0.5, ax=ax)

ax.set_xlabel("ER protein status")
ax.set_ylabel("Fibroblast percentage (aSMA_percent_PROI)")
ax.set_title(f"Fibroblast percentage: ER+ vs ER− (Luminal)\nMWU p = {p_asma:.2e}" if not np.isnan(p_asma)
             else "Fibroblast percentage: ER+ vs ER− (Luminal)\nMWU p = NA")

plt.tight_layout()
out_pdf = os.path.join(out_dir, "aSMA_vs_ESR1_boxplot_Luminal.pdf")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()
print("Saved:", out_pdf)


# ============================================================
# 7) Spearman correlation: ER signature vs aSMA
# ============================================================
aSMA = pd.to_numeric(adata_la.obs[stromal_value_col], errors="coerce").values

results_corr = []
for k in range(10, len(genes_in_adata) + 1):
    sig_k = np.nanmean(X_all[:, :k], axis=1)

    df_k = pd.DataFrame({"sig": sig_k, "aSMA": aSMA}).dropna()
    if len(df_k) < 5:
        r, p = np.nan, np.nan
    else:
        r, p = spearmanr(df_k["sig"], df_k["aSMA"])
    results_corr.append({"n_genes": k, "r": r, "p_value": p})

corr_df = pd.DataFrame(results_corr)
corr_df.to_csv(os.path.join(out_dir, "ER_signature_vs_aSMA_bestK_Spearman.csv"), index=False)

valid_corr = corr_df.dropna(subset=["p_value"])
if not valid_corr.empty:
    best = valid_corr.loc[valid_corr["p_value"].idxmin()]
    best_n = int(best["n_genes"])
    best_r = float(best["r"])
    best_p = float(best["p_value"])

    print(f"\nBest Spearman vs aSMA: n_genes={best_n}, r={best_r:.3f}, p={best_p:.3e}")

    sig_best = np.nanmean(X_all[:, :best_n], axis=1)
    plot_df3 = pd.DataFrame({"signature": sig_best, "aSMA": aSMA}).dropna()

    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    sns.regplot(data=plot_df3, x="signature", y="aSMA",
                scatter_kws=dict(alpha=0.6), ax=ax)

    ax.set_xlabel(f"Mean ER signature (first {best_n} genes)")
    ax.set_ylabel("aSMA_percent_PROI")
    ax.set_title("Spearman correlation (Luminal)")

    ax.text(
        0.05, 0.95,
        f"r = {best_r:.2f}\np = {best_p:.2e}",
        transform=ax.transAxes,
        ha="left", va="top"
    )

    plt.tight_layout()
    out_pdf = os.path.join(out_dir, "ER_signature_vs_aSMA_regplot_bestK_Luminal.pdf")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print("Saved:", out_pdf)
else:
    print("No valid Spearman p-values computed for bestK correlation search.")



