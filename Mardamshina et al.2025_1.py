#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Proteomics integration + TNBC analysis

Inputs:
- Expression Excel (MOESM11): sheet S9A/S9B/S9C etc.-Found in the supplementary information of the paper
  * columns = regions (e.g., BC25R7)
  * last column = gene symbols
- Metadata Excel (MOESM3): sheet S1A with Patient column (e.g., BC25)

Outputs (in out_dir):
- Clean matrix CSV (genes x regions)
- AnnData .h5ad (regions x genes) with metadata merged
- TNBC subset .h5ad (optional)
- TACSTD2 histogram (TNBC)
- TACSTD2 vs signature correlation plot PDF
- Quartile barplot (optional)
"""

from __future__ import annotations

import os
import re
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


# ============================================================
# 0) USER SETTINGS
# ============================================================
expr_file = r"file location\41467_2025_65477_MOESM11_ESM.xlsx"
expr_sheet = "S9B"  # <--- choose the sheet you want: S9A / S9B / S9C etc.

meta_file = r"file location\41467_2025_65477_MOESM3_ESM.xlsx"
meta_sheet = "S1A"
meta_patient_col = "Patient"         # MUST exist in metadata
meta_tnbc_col = "Pathologist DS"     # used for TNBC filtering
tnbc_value = "TN"                    # TNBC label

out_dir = r"folder location\Proteomics"
os.makedirs(out_dir, exist_ok=True)

# Output names
prefix = f"proteomics_{expr_sheet}"
out_clean_csv = os.path.join(out_dir, f"{prefix}_clean_matrix_genes_x_regions.csv")
out_h5ad = os.path.join(out_dir, f"{prefix}_with_{meta_sheet}_metadata.h5ad")
out_h5ad_tn = os.path.join(out_dir, f"{prefix}_TNBC_only.h5ad")

# Correlation plot output
out_corr_pdf = os.path.join(out_dir, "TACSTD2_signature_correlation.pdf")
out_hist_pdf = os.path.join(out_dir, "TACSTD2_histogram_TNBC.pdf")
out_quartile_pdf = os.path.join(out_dir, "Signature_by_TACSTD2_quartiles.pdf")

# Gene to correlate against signature (TROP2)
trop2_candidates = ["TACSTD2", "TROP2"]

# ============================================================
# 1) HELPERS
# ============================================================

def make_index_unique_numbered(idx) -> np.ndarray:
    """Make duplicate gene names unique by numbering in order: TP53, TP53_2, TP53_3, ..."""
    s = pd.Series(idx, dtype=str).astype(str).str.strip()
    s = s.replace({"nan": np.nan, "None": np.nan, "": np.nan})
    counts = s.groupby(s).cumcount()
    unique = s.where(counts == 0, s + "_" + (counts + 1).astype(str))
    return unique.values

def extract_patient(sample_id: str) -> str:
    """
    Extract patient from region ID.
    Example: BC25R7 -> BC25 (split at first 'R')
    """
    s = str(sample_id).strip()
    return s.split("R", 1)[0] if "R" in s else s

def make_h5ad_safe_column_names(cols) -> list[str]:
    """
    Make obs column names HDF5-safe (no /) and unique.
    """
    new_cols = []
    seen = set()
    for c in cols:
        c2 = str(c).replace("/", "_")
        c2 = re.sub(r"\s+", " ", c2).strip()
        c2 = c2.replace("\n", " ").replace("\r", " ")
        c2 = c2.replace(":", "_").replace("\\", "_")

        base = c2
        k = 2
        while c2 in seen:
            c2 = f"{base}_{k}"
            k += 1

        seen.add(c2)
        new_cols.append(c2)
    return new_cols

def to_dense_1d(x) -> np.ndarray:
    """Convert sparse or array-like to flat dense 1D numpy array."""
    if hasattr(x, "toarray"):
        return x.toarray().ravel()
    return np.asarray(x).ravel()

def safe_get_gene_vector(adata_obj: ad.AnnData, gene: str) -> np.ndarray:
    """Fetch gene vector from .X (preferred) or .raw if needed."""
    if gene in adata_obj.var_names:
        idx = np.where(adata_obj.var_names == gene)[0][0]
        return to_dense_1d(adata_obj.X[:, idx])
    if adata_obj.raw is not None and gene in adata_obj.raw.var_names:
        idx = np.where(adata_obj.raw.var_names == gene)[0][0]
        return to_dense_1d(adata_obj.raw.X[:, idx])
    raise KeyError(f"{gene} not found in adata.var_names or adata.raw.var_names")


# ============================================================
# 2) LOAD EXPRESSION (genes x regions) FROM EXCEL
# ============================================================
print(f"Loading expression sheet: {expr_sheet}")
df_raw = pd.read_excel(expr_file, sheet_name=expr_sheet)
print("Raw shape:", df_raw.shape)

# last column contains gene symbols
gene_col = df_raw.columns[-1]

genes = (
    df_raw[gene_col]
    .astype(str)
    .str.strip()
    .replace({"nan": np.nan, "None": np.nan, "": np.nan})
)

# keep only rows with valid gene names
df = df_raw.loc[genes.notna()].copy()
genes = genes.loc[genes.notna()].values

# set unique gene index
df.index = make_index_unique_numbered(genes)

# numeric-only conversion
df_num = df.apply(pd.to_numeric, errors="coerce")
df_num = df_num.dropna(axis=1, how="all").dropna(axis=0, how="all")

print("Final shape (genes x regions):", df_num.shape)
print(df_num.iloc[:5, :5])

# save clean matrix
df_num.to_csv(out_clean_csv)
print("Saved clean matrix:", out_clean_csv)


# ============================================================
# 3) BUILD AnnData (regions x genes) + MERGE METADATA
# ============================================================
# AnnData: obs = regions, var = genes
X = df_num.T.values
obs_names = df_num.columns.astype(str)   # region IDs (BCxxRyy)
var_names = df_num.index.astype(str)     # gene names

# load metadata
print(f"\nLoading metadata sheet: {meta_sheet}")
meta = pd.read_excel(meta_file, sheet_name=meta_sheet)

if meta_patient_col not in meta.columns:
    raise ValueError(
        f"Metadata sheet '{meta_sheet}' missing required '{meta_patient_col}'. "
        f"Columns: {list(meta.columns)}"
    )

meta[meta_patient_col] = meta[meta_patient_col].astype(str).str.strip()

# build obs table
obs = pd.DataFrame(index=obs_names)
obs["Patient"] = [extract_patient(x) for x in obs.index]

# deduplicate metadata per patient (keep first)
meta_dedup = meta.drop_duplicates(subset=[meta_patient_col], keep="first").copy()

# merge
obs = obs.merge(
    meta_dedup,
    how="left",
    left_on="Patient",
    right_on=meta_patient_col,
    suffixes=("", "_meta"),
)
obs.index = obs_names

# match summary
matched = obs[meta_patient_col].notna().sum() if meta_patient_col in obs.columns else obs[meta_patient_col].notna().sum()
print(f"Regions: {obs.shape[0]}")
print(f"Matched to metadata by Patient: {matched} / {obs.shape[0]}")

# list unmatched
if meta_patient_col in obs.columns:
    unmatched = obs.loc[obs[meta_patient_col].isna(), "Patient"].unique()
    if len(unmatched) > 0:
        print("Unmatched Patient IDs (first 20):", unmatched[:20])

# build AnnData
adata = ad.AnnData(
    X=X,
    obs=obs,
    var=pd.DataFrame(index=var_names),
)
adata.obs_names.name = "Region"
adata.var_names.name = "Gene"

# make obs columns HDF5-safe
adata.obs.columns = make_h5ad_safe_column_names(adata.obs.columns)

print("\nAnnData created:")
print(adata)
print("X shape (regions x genes):", adata.X.shape)

# save
adata.write_h5ad(out_h5ad)
print("Saved AnnData:", out_h5ad)


# ============================================================
# 4) TNBC SUBSET
# ============================================================
if meta_tnbc_col not in adata.obs.columns:
    print(f"\nWARNING: '{meta_tnbc_col}' not found in adata.obs. Skipping TNBC subset.")
    adata_tn = None
else:
    adata_tn = adata[adata.obs[meta_tnbc_col].astype(str).str.strip() == tnbc_value].copy()
    print("\nTNBC subset:")
    print(adata_tn)
    adata_tn.write_h5ad(out_h5ad_tn)
    print("Saved TNBC AnnData:", out_h5ad_tn)


# ============================================================
# 5) TACSTD2 HISTOGRAM (TNBC)
# ============================================================
if adata_tn is not None:
    gene_found = None
    for g in trop2_candidates:
        if g in adata_tn.var_names or (adata_tn.raw is not None and g in adata_tn.raw.var_names):
            gene_found = g
            break

    if gene_found is None:
        raise ValueError("TACSTD2/TROP2 not found in TNBC AnnData.")

    trop = safe_get_gene_vector(adata_tn, gene_found)

    plt.figure(figsize=(6, 4))
    plt.hist(trop[np.isfinite(trop)], bins=40)
    plt.title(f"Protein intensity histogram: {gene_found} (TNBC)")
    plt.xlabel("Protein intensity")
    plt.ylabel("Number of regions")
    plt.tight_layout()
    plt.savefig(out_hist_pdf, bbox_inches="tight")
    plt.close()
    print("Saved histogram:", out_hist_pdf)


# ============================================================
# 6) SIGNATURE vs TACSTD2 (TNBC)
# ============================================================
# Your TP21-like gene list (as provided)
gene_list = [
    "S100A9","SOX9","HSPH1","CD81","TMEM45A","RARRES1","SLPI","LCN2","HSPA6",
    "CP","S100A6","S100A8","CHI3L1","BAG3","IER5","LAMB3","PHLDA2","IL18",
    "NME2","GPRC5A","PMEPA1","HSPA1A","GABRE","CRYAB","MGP","HSP90AA1",
    "S100A11","ANXA1","PTN","HSPD1","CD59","ABL2","EGR4","S100A10","MUC20",
    "S100A4","ARPC1B","ITGB8","STEAP4","TPT1","DNAJB1","RDH10","IDS","HSPE1",
    "CSTB","KYNU","S100P","CLDN1","HSPA1B","CDK2AP1","KLF5","C15orf48","LAMP2",
    "DCBLD1","EGR3","UBC","PRRG4","ENC1","FN1","MAF","SDC4","LATS2","ZFAND2A",
    "SCPEP1","BRI3","TKT","CHI3L2","BCAP31","TNFRSF12A","SNCG","TGFBI","ANXA2",
    "TGM2","LTF","SFN","CAST","TNNC2","CLIC3","CMTM6","ATF3","CAMK2N1","TMSB10",
    "RCAN1","DEFB1","SERPINA1","DSC2","CD63","MAFK","RPL17","SSPN","SLC7A5",
    "CLTB","P4HB","SDCBP2","NQO1","ITM2B","CARD16","MGST1","SLC25A37","SSR4"
]

if adata_tn is not None:
    # choose gene symbol for TACSTD2/TROP2
    gene_found = None
    for g in trop2_candidates:
        if g in adata_tn.var_names:
            gene_found = g
            break
        if adata_tn.raw is not None and g in adata_tn.raw.var_names:
            gene_found = g
            break
    if gene_found is None:
        raise ValueError("TACSTD2/TROP2 not found in TNBC AnnData.")

    trop = safe_get_gene_vector(adata_tn, gene_found)

    # available signature genes
    available = [g for g in gene_list if g in adata_tn.var_names]
    missing = [g for g in gene_list if g not in adata_tn.var_names]

    print(f"\nSignature genes available: {len(available)}")
    if len(missing) > 0:
        print(f"Signature genes missing (n={len(missing)}): first 20 -> {missing[:20]}")

    if len(available) == 0:
        raise ValueError("None of the signature genes were found in adata_tn.var_names.")

    X_sig = adata_tn[:, available].X
    if hasattr(X_sig, "toarray"):
        X_sig = X_sig.toarray()
    else:
        X_sig = np.asarray(X_sig, dtype=float)

    signature = np.nanmean(X_sig, axis=1)

    df_corr = pd.DataFrame({"TACSTD2": trop, "signature": signature})
    df_corr = df_corr[np.isfinite(df_corr["TACSTD2"]) & np.isfinite(df_corr["signature"])].copy()

    # correlations
    pear_r, pear_p = pearsonr(df_corr["TACSTD2"], df_corr["signature"])
    spear_r, spear_p = spearmanr(df_corr["TACSTD2"], df_corr["signature"])

    print(f"\nCorrelation results (TNBC):")
    print(f"Pearson  r = {pear_r:.3f}, p = {pear_p:.3e}")
    print(f"Spearman r = {spear_r:.3f}, p = {spear_p:.3e}")

    # plot with your A4-fraction size + font size 8
    plt.rcParams["font.size"] = 8
    fig, ax = plt.subplots(figsize=(8.27 / 3.5, 11.69 / 4.7))  # your sizing

    ax.scatter(df_corr["TACSTD2"], df_corr["signature"], alpha=0.4)

    # regression line
    m, b = np.polyfit(df_corr["TACSTD2"], df_corr["signature"], 1)
    xs = np.linspace(df_corr["TACSTD2"].min(), df_corr["TACSTD2"].max(), 100)
    ax.plot(xs, m * xs + b)

    ax.set_xlabel(f"{gene_found} protein intensity", fontsize=8)
    ax.set_ylabel("Mean TP21 signature intensity", fontsize=8)
    ax.set_title("Proteomic TACSTD2 correlates with TP21 signature (TNBC)", fontsize=9)

    # Spearman annotation (bottom-right)
    mantissa, exponent = f"{spear_p:.1e}".split("e")
    exponent = int(exponent)
    ax.text(
        0.98, 0.02,
        f"$r = {spear_r:.2f}$\n$ p = {mantissa}\\times 10^{{{exponent}}}$",
        ha="right", va="bottom",
        transform=ax.transAxes,
        fontsize=8,
    )

    plt.savefig(out_corr_pdf, format="pdf", bbox_inches="tight")
    plt.close()
    print("Saved correlation PDF:", out_corr_pdf)

    # OPTIONAL: signature by TACSTD2 quartiles
    # (kept because you used this view)
    try:
        df_q = df_corr.copy()
        df_q["quartile"] = pd.qcut(
            df_q["TACSTD2"],
            q=4,
            labels=["Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"],
            duplicates="drop"
        )
        mean_by_q = df_q.groupby("quartile")["signature"].mean()
        print("\nMean signature by TACSTD2 quartile:")
        print(mean_by_q)

        plt.figure(figsize=(6, 4))
        plt.bar(mean_by_q.index.astype(str), mean_by_q.values)
        plt.ylabel("Mean TP21 signature intensity")
        plt.xlabel(f"{gene_found} quartile")
        plt.title("TP21 signature across TACSTD2 quartiles (TNBC)")
        plt.tight_layout()
        plt.savefig(out_quartile_pdf, bbox_inches="tight")
        plt.close()
        print("Saved quartile plot:", out_quartile_pdf)
    except Exception as e:
        print("Quartile plot skipped due to:", repr(e))

print("\nDONE.")
