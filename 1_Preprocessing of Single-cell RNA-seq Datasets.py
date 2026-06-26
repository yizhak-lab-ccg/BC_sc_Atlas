
#This code includes the preprocessing of each single-cell RNA-seq datasets- after downloading the GSE file and processing it when needed to h5ad file

import matplotlib.pyplot as plt
import os
import scanpy as sc
import anndata
import pandas as pd
import scrublet as scr
import scipy.sparse as sp
import scrublet as scr
import numpy as np

# 
save_file = r'path\16_DEMO_INPUT.h5ad' #Input- the DEMO file or another h5ad file
adata = sc.read_h5ad(save_file)
adata

# === Use counts from adata.X ===
X = adata.X

# Keep it sparse CSC and float32
if not sp.issparse(X):
    X = sp.csc_matrix(X)
else:
    X = X.tocsc()
X = X.astype(np.float32)
# === Run Scrublet ===
scrub = scr.Scrublet(counts_matrix=X, expected_doublet_rate=0.05, random_state=0)
doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2,    
    min_cells=3,     
    n_prin_comps=30
)

# === Store in .obs ===
adata.obs["doublet_score"] = doublet_scores
adata.obs["predicted_doublet"] = predicted_doublets.astype(bool)


# === 1) Define fixed cutoff and create a boolean flag ===
CUTOFF = 0.05
adata.obs["doublet_over_005"] = adata.obs["doublet_score"] > CUTOFF

# === 2) Summary BEFORE filtering ===
n_total = adata.n_obs
n_over = int(adata.obs["doublet_over_005"].sum())
print(f"Total cells BEFORE: {n_total:,}")
print(f"Cells with doublet_score > {CUTOFF}: {n_over:,} ({n_over/n_total:.2%})")
print(f"Max score BEFORE: {float(adata.obs['doublet_score'].max()):.4f}")

# === 3) Filter by the fixed cutoff ===
adata = adata[~adata.obs["doublet_over_005"]].copy()

# === 4) Verify AFTER filtering ===
n_after = adata.n_obs
print(f"Total cells AFTER:  {n_after:,}")
print("Any score > cutoff left?:", bool((adata.obs["doublet_score"] > CUTOFF).any()))
print(f"Max score AFTER: {float(adata.obs['doublet_score'].max()):.4f}")

#
sc.pl.highest_expr_genes(adata, n_top=20)

#
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


#
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


#
threshold = np.percentile(adata.obs["pct_counts_mt"], 90)  # 90th percentile cutoff
adata = adata[adata.obs["pct_counts_mt"] <= threshold].copy()

print(f"✅ Remaining cells: {adata.n_obs}")

# Violin plots after filtering

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


#
#manual checks for each dataset based on the graph
#adata = adata[adata.obs.n_genes_by_counts < #manual threshold if needed, :]
#adata = adata[adata.obs.pct_counts_mt <#manual threshold if needed, :]


#Check the data
adata

#Saving count data
adata.layers["counts"] = adata.X.copy()

#
sc.pp.normalize_total(adata, target_sum=1e4)

# 
sc.pp.log1p(adata)

#
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'patient')

print("Highly variable genes intersection: %d"%sum(adata.var.highly_variable_intersection))

print("Number of batches where gene is variable:")
print(adata.var.highly_variable_nbatches.value_counts())

var_genes_batch = adata.var.highly_variable_nbatches > 0


#
sc.pl.highly_variable_genes(adata)


#
adata.raw = adata


#
adata = adata[:, adata.var.highly_variable]


# 
sc.tl.pca(adata, svd_solver="arpack")


#
sc.pp.neighbors(adata)


#
sc.tl.umap(adata)
sc.pl.umap(adata,color=['patient'])


#
sc.external.pp.bbknn(adata, batch_key='patient')

#
sc.tl.umap(adata)
sc.pl.umap(adata,color=['patient'])

#
adata

#
save_file = r'17_DEMO_after_preprocessing.h5ad'
adata.write_h5ad(save_file)

