#Convert h5ad to csv for INFERCNV.R

import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp

folder_path = r'C:\Users\user\Desktop\New folder (2)'#path to folder contaning h5ad files for conversion to csv

for filename in os.listdir(folder_path):
    if not filename.endswith('.h5ad'):
        continue

    in_path = os.path.join(folder_path, filename)
    print(f"Reading: {in_path}")
    adata = sc.read_h5ad(in_path)

    # Choose the matrix you want (here: .X)
    X = adata.X

    # Build a DataFrame correctly (handles sparse or dense)
    if sp.issparse(X):
        X = X.tocsr()  # or .tocsc(), either is fine
        df = pd.DataFrame.sparse.from_spmatrix(
            X, index=adata.obs_names, columns=adata.var_names
        )
    elif isinstance(X, np.ndarray) and X.ndim == 2:
        df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    else:
        raise ValueError(f"Unsupported matrix type for {filename}: {type(X)}")

    # Save CSV next to the .h5ad
    csv_filename = os.path.splitext(filename)[0] + '.csv'
    csv_path = os.path.join(folder_path, csv_filename)
    df.to_csv(csv_path)
    print(f"Converted {filename} → {csv_filename} (shape {df.shape[0]}×{df.shape[1]})")

