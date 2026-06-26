import scanpy as sc
import os

# Define the folder containing your .h5ad files
folder_path = r"TUMOR datasets folder path"

# Get all .h5ad files in the folder
h5ad_files = [f for f in os.listdir(folder_path) if f.endswith('.h5ad')]

# Load the .h5ad files and skip those that do not have 'obs'
adata_list = []
for file in h5ad_files:
    try:
        adata = sc.read(os.path.join(folder_path, file))
        # Append to the list if 'obs' attribute exists
        if adata.obs is not None:
            adata_list.append(adata)
    except KeyError as e:
        print(f"Skipping {file} due to error: {e}")

# Concatenate the AnnData objects using outer join to keep all variables
if len(adata_list) > 0:
    adata = sc.concat(adata_list, join='outer')

    # Ensure all columns in .obs and .var are strings
    adata.obs.columns = adata.obs.columns.astype(str)
    adata.var.columns = adata.var.columns.astype(str)

    # Ensure that all keys in .uns are strings
    adata.uns = {str(key): value for key, value in adata.uns.items()}

    # Check for any problematic columns in obs and var and convert them to strings
    for column in adata.obs.columns:
        if not isinstance(adata.obs[column].dtype, str):
            adata.obs[column] = adata.obs[column].astype(str)

    for column in adata.var.columns:
        if not isinstance(adata.var[column].dtype, str):
            adata.var[column] = adata.var[column].astype(str)

    # Save the merged AnnData object
    adata.write('merged_adata_outer_fixed.h5ad')

    print(adata)
else:
    print("No valid .h5ad files were loaded.")

#Checking the data
adata.obs

#
adata.var

#
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="paper")

#
sc.pl.highly_variable_genes(adata)

#
adata.raw = adata

#
adata = adata[:, adata.var.highly_variable]


#
sc.tl.pca(adata)


#
sc.pl.pca(adata, color="GATA3")


#
sc.pl.pca_variance_ratio(adata, log=True)


#
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)


#
sc.tl.umap(adata)

# Replace "NA" with "UNKNOWN" in adata.obs["C.TYPE"]
adata.obs["C.TYPE"].replace("NA", "UNKNOWN", inplace=True)


# Ensure 'C.TYPE' recognizes 'UNKNOWN' as a category
adata.obs["C.TYPE"] = adata.obs["C.TYPE"].astype('category')
if 'UNKNOWN' not in adata.obs["C.TYPE"].cat.categories:
    adata.obs["C.TYPE"].cat.add_categories(['UNKNOWN'], inplace=True)

# Now plot the UMAP
sc.pl.umap(adata, color=["paper", "C.TYPE"])
sc.external.pp.bbknn(adata, batch_key='paper')  # running bbknn paper=data_set
sc.tl.umap(adata)
sc.pl.umap(adata,color=["paper", "C.TYPE"])

#
adata.write('tumor.h5ad')

###cNMF analysis was done following this step and is added in the"cnmf.zip"folder

