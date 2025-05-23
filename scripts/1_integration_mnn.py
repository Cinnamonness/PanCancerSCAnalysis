# Import libraries
import scanpy as sc
import anndata as ad
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import mnnpy

warnings.filterwarnings("ignore")

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)
sns.set(font_scale=1)
sns.set_style("ticks")

# Load the Data
adata_combined = sc.read_h5ad(
    "/home/karitskaya/karitskaya/own_data/adata_combined.h5ad"
)
# Select highly variable genes
sc.pp.highly_variable_genes(adata_combined, n_top_genes=5000, batch_key="dataset")
adata_combined = adata_combined[:, adata_combined.var.highly_variable]
# PCA and UMAP
sc.tl.pca(adata_combined, n_comps=50)
sc.pp.neighbors(adata_combined, use_rep="X_pca")
sc.tl.leiden(adata_combined)
# sc.tl.umap(adata_combined)
# Select batches
batch_categories = adata_combined.obs["Dataset"].unique()
datas = [
    adata_combined[adata_combined.obs["Dataset"] == batch] for batch in batch_categories
]
# Perform MNN
corrected = mnnpy.mnn_correct(*datas, batch_key="Dataset")
# Concat the Data
adata = corrected[0]
# PCA and UMAP after MNN
sc.pp.scale(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
# Save the file
adata.write("all_after_mnn.h5ad")
