# Import libraries
import celltypist
import scanpy as sc
from celltypist import models
import pandas as pd
import anndata as ad
import numpy as np

# File with clusters
adata = sc.read_h5ad(
    "/home/karitskaya/karitskaya/own_data/adata_combined_with_clusters.h5ad"
)

# Load models
models.download_models(
    model=[
        "Cells_Intestinal_Tract.pkl",
        "Human_Colorectal_Cancer.pkl",
        "Fetal_Human_Pancreas.pkl",
        "Adult_Human_PancreaticIslet.pkl",
        "Healthy_Human_Liver.pkl",
        "Adult_Human_Skin.pkl",
        "Developing_Human_Gonads.pkl",
    ]
)

# PANCREAS

pancreas_clusters = adata[
    (adata.obs["clusters"] == "44")
    | (adata.obs["clusters"] == "70")
    | (adata.obs["clusters"] == "5")
    | (adata.obs["clusters"] == "8")
]

# Annotation with 'Adult_Human_PancreaticIslet.pkl'
# An object with consists of prediction results
predictions_pancreas = celltypist.annotate(
    pancreas_clusters,
    over_clustering="clusters",
    model="/home/ukrainskaya/.celltypist/data/models/Fetal_Human_Pancreas.pkl",
    majority_voting=True,
)
pancreas_clusters.obs.loc[
    pancreas_clusters.obs["Cancer type"] == "Pancreas Normal", "final_cell_annotation"
] = (
    pancreas_clusters.obs["clusters"].astype(str)
    + "_"
    + pancreas_clusters.obs["Cancer type"].astype(str)
    + "_"
    + predictions_pancreas.predicted_labels["majority_voting"].astype(str)
)
pancreas_clusters.write(
    "/home/ukrainskaya/files/data/processed_data/pancreas_clusters.h5ad"
)
