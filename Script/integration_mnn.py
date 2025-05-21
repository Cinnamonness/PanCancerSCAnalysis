import scanpy as sc
import anndata as ad
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import mnnpy
warnings.filterwarnings('ignore')

# Настройка картинок
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)
sns.set(font_scale=1)
sns.set_style('ticks')

# Загрузка данных
adata_combined = sc.read_h5ad('/home/karitskaya/karitskaya/own_data/adata_combined.h5ad')
# Выделение высоковариабельных генов
sc.pp.highly_variable_genes(adata_combined, n_top_genes=5000, batch_key="dataset")
adata_combined = adata_combined[:, adata_combined.var.highly_variable]
# PCA и UMAP
sc.tl.pca(adata_combined, n_comps=50)
sc.pp.neighbors(adata_combined, use_rep="X_pca")
sc.tl.leiden(adata_combined)
#sc.tl.umap(adata_combined)
# Выделение в отдельные batches
batch_categories = adata_combined.obs['Dataset'].unique()
datas = [adata_combined[adata_combined.obs['Dataset'] == batch] for batch in batch_categories]
# Perform MNN
corrected = mnnpy.mnn_correct(*datas,
                              batch_key='Dataset')
#Объединение датасета
adata = corrected[0]
# PCA и UMAP после MNN
sc.pp.scale(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
# Сохранение файла
adata.write('all_after_mnn.h5ad')