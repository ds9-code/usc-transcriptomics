# -*- coding: utf-8 -*-
"""scRNA-SCPCA-LGG-Integrated.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1dn4bs-AfUk98Uu8xdFWqpG5yfftmxnzh
"""

!pip install scanpy

!pip install scikit-misc

import csv
import gzip
import os
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
import anndata as ad

from google.colab import drive
drive.mount('/content/drive')

def pp(dir_path, sample, grade):
    print("Processing " + dir_path)
    adata = sc.read_10x_mtx(path=dir_path, var_names='gene_ids', cache=True)
    adata.obs['Sample'] = grade + sample
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_genes(adata, min_cells=3)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]

    return adata

out = []
for dir in os.listdir('/content/drive/MyDrive/Transcriptomics/SCPC-2/'):
    out.append(pp('/content/drive/MyDrive/Transcriptomics/SCPC-2/' + dir, dir, 'LGG-'))

adata = sc.concat(out)

adata.X

adata.obs_names_make_unique()

adata

adata.obs

adata.write_h5ad('/content/drive/MyDrive/Transcriptomics/SCPC/lgg-combined.h5ad')

adata.shape

adata.X

"""# Preprocessing"""

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

adata.var.sort_values('n_cells_by_counts')

sc.pp.filter_genes(adata, min_cells=1000)

adata.var.sort_values('n_cells_by_counts')

adata.obs.sort_values('n_genes_by_counts')

sc.pl.highest_expr_genes(adata, n_top=20)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True)

upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
#upper_lim = 3000

upper_lim

adata = adata[adata.obs.n_genes_by_counts < upper_lim]

adata.obs

adata

"""# Normalization"""

adata.X.sum(axis = 1)

sc.pp.normalize_total(adata, target_sum=1e4) #normalize every cell to 10,000 UMI

adata.X.sum(axis = 1)

sc.pp.log1p(adata) #change to log counts

adata.X.sum(axis = 1)

adata.raw = adata

adata.raw.X

"""# Clustering"""

sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')

adata.var

sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts'])

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)

sc.pp.neighbors(adata, n_pcs = 30)

sc.tl.umap(adata)

sc.pl.umap(adata)

import locale
def getpreferredencoding(do_setlocale = True):
    return "UTF-8"
locale.getpreferredencoding = getpreferredencoding

!pip install leidenalg

sc.tl.leiden(adata, resolution = 0.5)

adata.obs

sorted_by_num_cells = adata.var['n_cells_by_counts'].sort_values(ascending=False)
sorted_by_num_cells

sc.pl.umap(adata, color=['leiden'])

sc.pl.umap(adata, color=["ENSG00000167996", "ENSG00000278996", "ENSG00000111640"])

sc.tl.rank_genes_groups(adata, "leiden", method="t-test")

sc.pl.rank_genes_groups(adata, n_genes=10)

sc.pl.rank_genes_groups_heatmap(adata, groups="5", n_genes=10, groupby="leiden")

# Obtain the top features for each cluster
result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
top_features = {}
n_top_genes = 10  # desired number of top genes per cluster
for group in groups:
    top_features[group] = result["names"][group][:n_top_genes]

# Print the top features for each cluster
for group, features in top_features.items():
    print(f"Cluster {group} top features:")
    for feature in features:
        print(feature)
    print()

# Access the marker genes results from rank_genes_groups
marker_genes = adata.uns['rank_genes_groups']

# Iterate over each group and print the marker genes
for group in marker_genes['names'].dtype.names:
    print(f"Group: {group}")
    print(marker_genes['names'][group][:10])  # Print the top 10 marker genes
    print("\n")

dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
dat.to_csv("/content/drive/MyDrive/Transcriptomics/SCPC/scanpy_integrated_lgg_nov30.csv")

!pip install openpyxl

!python /content/drive/MyDrive/Transcriptomics/SCSA/SCSA.py -d /content/drive/MyDrive/Transcriptomics/SCSA/whole_v2.db -i /content/drive/MyDrive/Transcriptomics/SCPC/scanpy_integrated_lgg_nov30.csv -s scanpy -g Human -f1.5 -p 0.01 -o /content/drive/MyDrive/Transcriptomics/SCPC/lgg_results_nov30.txt -m txt

# Read the text document into a DataFrame
df = pd.read_csv('/content/drive/MyDrive/Transcriptomics/SCPC/lgg_results_nov30.txt', sep='\t')

df

df['Z-score'] = df['Z-score'].fillna(0)

# Group the data by "Cluster" and find the cell type with the highest Z-score in each group
highest_zscores = df.groupby('Cluster')['Z-score'].idxmax()

# Extract the corresponding cell types for the highest Z-scores
cell_types_with_highest_zscores = df.loc[highest_zscores, 'Cell Type'].tolist()

cell_types_with_highest_zscores

cluster_to_cell_type = {cluster_num: cell_type for cluster_num, cell_type in enumerate(cell_types_with_highest_zscores)}
cluster_to_cell_type

adata.obs['leiden']

adata.obs['leiden'] = adata.obs['leiden'].astype(int)

adata.obs['cell_types'] = adata.obs['leiden'].map(cluster_to_cell_type)

adata.obs

sc.pl.umap(adata, color=['leiden', 'cell_types'])

from matplotlib.pyplot import rc_context
with rc_context({'figure.figsize': (8,8)}):
    sc.pl.umap(adata, color = ['cell_types'], frameon = False, s = 5, legend_loc = 'on data',
              legend_fontsize=12, legend_fontoutline=2)