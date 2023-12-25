# -*- coding: utf-8 -*-
"""Integrated-LGG-HGG.ipynb
Author: DS
"""

!pip install scanpy

!pip install scikit-misc

!pip install diffxpy

import locale
def getpreferredencoding(do_setlocale = True):
    return "UTF-8"
locale.getpreferredencoding = getpreferredencoding
!pip install leidenalg

!pip install openpyxl

from google.colab import drive
drive.mount('/content/drive')

import csv
import gzip
import os
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
import anndata as ad

hgg_adata = sc.read_h5ad('/content/drive/MyDrive/Transcriptomics/SCPC/hgg-combined.h5ad')
lgg_adata = sc.read_h5ad('/content/drive/MyDrive/Transcriptomics/SCPC/lgg-combined.h5ad')

adatas = {"LGG": lgg_adata, "HGG": hgg_adata}
adata = ad.concat(adatas, label="dataset_name", join = "outer")

adata.obs_names_make_unique()

adata.obs

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

adata

adata.obs

sc.pp.filter_genes(adata, min_cells = 100)

adata.obs.groupby('Sample').count()

adata.obs.groupby('dataset_name').count()

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata

sc.pl.highest_expr_genes(adata, n_top=20)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True)

sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')

sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var.highly_variable]

adata

sc.pp.regress_out(adata, ['total_counts'])

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)

sc.pp.neighbors(adata, n_pcs = 30)

sc.tl.umap(adata)

sc.pl.umap(adata)

sc.tl.leiden(adata, resolution = 0.5)
adata.obs

sc.pl.umap(adata, color = ['leiden', 'dataset_name'], frameon = False)

sc.tl.leiden(adata, resolution = 1)
sc.tl.rank_genes_groups(adata, 'leiden')

markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers

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
dat.to_csv("/content/drive/MyDrive/Transcriptomics/SCPC/integrated_lgg_hgg.csv")

!python /content/drive/MyDrive/Transcriptomics/SCSA/SCSA.py -d /content/drive/MyDrive/Transcriptomics/SCSA/whole_v2.db -i /content/drive/MyDrive/Transcriptomics/SCPC/integrated_lgg_hgg.csv -s scanpy -g Human -f1.5 -p 0.01 -o /content/drive/MyDrive/Transcriptomics/SCPC/integrated_cell_type.txt -m txt

# Read the text document into a DataFrame
df = pd.read_csv('/content/drive/MyDrive/Transcriptomics/SCPC/integrated_cell_type.txt', sep='\t')

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

df_adata = adata.obs
df_adata.to_csv('/content/drive/MyDrive/Transcriptomics/SCPC/final_adata.csv')

num_tot_cells = adata.obs.groupby(['Sample']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.n_genes_by_counts))
num_tot_cells

cell_type_counts = adata.obs.groupby(['Sample', 'dataset_name', 'cell_types']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]
cell_type_counts

cell_type_counts['total_cells'] = cell_type_counts.Sample.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.n_genes_by_counts / cell_type_counts.total_cells
cell_type_counts

cell_type_counts.to_csv('/content/drive/MyDrive/Transcriptomics/SCPC/cell_type_counts.csv')

import matplotlib.pyplot as plt
plt.figure(figsize = (10,4))
ax = sns.boxplot(data = cell_type_counts, x = 'cell_types', y = 'frequency', hue = 'dataset_name')
plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
plt.show()

sc.pl.umap(adata, color=['cell_types', 'dataset_name'], frameon = True, legend_loc = "on data")

adata.write_h5ad('/content/drive/MyDrive/Transcriptomics/SCPC/pre-dge-adata.h5ad')