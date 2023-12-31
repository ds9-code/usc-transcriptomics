# -*- coding: utf-8 -*-
"""Integrated DGE diffxpy.ipynb
Author: DS
"""

!pip show dask

!pip install scanpy

!pip install diffxpy

!pip install dask==2021.4.0

from google.colab import drive
drive.mount('/content/drive')

import diffxpy.api as de
import csv
import gzip
import os
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
import anndata as ad

adata = sc.read_h5ad('/content/drive/MyDrive/Transcriptomics/SCPC/pre-dge-adata.h5ad')

subset = adata[adata.obs['cell_types'].isin(['Microglial cell', 'Natural killer T (NKT) cell'])].copy()

subset = subset.raw.to_adata()

subset.X = subset.X.toarray()

len(subset.var)

subset.obs

res = de.test.wald(data=subset,
             formula_loc= '~ 1 + cell_types',
             factor_loc_totest='cell_types'
                  )

subset.obs.cell_types.unique()

dedf = res.summary().sort_values('log2fc', ascending = False).reset_index(drop = True)
dedf

most_up = dedf.iloc[0].gene
i = np.where(subset.var_names == most_up)[0][0]

a = subset[subset.obs.cell_types == 'Microglial cell'].X[:, i]
b = subset[subset.obs.cell_types == 'Natural killer T (NKT) cell'].X[:, i]
print(f"{most_up} expression:")
print(f"Microglial cell: {a.mean()}")
print(f"Natural killer T (NKT) cell: {b.mean()}")

dedf['log2fc'] = dedf['log2fc']*-1
dedf = dedf.sort_values('log2fc', ascending = False).reset_index(drop = True)
dedf

dedf = dedf[(dedf.qval < 0.05) & (abs(dedf.log2fc) > .5)]
dedf

dedf = dedf[dedf['mean'] > 0.15]
dedf

genes_to_show = dedf[-25:].gene.tolist() + dedf[:25].gene.tolist() #top 25 and bottom 25 from sorted df

sc.pl.heatmap(subset, genes_to_show, groupby='cell_types', swap_axes=True)

df_subset = subset.obs
df_subset.to_csv('/content/drive/MyDrive/Transcriptomics/SCPC/subset_dge.csv')