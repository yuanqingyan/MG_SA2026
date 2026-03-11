#!/usr/bin/env python
# coding: utf-8


import scvi
import scanpy as sc
import anndata as ad
import torch
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
pd.set_option('display.max_rows', None)
from copy import copy
import os, sys
import sc_utils
reds = copy(plt.cm.Reds)
reds.set_under("lightgray")
scvi.settings.seed = 1

projectName="Thymus_Manu"

DataFolder="./Projects/ManuscriptData/IntegrateThymus/data/toAnndata"
SaveFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData"
tempOut="./SC/Thymus/output/scVI/output"

NUThymus= scvi.data.read_h5ad(f"{DataFolder}/NUThymus_fromVDJ.h5ad")
Yasumizu= scvi.data.read_h5ad(f"{DataFolder}/Yasumizu.h5ad")
Xin= scvi.data.read_h5ad(f"{DataFolder}/Xin.h5ad")

adata=ad.concat([NUThymus,Yasumizu,Xin],join="outer")
adata.obs.index=adata.obs.index.astype(str)+"_"+adata.obs.GSE.astype(str)
adata.obs_names_make_unique()

del NUThymus,Yasumizu,Xin#


pd.crosstab(adata.obs.batchKey.astype(str), adata.obs.GSE.astype(str),dropna=False)
adata.write(f"{SaveFolder}/adata_{projectName}_concatenate.h5ad")

sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw=adata

sc.pp.highly_variable_genes(
  adata,
  n_top_genes=3000,
  subset=True,
  flavor="seurat",
  batch_key="batchKey")
print("Complete sc.pp.highly_variable_genes!")

scvi.model.SCVI.setup_anndata(adata,
                              layer="counts",
                              batch_key="batchKey")

model = scvi.model.SCVI(adata, n_layers=2, dropout_rate=0.2)
model
model.train(max_epochs=400)

model.save(f"{SaveFolder}/{projectName}_adata_model")

latent= model.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=10e4)

#############################################################################################
#run PCA then generate UMAP plots
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)

adata.write(f"{tempOut}/{projectName}.gpu_scvi_umap.h5ad")

sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3)
adata.write(f"{tempOut}/{projectName}.gpu_scvi_umap_bathC.h5ad")

############################################################################################################
###########################    some plots --- before vs after scVI (Fig.S1A)     ###########################
############################################################################################################
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["figure.figsize"] = (10, 10)

adata_B=scvi.data.read_h5ad(f"{tempOut}/{projectName}.gpu_scvi_umap.h5ad")
adata=scvi.data.read_h5ad(f"{tempOut}/{projectName}.gpu_scvi_umap_bathC.h5ad")

fig, axs=plt.subplots(nrows=1, ncols=2, figsize=(18,4), constrained_layout=True)
sc.pl.umap(adata_B,color="GSE",ax=axs[0],frameon=False,add_outline=False,show=False,palette="tab20",size=2)
sc.pl.umap(adata,color="GSE",ax=axs[1],frameon=False,add_outline=False,show=False,palette="tab20",size=2)
plt.savefig(f"{tempOut}/{projectName}_umap.Before.AfterBC.pdf")

del adata_B
############################################################################################################
###########################          some  gene plots  (Fig.S1C)                 ###########################
############################################################################################################
sc.settings.vector_friendly = True

fig, axs=plt.subplots(nrows=2, ncols=2, figsize=(12,12), constrained_layout=True)
sc.pl.umap(adata,color=["EPCAM"],ax=axs[0,0],title="EPCAM",frameon=False, show=False,cmap=reds,size=3, vmin=0.001)
sc.pl.umap(adata,color=["PECAM1"],ax=axs[0,1],title="PECAM1",frameon=False, show=False,cmap=reds,size=3, vmin=0.001)
sc.pl.umap(adata,color=["PTPRC"],ax=axs[1,0],title="PTPRC",frameon=False, show=False,cmap=reds,size=3, vmin=0.001)
sc.pl.umap(adata,color=["COL1A1"],ax=axs[1,1],title="COL1A1",frameon=False, show=False,cmap=reds,size=3, vmin=0.001)
plt.savefig(f"{tempOut}/{projectName}_umap_4Markers.pdf")


sc.tl.leiden(adata, key_added="leiden_scVI", resolution=1)
adata.write(f"{SaveFolder}/{projectName}_scVI.res1.h5ad")

sc.set_figure_params(figsize=(12, 8))
sc.pl.umap(adata,color=["leiden_scVI"],frameon=False, show=False,palette="Set1", legend_loc="on data")
plt.savefig(f"{tempOut}/{projectName}_umap_leiden_scVI.pdf")
