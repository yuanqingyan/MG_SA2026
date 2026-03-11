import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys
import subprocess
import fnmatch
import glob
import bin2cell as b2c
import cv2
import sc_utils
import tensorflow
import celltypist
from celltypist import models
import anndata as ad
import squidpy as sq
import cellcharter as cc
import scvi
from lightning.pytorch import seed_everything
from copy import copy
reds = copy(plt.cm.Reds)

seed_everything(12345)
scvi.settings.seed = 12345

#sampleName='12-thymoma'#'6-Thymus'

os.makedirs("./output/B2C/stardist", exist_ok=True)
outPath="./output/B2C/stardist"

B2CFolderPath='./Apr2025_H5ADRDS/B2C_H5AD'
pathRaw='./VisiumHD-set2'
SaveFolderPath="./CellTypist_H5AD"
FigOut="./output/QC_Fig"

modelPath="./Projects/ManuscriptData/IntegrateThymus/data/H5AD"

os.chdir(FigOut)

sampleIn=['6-Thymus','17-Thymus']
temp6=sc.read_h5ad(f"{SaveFolderPath}/6_thymus_addNiche.h5ad")
temp17=sc.read_h5ad(f"{SaveFolderPath}/17_thymus_addNiche.h5ad")

temp6.obs['FNich2']= temp6.obs['F_Niche'].apply(lambda x: 'Medulla_Stroma' if x == 'Medulla_B' else x)
temp17.obs['FNich2']=temp17.obs['F_Niche']

allAda=[temp6,temp17]

fdata = ad.concat(allAda, 
                  axis=0, 
                  merge='same', 
                  pairwise=True, 
                  index_unique='_')
fdata



##CellCharter’s spatial clustering
import anndata as ad
import squidpy as sq
import cellcharter as cc

fdata.uns['spatial_fov'] = {s: {} for s in fdata.obs['SampleName'].unique()}
fdata.obs['sample'] = pd.Categorical(fdata.obs['SampleName'])
fdata

fdata.layers["counts"] = fdata.X.copy()
adata=fdata[fdata.obs.FNich2=="Medulla_Stroma"]
adata=adata.copy()

scvi.model.SCVI.setup_anndata(
  adata, 
  layer="counts", 
  batch_key='sample',
)

model = scvi.model.SCVI(adata)
model.train(early_stopping=True, 
            enable_progress_bar=True)

adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
adata.obsm['spatial_fov']=adata.obsm['spatial_cropped150_buffer']
sq.gr.spatial_neighbors(adata, 
                        library_key='sample', 
                        coord_type='generic', 
                        delaunay=True, 
                        spatial_key='spatial_fov', 
                        percentile=99)

cc.gr.aggregate_neighbors(adata, 
                          n_layers=3, 
                          use_rep='X_scVI', 
                          out_key='X_cellcharter', 
                          sample_key='sample')

autok = cc.tl.ClusterAutoK(
  n_clusters=(2,20), 
  max_runs=10,
  convergence_tol=0.001
)
autok.fit(adata, use_rep='X_cellcharter')
cc.pl.autok_stability(autok)
adata.obs['cluster_cellcharter'] = autok.predict(adata, use_rep='X_cellcharter',k=15)

selected_genesF=["CD79A","MS4A1","CXCL13",'CXCR4',
                 'BCL6','ICOS','CXCR5','CD83','PRDM1',##GC
                 "FN1","COL1A1",'FAS','MKI67']


sc.pl.matrixplot(
  adata,
  var_names=selected_genesF, 
  groupby='cluster_cellcharter',
  dendrogram=False,
  cmap="Reds",
  standard_scale="var",
  colorbar_title="Mean expression",
  show=False
)

adata.obs['Germinal']=adata.obs['cluster_cellcharter'].apply(lambda x: 'GC' if x in [1,10] else 'noGC')

sc.pl.matrixplot(
  adata,
  var_names=selected_genesF, 
  groupby='Germinal',
  dendrogram=False,
  cmap="Reds",
  standard_scale="var",
  colorbar_title="Mean expression",
  show=False
)

adata.write(f'{SaveFolderPath}/Merg6_17_withGC.h5ad')



