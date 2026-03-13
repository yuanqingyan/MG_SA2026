import scvi
import scanpy as sc
import torch
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import sc_utils
import matplotlib as mpl
pd.set_option('display.max_rows', None)
from copy import copy
reds = copy(plt.cm.Reds)
reds.set_under("lightgrey")
import os
import sys
from scipy import io
import subprocess
import math
import celltypist
from celltypist import models

print("Last run with scvi-tools version:", scvi.__version__)

projectName="ThymomaB_May2025"

SaveFolder="./Thymus/h5ad"
tempOut="./Thymus/h5ad"

adata=scvi.data.read_h5ad("./Thymus/h5ad/ReferenceB.h5ad")

#####build celltypise model
adataOut = adata
io.mmwrite(f'{tempOut}/{projectName}_forCellTypist.mtx', adataOut.X)


with open(f'{tempOut}/{projectName}_gene.tsv', 'w') as f:
    for item in adataOut.var_names:
        f.write(item + '\n')
#zcat features.tsv.gz | cut -f 1 > ../gene.tsv


reference_adata = adata
cell_type_column = 'cell_type__custom'

# Train CellTypist model
model=celltypist.train(X=f'{tempOut}/{projectName}_forCellTypist.mtx',
                 genes=f"{tempOut}/{projectName}_gene.tsv",
                 labels=reference_adata.obs[cell_type_column],
                 n_jobs = 10,
                 max_iter = 100)

model.write(f"{SaveFolder}/{projectName}.pkl")



##################################################################
##################################################################
##################################################################
nuData=scvi.data.read_h5ad("./Thymus/h5ad/Thymoma_All.B.h5ad")

ct_predictions = celltypist.annotate(nuData,
                                     model =f"{SaveFolder}/{projectName}.pkl",
                                     majority_voting = True)
nuData0 = ct_predictions.to_adata()

nuData0.obs.rename(columns={'predicted_labels': 'ThymomaLabel',
  'over_clustering':'ThymomaOver_clustering',
  'majority_voting':'ThymomaMajority_voting',
  'conf_score': 'ThymomaConfScore'}, inplace=True)

atlas_predictions = celltypist.annotate(nuData0,
                                        model = 'Developing_Human_Thymus.pkl',
                                        majority_voting = True)
nuData1 = atlas_predictions.to_adata()
nuData1.obs.rename(columns={'predicted_labels': 'AtlasLabel',
  'over_clustering':'AtlasOver_clustering',
  'majority_voting':'AtlasMajority_voting',
  'conf_score': 'AtlasConfScore'}, inplace=True)

IMHigh = celltypist.annotate(nuData1,
                             model = 'Immune_All_High.pkl',
                             majority_voting = True)
nuData2 = IMHigh.to_adata()
nuData2.obs.rename(columns={'predicted_labels': 'ImHLabel',
  'over_clustering':'ImHOver_clustering',
  'majority_voting':'ImHMajority_voting',
  'conf_score': 'ImHConfScore'}, inplace=True)
IMLow = celltypist.annotate(nuData2,
                            model = 'Immune_All_Low.pkl',
                            majority_voting = True)
nuData = IMLow.to_adata()
nuData.obs.rename(columns={'predicted_labels': 'ImLLabel',
  'over_clustering':'ImLOver_clustering',
  'majority_voting':'ImLMajority_voting',
  'conf_score': 'ImLConfScore'}, inplace=True)

nuData.write(f"{SaveFolder}/{projectName}_Bcells.h5ad")

nuData.obs.to_csv(f"{SaveFolder}/{projectName}_Bcells_predictFromCellTypise.csv")

##################################################################
#################   convert to seurat object   ###################
##################################################################
###to seurate object
import scanpy as sc
from scipy import io
import sys
import os
import subprocess
import pandas as pd


out_dir = './Thymus/h5ad/B'
outMetaPath=f'{out_dir}/outMeta'
outDataPath=f'{out_dir}/outData'

if not os.path.exists(out_dir):
  os.makedirs(out_dir)
if not os.path.exists(outMetaPath):
  os.makedirs(outMetaPath)
if not os.path.exists(outDataPath):
  os.makedirs(outDataPath)


nuDataOut = nuData.copy()
with open(outDataPath + '/barcodes.tsv', 'w') as f:
  for item in nuDataOut.obs_names:
  f.write(item + '\n')

with open(outDataPath + '/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in nuDataOut.var_names]:
  f.write(item + '\n')

io.mmwrite(outDataPath +'/matrix', nuDataOut.X.T)

CmdGzip=["gzip %s/*" %outDataPath]
subprocess.call(CmdGzip,shell=True)


nuDataOut.obs.to_csv(outMetaPath + '/metadata.csv')
pd.DataFrame(adataOut.obsm["X_umap"], index=nuDataOut.obs_names).to_csv(outMetaPath + '/X_umap.csv')
pd.DataFrame(adataOut.obsm["X_pca"], index=nuDataOut.obs_names).to_csv(outMetaPath + '/X_pca.csv')
