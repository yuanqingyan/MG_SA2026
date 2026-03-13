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

# sc.settings.vector_friendly = True

projectName="ThymomaEP"

SaveFolder="./Thymus/h5ad"
tempOut="./Thymus/h5ad"

#### build the model
adata=scvi.data.read_h5ad("./Thymus/h5ad/ThymomaEP.h5ad")
adataOut = adata
io.mmwrite(f'{tempOut}/{projectName}_forCellTypist.mtx', adataOut.X) 


with open(f'{tempOut}/{projectName}_gene.tsv', 'w') as f:
  for item in adataOut.var_names:
  f.write(item + '\n')
#zcat features.tsv.gz | cut -f 1 > ../gene.tsv

reference_adata = adata
cell_type_column = 'cell_type__custom'


model=celltypist.train(X=f'{tempOut}/{projectName}_forCellTypist.mtx',
                       genes=f"{tempOut}/{projectName}_gene.tsv",
                       labels=reference_adata.obs[cell_type_column],
                       n_jobs = 10,
                       max_iter = 100)

model.write(f"{tempOut}/{projectName}.pkl")

nuData=scvi.data.read_h5ad("./Thymus/h5ad/NUThymusEP.h5ad")
nuData

sc.pp.pca(nuData, n_comps=50)
ct_predictions = celltypist.annotate(nuData,
                                     model =f"{tempOut}/{projectName}.pkl",
                                     majority_voting = True)
nuData0 = ct_predictions.to_adata()

nuData0.obs.rename(columns={'predicted_labels': 'ThymomaLabel',
  'over_clustering':'ThymomaOver_clustering',
  'majority_voting':'ThymomaMajority_voting',
  'conf_score': 'ThymomaConfScore'}, inplace=True)
nuData0
atlas_predictions = celltypist.annotate(nuData0,
                                        model = 'Developing_Human_Thymus.pkl',
                                        majority_voting = True)
nuData1 = atlas_predictions.to_adata()
nuData1.obs.rename(columns={'predicted_labels': 'AtlasLabel',
  'over_clustering':'AtlasOver_clustering',
  'majority_voting':'AtlasMajority_voting',
  'conf_score': 'AtlasConfScore'}, inplace=True)

nuData1.write(f"{SaveFolder}/{projectName}_EP_predictFromCellTypist.h5ad")
nuData1.obs.to_csv(f"{SaveFolder}/{projectName}_EP_predictFromCellTypist.csv")
