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
import math


print("Last run with scvi-tools version:", scvi.__version__)

projectName="Thymus_Manu"

DataFolder="./Projects/ManuscriptData/IntegrateThymus/data/toAnndata"
SaveFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData"
tempOut="./SC/Thymus/output/scVI/output"

adata=scvi.data.read_h5ad(f"{tempOut}/{projectName}_population.h5ad")


import scanpy as sc
from scipy import io
import sys
import os
import subprocess
import pandas as pd

out_dir = f'{tempOut}/Seurat'
outMetaPath=f'{out_dir}/outMeta'
outDataPath=f'{out_dir}/outData'

if not os.path.exists(out_dir):
  os.makedirs(out_dir)
if not os.path.exists(outMetaPath):
  os.makedirs(outMetaPath)   
if not os.path.exists(outDataPath):
  os.makedirs(outDataPath)   


adataOut = adata.raw.to_adata() #

with open(outDataPath + '/barcodes.tsv', 'w') as f:
  for item in adataOut.obs_names:
  f.write(item + '\n')

with open(outDataPath + '/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in adataOut.var_names]:
  f.write(item + '\n')

io.mmwrite(outDataPath +'/matrix', adataOut.X.T)

CmdGzip=["gzip %s/*" %outDataPath]
subprocess.call(CmdGzip,shell=True)


adataOut.obs.to_csv(outMetaPath + '/metadata.csv')
pd.DataFrame(adataOut.obsm["X_scVI"], index=adataOut.obs_names).to_csv(outMetaPath + '/X_scVI.csv')
pd.DataFrame(adataOut.obsm["X_umap"], index=adataOut.obs_names).to_csv(outMetaPath + '/X_umap.csv')
pd.DataFrame(adataOut.obsm["X_pca"], index=adataOut.obs_names).to_csv(outMetaPath + '/X_pca.csv')


