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

sampleIn=['6-Thymus','17-Thymus', '12-thymoma']
allAda=[sc.read_h5ad(f"{SaveFolderPath}/{i}_celltypist.H5AD") for i in sampleIn]
len(allAda)

for i in range(len(sampleIn)):
  allAda[i].obs['SampleName']=sampleIn[i]
allAda[0].obs.head(3) 
cutoff=0.7
Dat=[allAda[0][allAda[0].obs.ConfScore>=cutoff],
     allAda[1][allAda[1].obs.ConfScore>=cutoff],
     allAda[2][allAda[2].obs.ConfScore>=cutoff]]

sns.boxplot(x='PredLabel', y='ConfScore', data=Dat[0].obs, palette='viridis')
plt.xticks(rotation=90, ha='center')
plt.tight_layout()
plt.show()



### niche analysis using cell type
sys.path.insert(0, './VisiumHD')
import S7_sourceIdentifyCellNich2

MergNich6=S7_sourceIdentifyCellNich2.IdentifyNich(adata=Dat[0], 
                                         group_by='PredLabel', 
                                         cluster_name="niches", 
                                         neighbors_k=60, 
                                         niches_k=8)
selected_genes3=["RAG1","DNTT","CCR9","CCL25", ###cortex
                 "CCR7", "AIRE", "CCL19",'CCL21', ###Medulla
                 "CD79A","MS4A1",'IGKC','MKI67',"CXCL13","TNFSF13B",'CD70', ##GC
                 "FN1","COL1A1"]

sc.pl.matrixplot(
  MergNich6,
  var_names=selected_genes3, # Genes to plot (e.g., top marker genes or specific genes)
  groupby='niches',
  dendrogram=False,
  cmap="Reds",
  standard_scale="var",
  colorbar_title="Mean expression",
)

Nich_dict6 = {'Cortex':[3,6],
  'Medulla_Stroma':[0,2,4,7],
  'Medulla_B':[1,5]}
MergNich6.obs['F_Niche'] = np.nan
for i in Nich_dict6.keys():
  ind = pd.Series(MergNich6.obs['niches']).isin(Nich_dict6[i])
MergNich6.obs.loc[ind,'F_Niche'] = i
MergNich6.obs['F_Niche'] = MergNich6.obs['F_Niche'].astype(str).astype('category')
pd.crosstab(MergNich6.obs['F_Niche'], MergNich6.obs['niches'])
MergNich6.write(f"{SaveFolderPath}/6_thymus_addNiche.h5ad")


MergNich17=S7_sourceIdentifyCellNich2.IdentifyNich(adata=Dat[1], 
                                          group_by='PredLabel', 
                                          cluster_name="niches", 
                                          neighbors_k=100, 
                                          niches_k=10)###

testForSub=MergNich17[MergNich17.obs.niches.isin([5])].copy()
subC5_17=S7_sourceIdentifyCellNich2.IdentifyNich(adata=testForSub, 
                                        group_by='PredLabel', 
                                        cluster_name="SubNich", 
                                        neighbors_k=30, 
                                        niches_k=6)

MeduBID=testForSub[testForSub.obs.SubNich.isin([0,1,2,3])].obs.index
len(MeduBID)


Nich_dict17 = {'Cortex':[0,1,2,3,4,6,7,8],
  'Medulla_Stroma':[5,9]}
MergNich17.obs['F_Niche'] = np.nan
for i in Nich_dict17.keys():
  ind = pd.Series(MergNich17.obs['niches']).isin(Nich_dict17[i])
MergNich17.obs.loc[ind,'F_Niche'] = i
MergNich17.obs['F_Niche'] = MergNich17.obs['F_Niche'].astype('string')
MergNich17.obs['F_Niche'] = MergNich17.obs['F_Niche'].astype(str).astype('category')
MergNich17.write(f"{SaveFolderPath}/17_thymus_addNiche.h5ad")

