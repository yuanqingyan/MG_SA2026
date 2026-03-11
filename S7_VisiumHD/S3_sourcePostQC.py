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
from copy import copy
reds = copy(plt.cm.Reds)


help_menu='''\nUsage: python sourcePost1_QCCelltypist.py -ID xxxx
**parameters**
'''

args=sys.argv
if '-h' in args or '-help' in args or len(args)==1:
  print (help_menu)
sys.exit(0)

if '-ID' not in args:
  sys.exit('SampleID needed')
else:
  i=sys.argv.index('-ID')
sampleName=sys.argv[i+1]



os.makedirs("./output/B2C/stardist", exist_ok=True)
outPath="./output/B2C/stardist"

B2CFolderPath='./Apr2025_H5ADRDS/B2C_H5AD'
pathRaw='./VisiumHD-set2'
SaveFolderPath="./CellTypist_H5AD"
FigOut="./output/QC_Fig"

modelPath="./Projects/ManuscriptData/IntegrateThymus/data/H5AD"

os.chdir(FigOut)



pCount="./VisiumHD/Apr2025/pThymus_All.txt"
sourcePath="./VisiumHD/Apr2025/Bin2Cell/sourceFile"

tmp='./Spatial/SlurmF'



###################################
##############  read ddata  #######
###################################
adata=sc.read_h5ad(f'{B2CFolderPath}/{sampleName}.b2c.h5ad')
adata

###################################
######### rename cell ID    #######
###################################
newCellID=[f"{sampleName}_" + str(i) for i in adata.obs.index]
adata.obs.index=newCellID

###################################
#########   plot hist       #######
###################################
plt.hist(adata.obs['bin_count'], bins=200, range=(0, 200),color='skyblue')
plt.savefig(f"{FigOut}/Fig1_{sampleName}_QC_hist.pdf")
plt.show()

plt.boxplot(adata.obs['bin_count'], vert=True)
plt.savefig(f"{FigOut}/Fig1_{sampleName}_QC_boxplot.pdf")
plt.show()

###################################
#########   filtering 1      ######
###################################
fdata = adata[adata.obs['bin_count']>=5] #
fdata = fdata[fdata.obs['bin_count']<=500] #
fdata.X.data = np.round(fdata.X.data)
fdata.raw = fdata.copy()
fdata

###################################
#########  add % mitoch      ######
###################################
fdata.var["mt"] = fdata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(fdata, qc_vars=["mt"], inplace=True)


###################################
#########  some plot         ######
###################################
fig, axs = plt.subplots(2, 2, figsize=(8, 8))
sns.histplot(fdata.obs["total_counts"], kde=False, ax=axs[0,0])
sns.histplot(
  fdata.obs["total_counts"][fdata.obs["total_counts"] < 500],
  kde=False,
  bins=40,
  ax=axs[0,1],
)

sns.histplot(fdata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1,0])
sns.histplot(
  fdata.obs["n_genes_by_counts"][fdata.obs["n_genes_by_counts"] < 200],
  kde=False,
  bins=60,
  ax=axs[1,1],
)
plt.savefig(f"{FigOut}/Fig2_{sampleName}_QC_hist.pdf")



###################################
######### filtering           #####
###################################
sc.pp.filter_cells(fdata, min_counts=80)#
sc.pp.filter_cells(fdata, max_counts=40000)
fdata = fdata[fdata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {fdata.n_obs}")
sc.pp.filter_genes(fdata, min_cells=3)
sns.histplot(fdata.obs["log1p_n_genes_by_counts"], kde=False)

###################################
######### normalization       #####
###################################
sc.pp.normalize_total(fdata, target_sum=1e4,inplace=True)##target_sum=1e4
sc.pp.log1p(fdata)
sc.pp.highly_variable_genes(fdata, flavor="seurat", n_top_genes=5000)

if sampleName not in ['17-Thymus']:
  b2c_predictions = celltypist.annotate(fdata, 
                                        model =f'{modelPath}/{sampleName}.pkl', 
                                        majority_voting = True)
fdata = b2c_predictions.to_adata()
fdata.obs.rename(columns={'predicted_labels': 'SampPredLabel', 
  'over_clustering':'SampOver_clustering',
  'majority_voting':'SampMajority_voting',
  'conf_score': 'SampConfScore'}, inplace=True)

b2c_pred2 = celltypist.annotate(fdata, 
                                model =f'{modelPath}/NuThy.pkl', 
                                majority_voting = True)
fdata = b2c_pred2.to_adata()
fdata.obs.rename(columns={'predicted_labels': 'PredLabel', 
  'over_clustering':'Over_clustering',
  'majority_voting':'Majority_voting',
  'conf_score': 'ConfScore'}, inplace=True)


# ###################################
# ######### final plotting        ###
# ###################################

sc.pl.spatial(fdata, 
              img_key="mpp0.5", 
              basis='spatial_cropped150_buffer', 
              color=["PredLabel"], 
              size=8,
              save=f'Fig5_{sampleName}_QC_Spat.celltypist.pdf')


# ###################################
# ######### save h5ad             ###
# ###################################
fdata.write_h5ad(f"{SaveFolderPath}/{sampleName}_celltypist.H5AD")



