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


os.makedirs("./output/B2C/stardist", exist_ok=True)
outPath="./output/B2C/stardist"

B2CFolderPath='./Apr2025_H5ADRDS/B2C_H5AD'
pathRaw='./VisiumHD-set2'
SaveFolderPath="./CellTypist_H5AD"
FigOut="./output/QC_Fig"

modelPath="./Projects/ManuscriptData/IntegrateThymus/data/H5AD"

os.chdir(FigOut)

adataFile='./Projects/ManuscriptData/IntegrateThymus/publicVisium/25052546/res.cxg.h5ad'
adata=sc.read_h5ad(f"{adataFile}")

color_palette=["#377EB8", "#4DAF4A", "#984EA3", "#FF4500", "#17BECF", "#F781BF"]
len(set(color_palette))

for item in color_palette:
  if color_palette.count(item) > 1:
  print(item)

import itertools
pltCT=['cortex', 'junction', 'medulla', 'medulla_FN1', 'medulla_GC', 'stroma']

color_mapping = dict(zip(pltCT, color_palette))

adata.obs['niche_annot'].cat.reorder_categories(list(color_mapping.keys()))

adata.uns['niche_annot_colors'] = list(color_mapping.values())

##### Fig. 5E ###

plt.rcParams['figure.dpi']=600
plt.rcParams["figure.figsize"] = (10, 8)

sc.pl.umap(adata,
           color="niche_annot",
           frameon=False,
           add_outline=False,
           show=False,
           size=5,
           legend_loc="on data",
           legend_fontweight='light')
plt.savefig(f"{FigOut}/Public__umap_Niche.pdf",bbox_inches='tight')


sc.pl.umap(adata,
           color="niche_annot",
           frameon=False,
           add_outline=False,
           show=False,
           size=5,
           legend_fontsize='small',
           legend_fontweight='light')
plt.savefig(f"{FigOut}/Public__umap_Niche.tiff",bbox_inches='tight')


## Fig. 5E
sc.pl.umap(adata, 
           color=['MS4A1'], 
           wspace=0.4,
           vmax=3,
           vmin=0.5,
           size=5,
           show=False,
           frameon=False,
           cmap=reds)
plt.savefig(f"{FigOut}/Public__FP_MS4A1.pdf",bbox_inches='tight')

sc.pl.umap(adata, 
           color=['MS4A1'], 
           wspace=0.4,
           vmax=3,
           vmin=0.5,
           size=5,
           show=False,
           frameon=False,
           cmap=reds)
plt.savefig(f"{FigOut}/Public__FP_MS4A1.tiff",bbox_inches='tight')

sc.pl.umap(adata, 
           color=['CD83'], 
           wspace=0.4,
           vmax=3,
           size=5,
           vmin=0.5,
           show=False,
           frameon=False,
           cmap=reds)
plt.savefig(f"{FigOut}/Public__FP_CD83.pdf",bbox_inches='tight')

sc.pl.umap(adata, 
           color=['CD83'], 
           wspace=0.4,
           vmax=3,
           size=5,
           vmin=0.5,
           show=False,
           frameon=False,
           cmap=reds)
plt.savefig(f"{FigOut}/Public__FP_CD83.tiff",bbox_inches='tight')

sc.pl.umap(adata, 
           color=['CD70'], 
           wspace=0.4,
           vmax=1.2,
           size=5,
           vmin=0.3,
           show=False,
           frameon=False,
           cmap=reds)
plt.savefig(f"{FigOut}/Public__FP_CD70.pdf",bbox_inches='tight')


sc.pl.umap(adata, 
           color=['CD70'], 
           wspace=0.4,
           vmax=1.2,
           size=5,
           vmin=0.3,
           show=False,
           frameon=False,
           cmap=reds)
plt.savefig(f"{FigOut}/Public__FP_CD70.tiff",bbox_inches='tight')


#### Fig. 5F

for i in range(0,len(AllLib)):
	LabName=AllLib[i]
	adataSet=adata[adata.obs.sample_id==LabName]
	
	items_not_in2 = set(adataSet.obs["niche_annot"])-set(pltCT)
	if(len(items_not_in2)>0):
		print(f"Wrong in {LabName}")    
	
	
	adataCT=adataSet.obs["niche_annot"].unique()
	cm_in = {k: color_mapping[k] for k in color_mapping if k in adataCT}
	ctCat=list(cm_in.keys())
	ColCat=list(cm_in.values())
	
	adataSet.obs["niche_annot"] = pd.Categorical(
		adataSet.obs["niche_annot"], 
		categories=ctCat, 
		ordered=True  # Optional: Set to True if the order matters
		)
	sc.pl.spatial(adataSet, 
                  library_id=LabName,
                  img_key="lowres",
                  color=['niche_annot'], 
                  basis='spatial',
                  palette=ColCat,
                  size=1,
                  save=f'Public_Nich_{LabName}.pdf')
	
	
	
	
sc.pl.spatial(adata[adata.obs.sample_id=='TMG_1_1'], 
              library_id="TMG_1_1",
              img_key="lowres",
              color=['CD70'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=0.8,
              vmin=0.1,
              frameon=False,
              size=1,
             save=f'Public_CD70_TMG1_1_GC.pdf')
    
sc.pl.spatial(adata[adata.obs.sample_id=='TMG_1_1'], 
              library_id="TMG_1_1",
              img_key="lowres",
              color=['MS4A1'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=2,
              vmin=0.1,
              frameon=False,
              size=1,
              save=f'Public_MS4A1_TMG1_1_GC.pdf')


sc.pl.spatial(adata[adata.obs.sample_id=='TP_1_1'], 
              library_id="TP_1_1",
              img_key="lowres",
              color=['CD70'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=0.8,
              vmin=0.1,
              frameon=False,
              size=1,
              save=f'Public_CD70_TP_1_1_NoGC.pdf')



sc.pl.spatial(adata[adata.obs.sample_id=='TP_1_1'], 
              library_id="TP_1_1",
              img_key="lowres",
              color=['MS4A1'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=2,
              vmin=0.1,
              frameon=False,
              size=1,
              save=f'Public_MS4A1_TP_1_1_NoGC.pdf')



#### Fig. 6F
sc.pl.umap(adata[adata.obs.MG_status.isin(["nonMG"])], 
           color=['ADAM28'], 
           wspace=0.4,
           vmax=3,
           vmin=0.5,
           size=5,
           show=False,
           frameon=False,
           cmap=reds)
plt.savefig(f"{FigOut}/Public__FP_ADAM28_noMG.pdf",bbox_inches='tight')

sc.pl.umap(adata[adata.obs.MG_status.isin(["nonMG"])==False], 
           color=['ADAM28'], 
           wspace=0.4,
           vmax=3,
           vmin=0.5,
           size=5,
           show=False,
           frameon=False,
           cmap=reds)
plt.savefig(f"{FigOut}/Public__FP_ADAM28_MG.pdf",bbox_inches='tight')


TMG1=adata[adata.obs.sample_id=='TMG_1_1']
TP1=adata[adata.obs.sample_id=='TP_1_1']

sc.pl.spatial(TMG1, 
              library_id="TMG_1_1",
              img_key="lowres",
              color=['ADAM28'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=2,
              vmin=0.1,
              frameon=False,
              size=1,
              show=False)
plt.savefig(f"{FigOut}/TMG1_ADAM28.tiff",bbox_inches='tight')

sc.pl.spatial(TP1, 
              library_id="TP_1_1",
              img_key="lowres",
              color=['ADAM28'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=2,
              vmin=0.1,
              frameon=False,
              size=1,
              show=False)
plt.savefig(f"{FigOut}/TP1_ADAM28.tiff",bbox_inches='tight')


### Fig. 6D

sc.pl.spatial(TMG1, 
              library_id="TMG_1_1",
              img_key="lowres",
              color=['TNFRSF17'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=1,
              vmin=0.1,
              frameon=False,
              size=1,
              show=False)
plt.savefig(f"{FigOut}/TMG1_TNFRSF17.tiff",bbox_inches='tight')


sc.pl.spatial(TP1, 
              library_id="TP_1_1",
              img_key="lowres",
              color=['TNFRSF17'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=1,
              vmin=0.1,
              frameon=False,
              size=1,
              show=False)
plt.savefig(f"{FigOut}/TP1_TNFRSF17.tiff",bbox_inches='tight')

 ## Fig. 6B
sc.pl.spatial(TMG1, 
              library_id="TMG_1_1",
              img_key="lowres",
              color=['CD27'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=1.5,
              vmin=0.1,
              frameon=False,
              size=1,
              show=False)
plt.savefig(f"{FigOut}/TMG1_CD27.tiff",bbox_inches='tight')


sc.pl.spatial(TP1, 
              library_id="TP_1_1",
              img_key="lowres",
              color=['CD27'], 
              basis='spatial',
              palette="tab20",
              #cmap=reds,
              vmax=1.5,
              vmin=0.1,
              frameon=False,
              size=1,
              show=False)
plt.savefig(f"{FigOut}/TP1_CD27.tiff",bbox_inches='tight')


