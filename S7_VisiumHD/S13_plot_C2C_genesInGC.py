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


adaWithGC=sc.read_h5ad(f"{SaveFolderPath}/Merg6_17_withGC.h5ad")
gcID=adaWithGC[adaWithGC.obs.Germinal=='GC'].obs.index.str.replace('_0$', '', regex=True)


for i in range(len(allAda)):
  allAda[i].obs['Germinal']=allAda[i].obs['FNich2'].astype('str')
allAda[i].obs.loc[allAda[i].obs.index.isin(gcID), 'Germinal'] = "GC"

pd.crosstab(allAda[0].obs.Germinal, allAda[0].obs.FNich2)


allGC=allAda

allGC[0].obs.Germinal.value_counts()

color_GC=["#377EB8","#000000","#E41A1C"]
pltGC=['Cortex','GC','Medulla_Stroma']

color_mapGC = dict(zip(pltGC, color_GC))

for i in range(len(allGC)):
  allGC[i].obs["Germinal"] = pd.Categorical(allGC[i].obs["Germinal"], categories=pltGC,ordered=True)

mask1 = ((allAda[0].obs['array_row'] >=1040) & 
           (allAda[0].obs['array_row'] <= 1340) & 
           (allAda[0].obs['array_col'] >= 1700) & 
           (allAda[0].obs['array_col'] <= 2000)
)

bdata1 = allAda[0][mask1] 


mask2 = ((allAda[0].obs['array_row'] >=785) & 
           (allAda[0].obs['array_row'] <= 1085) & 
           (allAda[0].obs['array_col'] >= 1340) & 
           (allAda[0].obs['array_col'] <= 1640)
)

bdata2 = allAda[0][mask2]




#### Fig. 6A

sc.pl.spatial(allAda[0], 
              color=["CD27"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=2,
              show=False,
              vmax=2,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Thy6_CD27.pdf",bbox_inches='tight')

sc.pl.spatial(allAda[1], 
              color=["CD27"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=2,
              show=False,
              vmax=2,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Thy17_CD27.pdf",bbox_inches='tight')


sc.pl.spatial(bdata1, 
              color=["CD27"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=2,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom1_CD27.pdf",bbox_inches='tight')

sc.pl.spatial(bdata2, 
              color=["CD27"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=2,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom2_CD27.pdf",bbox_inches='tight')



#### Fig. 6C
sc.pl.spatial(allAda[0], 
              color=["TNFRSF17"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=2,
              show=False,
              vmax=1,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Thy6_TNFRSF17.pdf",bbox_inches='tight')

sc.pl.spatial(allAda[1], 
              color=["TNFRSF17"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=2,
              show=False,
              vmax=1,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Thy17_TNFRSF17.pdf",bbox_inches='tight')


sc.pl.spatial(bdata1, 
              color=["TNFRSF17"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=1,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom1_TNFRSF17.pdf",bbox_inches='tight')

sc.pl.spatial(bdata2, 
              color=["TNFRSF17"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=1,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom2_TNFRSF17.pdf",bbox_inches='tight')

#### Fig. 6E
sc.pl.spatial(allAda[0], 
              color=["ADAM28"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=2,
              show=False,
              vmax=1.5,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Thy6_ADAM28.pdf",bbox_inches='tight')

sc.pl.spatial(allAda[1], 
              color=["ADAM28"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=2,
              show=False,
              vmax=1.5,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Thy17_ADAM28.pdf",bbox_inches='tight')


sc.pl.spatial(bdata1, 
              color=["ADAM28"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=1.5,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom1_ADAM28.pdf",bbox_inches='tight')

sc.pl.spatial(bdata2, 
              color=["ADAM28"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=1.5,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom2_ADAM28.pdf",bbox_inches='tight')



