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
adaWithGC.obs.head(3)

gcID=adaWithGC[adaWithGC.obs.Germinal=='GC'].obs.index.str.replace('_0$', '', regex=True)

for i in range(len(allAda)):
  allAda[i].obs['Germinal']=allAda[i].obs['FNich2'].astype('str')
allAda[i].obs.loc[allAda[i].obs.index.isin(gcID), 'Germinal'] = "GC"

### matrixplot to show the identification of GC
selected_genesF=["CD79A","MS4A1","CXCL13",'CXCR4',
                 'BCL6','ICOS','CXCR5','CD83','PRDM1',##GC
                 "FN1","COL1A1",'FAS','MKI67']


sc.pl.matrixplot(
  adaWithGC,
  var_names=selected_genesF, 
  groupby='cluster_cellcharter',
  dendrogram=False,
  cmap="Reds",
  standard_scale="var",
  colorbar_title="Mean expression",
  show=False
)


#### matrixplot to show the identification of GC
## Fig. 5C
selected_genesF=["CD79A","MS4A1","CXCL13",'CXCR4','CXCR5','MKI67']


sc.pl.matrixplot(
  adaWithGC,
  var_names=selected_genesF, # Genes to plot (e.g., top marker genes or specific genes)
  groupby='Germinal',
  dendrogram=False,
  cmap="Reds",
  standard_scale="var",
  colorbar_title="Mean expression",
  show=False
)

plt.rcParams['pdf.fonttype'] = 42
plt.savefig(f"{FigOut}/matrixplt_GC.pdf",bbox_inches='tight')


mergDat=ad.concat(allAda,join='outer',axis=0)
mergDat.obs['newCN']=mergDat.obs['SampleName'].astype('str')+mergDat.obs['FNich2'].astype('str')


selected_genesF=["RAG1","DNTT","CCR9","CCL25", ###cortex
                 "CCR7", "AIRE", "CCL19", ###Medulla
                 "CD79A","MS4A1",'IGKC',"CXCL13", ##GC
                 "FN1","COL1A1"]


sc.pl.matrixplot(
  mergDat,
  var_names=selected_genesF, 
  groupby='newCN',
  dendrogram=False,
  cmap="Reds",
  standard_scale="var",
  colorbar_title="Mean expression",
  show=False
)

### Fig. 5A
selected_genesF=["RAG1","DNTT","CCR9","CCL25", ###cortex
                 "CCR7", "AIRE", "CCL19", ###Medulla
                 "CD79A","MS4A1",'IGKC',"CXCL13", ##GC
                 "FN1","COL1A1"]


sc.pl.matrixplot(
  mergDat,
  var_names=selected_genesF, # Genes to plot (e.g., top marker genes or specific genes)
  groupby='FNich2',
  dendrogram=False,
  cmap="Reds",
  standard_scale="var",
  colorbar_title="Mean expression",
  show=False
)

plt.rcParams['pdf.fonttype'] = 42
plt.savefig(f"{FigOut}/matrixplt_Niches.pdf",bbox_inches='tight')


color_palette=["#377EB8","#E41A1C"]
pltN=['Cortex','Medulla_Stroma']

color_palette17=["#377EB8","#E41A1C"]
pltN17=['Cortex','Medulla_Stroma']

color_mapping = dict(zip(pltN, color_palette))
color_mapping17 = dict(zip(pltN17, color_palette17))

allAda[0].obs["FNich2"] = pd.Categorical(allAda[0].obs["FNich2"], categories=pltN,ordered=True)
allAda[1].obs["FNich2"] = pd.Categorical(allAda[1].obs["FNich2"], categories=pltN17,ordered=True)

## Fig. 5B
sc.pl.spatial(allAda[1], 
              img_key="mpp0.5", 
              basis='spatial_cropped150_buffer', 
              color=["FNich2"], 
              size=2,
              show=False,
              frameon=False,
              palette=color_mapping17)
plt.savefig(f"{FigOut}/Niche_Thy17.tiff",bbox_inches='tight')


## Fig. 5B
sc.pl.spatial(allAda[0], 
              img_key="mpp0.5", 
              basis='spatial_cropped150_buffer', 
              color=["FNich2"], 
              size=2,
              show=False,
              frameon=False,
              palette=color_mapping)
plt.savefig(f"{FigOut}/Niche_Thy6.tiff",bbox_inches='tight')


##### Fig. 5B
sc.pl.spatial(allAda[0], 
              img_key="mpp0.5", 
              basis='spatial_cropped150_buffer', 
              color=["MS4A1"], 
              size=2,
              show=False,
              vmax=3,
              vmin=0.5,
              frameon=False,#cmap=reds
)
plt.savefig(f"{FigOut}/MS4A1_Thy6.pdf",bbox_inches='tight')

sc.pl.spatial(allAda[1], 
              img_key="mpp0.5", 
              basis='spatial_cropped150_buffer', 
              color=["MS4A1"], 
              size=2,
              show=False,
              vmax=3,
              vmin=0.5,
              frameon=False,#cmap=reds
)
plt.savefig(f"{FigOut}/MS4A1_Thy17.pdf",bbox_inches='tight')



sc.pl.spatial(allAda[0], 
              img_key="mpp0.5", 
              basis='spatial_cropped150_buffer', 
              color=["CCL25"], 
              size=2,
              show=False,
              vmax=3,
              vmin=0.5,
              frameon=False,#cmap=reds
)
plt.savefig(f"{FigOut}/CCL25_Thy6.pdf",bbox_inches='tight')

sc.pl.spatial(allAda[1], 
              img_key="mpp0.5", 
              basis='spatial_cropped150_buffer', 
              color=["CCL25"], 
              size=2,
              show=False,
              vmax=3,
              vmin=0.5,
              frameon=False,#cmap=reds
)
plt.savefig(f"{FigOut}/CCL25_Thy17.pdf",bbox_inches='tight')



allAda[0].write(f'{FigOut}/Thy6_withNich.h5ad')

allGC=allAda

color_GC=["#377EB8","#000000","#E41A1C"]
pltGC=['Cortex','GC','Medulla_Stroma']

color_mapGC = dict(zip(pltGC, color_GC))

for i in range(len(allGC)):
  allGC[i].obs["Germinal"] = pd.Categorical(allGC[i].obs["Germinal"], categories=pltGC,ordered=True)

######## Fig 5C#########
sc.pl.spatial(allGC[i], 
              img_key="mpp0.5", 
              basis='spatial_cropped150_buffer', 
              color=["Germinal"], 
              size=2,
              show=False,
              frameon=False,
              palette=color_mapGC)
plt.savefig(f"{FigOut}/GC_samp_{i}.pdf",bbox_inches='tight')



mask1 = ((allAda[0].obs['array_row'] >=1040) & 
           (allAda[0].obs['array_row'] <= 1340) & 
           (allAda[0].obs['array_col'] >= 1700) & 
           (allAda[0].obs['array_col'] <= 2000)
)

bdata1 = allAda[0][mask1] # read a crop of the object 

sc.set_figure_params(figsize=[10,10],dpi=100)
sc.pl.spatial(bdata1, 
              color=["Germinal"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              frameon=False,
              palette=color_mapGC)

plt.savefig(f"{FigOut}/Zoom1_Thy6.pdf",bbox_inches='tight')

sc.pl.spatial(bdata1, 
              color=["CXCL13"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=0.8,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom1_CXCL13.pdf",bbox_inches='tight')


sc.pl.spatial(bdata1, 
              color=["CD70"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=0.8,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom1_CD70.pdf",bbox_inches='tight')




mask2 = ((allAda[0].obs['array_row'] >=785) & 
           (allAda[0].obs['array_row'] <= 1085) & 
           (allAda[0].obs['array_col'] >= 1340) & 
           (allAda[0].obs['array_col'] <= 1640)
)

bdata2 = allAda[0][mask2] # read a crop of the object 

sc.set_figure_params(figsize=[10,10],dpi=100)
sc.pl.spatial(bdata2, 
              color=["Germinal"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              frameon=False,
              palette=color_mapGC)
plt.savefig(f"{FigOut}/Zoom2_Thy6.pdf",bbox_inches='tight')

sc.pl.spatial(bdata2, 
              color=["CD70"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=0.8,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom2_CD70.pdf",bbox_inches='tight')


sc.pl.spatial(bdata2, 
              color=["CXCL13"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=0.8,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom2_CXCL13.pdf",bbox_inches='tight')



sc.pl.spatial(bdata1, 
              color=["CD79A"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=1,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom1_CD79A.pdf",bbox_inches='tight')

sc.pl.spatial(bdata2, 
              color=["CD79A"],
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              size=1,
              show=False,
              vmax=1,
              vmin=0.1,
              frameon=False)
plt.savefig(f"{FigOut}/Zoom2_CD79A.pdf",bbox_inches='tight')



