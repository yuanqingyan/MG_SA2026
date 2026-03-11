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


allGC=allAda#[allAda[0][allAda[0].obs["FNich2"]=="Medulla_Stroma"],allAda[1][allAda[1].obs["FNich2"]=="Medulla_Stroma"]]

allGC[0].obs.Germinal.value_counts()

color_GC=["#377EB8","#000000","#E41A1C"]
pltGC=['Cortex','GC','Medulla_Stroma']

color_mapGC = dict(zip(pltGC, color_GC))

for i in range(len(allGC)):
  allGC[i].obs["Germinal"] = pd.Categorical(allGC[i].obs["Germinal"], categories=pltGC,ordered=True)


rename_dict = {
  "DPT":"DPT",
  "CD4.Naive":"CD4T.Naive",
  "mcTEC":"Ep",
  "MemB.Rib+Mit+":"B",
  "Plasma":"Plasma",
  "Macrophage":"Mye",
  "Prolf.T":"Prolif.T",
  "KRT15+.mTEC":"Ep",
  "Fibroblast":"St",
  "NaivB":"B",
  "VSMC":"St",
  "CD8.Tem":"CD8T",
  "CD8a/b(entry)":"CD8T",
  "CD8T.Naive":"CD8T",
  "MemB.Rib+":"B",
  "T(agonist)":"T(agonist)",
  "SLPI+.mTEC":"Ep",
  "DNT":"DNT",
  "Activated.En":"En",
  "Vascular":"En",
  "Endo":"En",
  "CD8a/a":"CD8T",
  "CD4.Tcm(Th2)":"CD4T",
  "DC":"Mye",
  "Pericyte":"St",
  "Treg":"CD4T",
  "NK":"NK/ILC",
  "MemB":"B",
  "CD4.Tcm(Th0)":"CD4T",
  "CD4.Tcm(Th17)":"CD4T",
  "cTEC":"Ep",
  "CD8.Temra":"CD8T",
  "CD8.Tcm":"CD8T",
  "CD8.Trm":"CD8T",
  "CRTAM+.gdT":"gdT",
  "gdT":"gdT",
  "MemB.Mit+":"B",
  "CD4.Tem(Th1)":"CD4T",
  "Tfh":"CD4T",
  "MemB.ISG":"B",
  "CD4.Tem(Th1/Th17)":"CD4T",
  "Neutrophils":"Mye",
  "ILC":"NK/ILC",
  "Lymphatic":"En",
  'Switch':'B'
}


for i in range(len(allGC)):
  allGC[i].obs["NewCT"] =[rename_dict[elem] if elem in rename_dict else elem for elem in allGC[i].obs["PredLabel"]]
allGC[i].obs['NewCT'] = allGC[i].obs['NewCT'].astype('category')
allGC[i].obsm['spatialOrig'] = allGC[i].obsm['spatial']
allGC[i].obsm['spatial'] = allGC[i].obsm['spatial_cropped150_buffer']

pd.crosstab(allGC[0].obs.PredLabel, allGC[0].obs.NewCT)

MG=allGC[0].copy()
NoMG=allGC[1].copy()

shared_intervals = np.linspace(0, 1000, 50)
sq.gr.spatial_neighbors(MG, 
                        coord_type="generic", 
                        delaunay=True)

sq.gr.centrality_scores(MG, 
                        cluster_key="NewCT")

sq.gr.co_occurrence(MG, 
                    cluster_key="NewCT",
                    interval=shared_intervals)

sq.gr.nhood_enrichment(MG,
                       cluster_key="NewCT")


##################################################
sq.gr.spatial_neighbors(NoMG, 
                        coord_type="generic", 
                        delaunay=True)

sq.gr.centrality_scores(NoMG, 
                        cluster_key="NewCT")

sq.gr.co_occurrence(NoMG, 
                    cluster_key="NewCT",
                    interval=shared_intervals)

sq.gr.nhood_enrichment(NoMG,
                       cluster_key="NewCT")

MG.write(f"{SaveFolderPath}/MG_visiumHD_Coo.h5ad")
NoMG.write(f"{SaveFolderPath}/NoMG_visiumHD_Coo.h5ad")

MG=sc.read_h5ad(f"{SaveFolderPath}/MG_visiumHD_Coo.h5ad")
NoMG=sc.read_h5ad(f"{SaveFolderPath}/NoMG_visiumHD_Coo.h5ad")

from matplotlib import rcParams
rcParams["pdf.fonttype"] = 42  # Ensures text is stored as text, not paths
rcParams["ps.fonttype"] = 42


colCT=["#D85085",  "#845172", '#FFD700', "#3983AC", "#999999", "#FFB315",
       "#E41A1C","#9B4F9D","#F07816",'#0000FF','#00FF00',"#00BFFF",
       "#9400D3", '#000080','#AAFFC3']
facCT=['B','CD4T','CD4T.Naive','CD8T','DNT','DPT',
       'En','Ep','NK/ILC','Mye','Plasma','Prolif.T',
       'St','T(agonist)','gdT']

CTcol=pd.DataFrame({"CT": facCT,'colSch':colCT})


MG.obs['NewCT']=pd.Categorical(MG.obs['NewCT'],
                               categories=CTcol.CT,
                               ordered=True)

MG.uns['NewCT_colors']=CTcol.colSch

NoMG.obs['NewCT']=pd.Categorical(NoMG.obs['NewCT'],
                                 categories=CTcol.CT,
                                 ordered=True)

NoMG.uns['NewCT_colors']=CTcol.colSch


vmin=-5
vmax=5

sq.pl.nhood_enrichment(
  MG, 
  cluster_key="NewCT", 
  vmin=vmin, 
  vmax=vmax, 
  cmap="vlag",
  title="MG",
  show=False
)

#### Fig. S10A
pdf_path = f"{FigOut}/MG_Coo_plot_V.pdf"
plt.savefig(pdf_path, format="pdf", bbox_inches="tight")


sq.pl.nhood_enrichment(
  NoMG, 
  cluster_key="NewCT", 
  vmin=vmin, 
  vmax=vmax, 
  cmap="vlag",
  title="NoMG",
  show=False
)

pdf_path = f"{FigOut}/NoMG_Coo_plot_v.pdf"
plt.savefig(pdf_path, format="pdf", bbox_inches="tight")


