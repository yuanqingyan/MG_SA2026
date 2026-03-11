#!/usr/bin/env python
import scvelo as scv
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import scanpy as sc
import cellrank as cr
import sys
import matplotlib.pyplot as plt

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2


############################
dataPath='./Projects/ManuscriptData/IntegrateThymus/data/Velocyte'
FigOut="./Thymus/FigOut"
h5adName='MG_Bcell'
############################



pCount_Dec="./SC/Thymus/Dec2023_nucore/pCount.txt" 
BharatCountFolder_Dec="./Projects/SC/Thymus/Diego/Dec2023_nucore/Count"

pCount_June="./SC/Thymus/June2023_nucore/pCount.txt" 
BharatCountFolder_June="./Projects/SC/Thymus/Diego/June2023_nucore/Count"

pCount_May="./SC/Thymus/May242023/pCount.txt" ###Sample ID
BharatCountFolder_May="./Projects/SC/Thymus/Diego/May242023/deepSeq" #count data Folder

pCount_Oct="./SC/Thymus/Oct2023_nucore/pCount.txt" 
BharatCountFolder_Oct="./Projects/SC/Thymus/Diego/Oct2023_nucore/Count"

pCount_Sep="./SC/Thymus/Sep2023_nucore/pCount.txt" 
BharatCountFolder_Sep="./Projects/SC/Thymus/Diego/Sep2023_nucore/Count"


sampleName=[]
loomPath=[]
infile_May=open("%s" %pCount_May,"r").readlines()
infile_June=open("%s" %pCount_June,"r").readlines()
infile_Sep=open("%s" %pCount_Sep,"r").readlines()
infile_Oct=open("%s" %pCount_Oct,"r").readlines()
infile_Dec=open("%s" %pCount_Dec,"r").readlines()

for line in infile_May:
  name=line.strip().strip('"').split("\t")[0]
sampleName.append(name)
loomPath.append(BharatCountFolder_May)

for line in infile_June:
  name=line.strip().strip('"').split("\t")[0]
sampleName.append(name)  
loomPath.append(BharatCountFolder_June)

for line in infile_Sep:
  name=line.strip().strip('"').split("\t")[0]
sampleName.append(name)  
loomPath.append(BharatCountFolder_Sep)

for line in infile_Oct:
  name=line.strip().strip('"').split("\t")[0]
sampleName.append(name)  
loomPath.append(BharatCountFolder_Oct)

for line in infile_Dec:
  name=line.strip().strip('"').split("\t")[0]
sampleName.append(name)
loomPath.append(BharatCountFolder_Dec)
print (sampleName)

loomData=[scv.read('%s/%s/velocyto/%s.loom' % (loomPath[isample],sampleName[isample],sampleName[isample]), cache=True) for isample in range(0,len(sampleName))]


for inumber in range(len(sampleName)):
  loomData[inumber].obs['origidentid']= [bc.split(':')[0] for bc in loomData[inumber].obs.index.tolist()]
tempSam=sampleName[inumber]
print(gse_dict)

tempBar=[bc.split(':')[1] for bc in loomData[inumber].obs.index.tolist()]
tempBarcodes = [bc[0:len(bc)-1] + '-1' for bc in tempBar]
loomData[inumber].obs['TempBarcodes']=tempBarcodes
conCol=loomData[inumber].obs['TempBarcodes'].astype(str)+"_"+loomData[inumber].obs['indind'].astype(str)+"_"+loomData[inumber].obs['GSEName'].astype(str)
loomData[inumber].obs.index = conCol

loomData[inumber].var_names_make_unique()

loomDataCon = loomData[0].concatenate(loomData[1:])
loomDataCon.write(f'{dataPath}/AllThymusData.h5ad')


adata = sc.read_h5ad(f'{FigOut}/forRNAVelocity.h5ad')
adata = scv.utils.merge(adata, loomDataCon)

adata.write(f'{dataPath}/B_merge_loom.h5ad')


###########################################################################
###########################################################################
###########################################################################

adata3=sc.read_h5ad(f'{dataPath}/B_merge_loom.h5ad')
adata3

sc.pp.neighbors(adata3, n_neighbors=80, use_rep='X_pca')
scv.pp.filter_and_normalize(adata3,min_shared_counts=20,n_top_genes=5000) 
scv.pp.moments(adata3,n_pcs=30,n_neighbors=80)

scv.tl.recover_dynamics(adata3)
scv.tl.velocity(adata3, mode='dynamical')
scv.tl.velocity_graph(adata3)
adata3.write('%s/%s_dynamical.h5ad' %(dataPath,h5adName))


###########################################################################
###########################################################################
###########################################################################

adata3=sc.read_h5ad('%s/%s_dynamical.h5ad' %(dataPath,h5adName))
plt.rcParams["figure.figsize"] = (10, 10)
scv.pl.velocity_embedding_stream(adata3, 
                                 basis='umap',
                                 color='Bsub', 
                                 title='')# 


BcellLevel=["Centroblast_GC","Centrocyte_GC","Early_GC","Mature_GC",
            "NaivB", "Naiv_Cir_1", "Naiv_Cir_2", "Naiv_Cir_3", "Naiv_Cir_4", "Naiv_Cir_5",
            "Un_switch","Switch",
            "MemB_Cir","MemB",
            "MemB.ISG",
            "MemB.Rib+Mit+","MemB.Rib+","MemB.Mit+","MemB.EBI3+",
            "Plasma"]
celltypeColor=["#8A9D12", "#00FF00", "#0000FF", "#FF0000", "#11FFFF", "#FF00FF", "#FFA500", "#800080", "#FFC0CB", 
               "#A52A2A", "#B0F111", "#308080", "#FFD700", "#4682B4", "#DC143C", "#00CED1", "#ADFF2F", "#FF1493", 
               "#2E8B57", "#7B68EE"]
adata3.obs['Bsub'].cat.reorder_categories(BcellLevel, inplace=True)
adata3.uns['Bsub_colors'] = celltypeColor

############### Fig. S5 ###########
scv.pl.velocity_embedding_grid(adata3, 
                               basis='umap', 
                               color='Bsub', 
                               save='BcellSub.embedding_grid.pdf', 
                               title='', 
                               scale=0.2)#


