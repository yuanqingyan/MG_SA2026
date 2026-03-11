#!/usr/bin/env python
# coding: utf-8

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
sys.path.insert(0, './code/sourcePython/CellTypeMarker')  #
import Thymus


print("Last run with scvi-tools version:", scvi.__version__)

# sc.settings.vector_friendly = True

projectName="Thymus_Manu"

DataFolder="./Projects/ManuscriptData/IntegrateThymus/data/toAnndata"
SaveFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData"
tempOut="/./SC/Thymus/output/scVI/output"

FigOut=f"{tempOut}/GeneMarker"
if not os.path.exists(FigOut):
  os.makedirs(FigOut)


adata=scvi.data.read_h5ad(f"{SaveFolder}/{projectName}_scVI.res1.h5ad")



Pop_dict = {'Epithelial': ['21','7','9','17','30'],
  'Endothelial':['20'],
  'Stromal':['19'],
  'Immune':['0','1','2','3','4','5','6','8','10','11','12','13','14','15','16','18','22','23','24','25','26','27','28','29','31']}
adata.obs['Population'] = np.nan

for i in Pop_dict.keys():
  ind = pd.Series(adata.obs.leiden_scVI).isin(Pop_dict[i])
adata.obs.loc[ind,'Population'] = i


Pop2_dict = {'Structual': ['21','7','9','17','30','20','19'],
  'IM':['0','1','2','3','4','5','6','8','10','11','12','13','14','15','16','18','22','23','24','25','26','27','28','29','31']}
adata.obs['Population2'] = np.nan


for i in Pop2_dict.keys():
  ind = pd.Series(adata.obs.leiden_scVI).isin(Pop2_dict[i])
adata.obs.loc[ind,'Population2'] = i

adata.write(f"{tempOut}/{projectName}_population.h5ad")


############################################################################################
Study_dict = {'NU': ['ThymusMay242023','Thymus_June2023_GEMonly','Thymus_Sep2023_GEMonly',
                     'Thymus_Dec2023_GEMonly','Thymus_Oct2023_GEMonly'],
  'Xin':['Xin'],
  'Yasumizu':['Yasumizu']}
adata.obs['Study'] = np.nan

for i in Study_dict.keys():
  ind = pd.Series(adata.obs.GSE).isin(Study_dict[i])
adata.obs.loc[ind,'Study'] = i
pd.crosstab(adata.obs.GSE,adata.obs.Study)

############################################################################################
adata.write(f"{tempOut}/{projectName}_popu_meta.h5ad") 


for ct in adata.obs.Population2.unique():
  temp=adata[adata.obs.Population2.isin([ct])]
temp.write(f"{tempOut}/Pop.{ct}.h5ad")


