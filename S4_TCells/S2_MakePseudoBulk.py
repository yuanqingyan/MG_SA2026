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


projectName="Thymus_Manu_bulkT"

DataFolder="./Projects/ManuscriptData/IntegrateThymus/data/toAnndata"
SaveFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData"
tempOut="./SC/Thymus/output/scVI/output"

adata=scvi.data.read_h5ad("./Thymus/FigOut/forPseudoBulkR.h5ad")
adata.obs['TSub_F'] = adata.obs['TSub_F'].replace('CD8a/b(entry)', 'CD8ab(entry)')
adata.obs['TSub_F'] = adata.obs['TSub_F'].replace('CD8a/a', 'CD8aa')
adata.obs['TSub_F'] = adata.obs['TSub_F'].replace('CD4.Tem(Th1/Th17)', 'CD4.Tem(Th1Th17)')

CT=adata.obs
adata_conc=scvi.data.read_h5ad(f"{SaveFolder}/adata_Thymus_Manu_concatenate.h5ad")
adata_conc.obs = adata_conc.obs.merge(CT, left_index=True, right_index=True, how='left')
rawCount=adata_conc[adata_conc.obs.index.isin(CT.index)]

meta=CT

savePath=f"{SaveFolder}/PseudoBulk/SubTcell"
CUTOFF = 25

for i in meta.TSub_F.unique():
    ct=i
    cells = meta.index[meta.TSub_F == ct] 
    patient_status = meta.loc[cells, ['ptDis2']].drop_duplicates() 
    patient_status = patient_status[['ptDis2']].drop_duplicates()
    print(patient_status.shape)
    samples = patient_status["ptDis2"].unique() 
    print(len(samples))
    sample_values = [] 
    sample_ncells = []  
    filtered_samples = [] #
    n_cells = []
    for s in samples:
        s_cells = rawCount.obs_names.isin(cells) & (rawCount.obs["ptDis2"] == s) 
        if s_cells.sum() >= CUTOFF:
            sample_values.append(rawCount.X[s_cells, :].sum(axis=0).A[0])
            sample_ncells.append((rawCount.X[s_cells, :] > 0).sum(axis=0).A[0])
            filtered_samples.append(s)
            n_cells.append(s_cells.sum())
            
    sample_values = pd.DataFrame(sample_values, index=filtered_samples, columns=rawCount.var_names).T
    sample_ncells = pd.DataFrame(sample_ncells, index=filtered_samples, columns=rawCount.var_names).T

    sample_values.to_csv("%s/sampleValues_%s.txt" %(savePath,ct), sep="\t")
    sample_ncells.to_csv("%s/sample_NCells_%s.txt" %(savePath,ct), sep="\t")
    
    patient_status = patient_status.loc[patient_status["ptDis2"].isin(filtered_samples), :]
    print(n_cells)
    patient_status["n_cells"] = n_cells
    patient_status.to_csv("%s/sampleMeta_%s.csv" %(savePath,ct))
    print(f"{ct} done, {len(filtered_samples)}/{len(samples)}")

