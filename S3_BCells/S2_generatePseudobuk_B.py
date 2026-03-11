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


adata2=sc.read_h5ad(f"{SaveFolder}/{projectName}_AddCT.h5ad")
adataClean=adata2[~(adata2.obs['SubCT'].isna() | (adata2.obs['SubCT'] == 'Doublet'))].copy()

B=adataClean[adataClean.obs.SubCT.isin(["B",'Prolif.B'])]
pd.crosstab(B.obs.SubCT,B.obs.GSE)
B.obs['PatientDis'] = B.obs['Info_SampleType'].astype(str) + '_' + B.obs['Info_Patient'].astype(str)
B.obs['PatientDis']=B.obs['PatientDis'].astype('category')
B.obs.PatientDis.value_counts()

adata_conc=scvi.data.read_h5ad(f"{SaveFolder}/adata_Thymus_Manu_concatenate.h5ad")
adata_conc.obs = adata_conc.obs.merge(CT, left_index=True, right_index=True, how='left')
adata_conc.obs['PatientDis'] = adata_conc.obs['Info_SampleType'].astype(str) + '_' + adata_conc.obs['Info_Patient'].astype(str)
adata_conc.obs['PatientDis']=adata_conc.obs['PatientDis'].astype('category')

meta=B.obs
rawCount=adata_conc[adata_conc.obs.index.isin(meta.index)]

savePath=f"{SaveFolder}/PseudoBulk"
CUTOFF = 25

for i in range(0,1):
    ct='Bcell'
    cells = meta.index
    patient_status = meta.loc[cells, ['PatientDis']].drop_duplicates() 
    patient_status = patient_status[['PatientDis']].drop_duplicates()
    print(patient_status.shape)
    samples = patient_status["PatientDis"].unique() 
    print(len(samples))
    sample_values = [] 
    sample_ncells = [] 
    filtered_samples = [] 
    n_cells = []
    for s in samples:
        s_cells = rawCount.obs_names.isin(cells) & (rawCount.obs["PatientDis"] == s) 
        if s_cells.sum() >= CUTOFF:
            sample_values.append(rawCount.X[s_cells, :].sum(axis=0).A[0])
            sample_ncells.append((rawCount.X[s_cells, :] > 0).sum(axis=0).A[0])
            filtered_samples.append(s)
            n_cells.append(s_cells.sum())
            
    sample_values = pd.DataFrame(sample_values, index=filtered_samples, columns=rawCount.var_names).T
    sample_ncells = pd.DataFrame(sample_ncells, index=filtered_samples, columns=rawCount.var_names).T

    sample_values.to_csv("%s/sampleValues_%s.txt" %(savePath,ct), sep="\t")
    sample_ncells.to_csv("%s/sample_NCells_%s.txt" %(savePath,ct), sep="\t")
    
    patient_status = patient_status.loc[patient_status["PatientDis"].isin(filtered_samples), :]
    print(n_cells)
    patient_status["n_cells"] = n_cells
    patient_status.to_csv("%s/sampleMeta_%s.csv" %(savePath,ct))
    print(f"{ct} done, {len(filtered_samples)}/{len(samples)}")


