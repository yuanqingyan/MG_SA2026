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
from scipy import io
import subprocess
import math
import celltypist
from celltypist import models
from scipy.sparse import csr_matrix

print("Last run with scvi-tools version:", scvi.__version__)

adataAll=scvi.data.read_h5ad("./Projects/ManuscriptData/IntegrateThymus/data/H5AD/seuThymoma___withAllCellType.h5ad")
adataAll


mpl.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adataAll,
           color="FinalCT",
           frameon=False,
           add_outline=False,
           show=False,
           palette="tab20",
           size=2, 
           legend_loc="on data")

def makeCellTypistModel(adataInTest, modelName='test', modelPath='test', tempOut='test',cell_type_column='FinalCT'):
    os.makedirs(modelPath, exist_ok=True)
    os.makedirs(tempOut, exist_ok=True)
    
    matrix = csr_matrix(adataInTest.X)
    io.mmwrite(f'{tempOut}/{modelName}_forCellTypist.mtx', matrix) 
    
    with open(f'{tempOut}/{modelName}_gene.tsv', 'w') as f:
        for item in adataInTest.var_names:
            f.write(item + '\n')
     reference_adata = adataInTest.copy()
    
    model=celltypist.train(X=f'{tempOut}/{modelName}_forCellTypist.mtx', 
                 genes=f"{tempOut}/{modelName}_gene.tsv",
                 labels=reference_adata.obs[cell_type_column],
                 use_GPU=True,
                 n_jobs = 30, 
                 max_iter = 100)
    model.write(f"{modelPath}/{modelName}.pkl")


ThymusThymoma=adataAll[adataAll.obs.Disease_Fig.isin(['Thymoma_MG','Thymoma_no_MG','Thymus_MG','Thymus_no_MG'])]
NUThy=ThymusThymoma[ThymusThymoma.obs.Study.isin(['NU'])]

makeCellTypistModel(adataInTest=NUThy.copy(), 
                    modelName='NUThy', 
                    modelPath='./Projects/ManuscriptData/IntegrateThymus/data/H5AD', 
                    tempOut='./output/CellTypistModel',
                    cell_type_column='FinalCT')


