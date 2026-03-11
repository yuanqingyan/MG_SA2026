import pandas as pd
import scanpy as sc
import plotnine as p9
import torch
import tensorly as tl
import liana as li
import cell2cell as c2c
import decoupler as dc # needed for pathway enrichment
import pickle
import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict

use_gpu = True

if use_gpu:
  device = "cuda" if torch.cuda.is_available() else "cpu"
if device == "cuda":
  tl.set_backend('pytorch')
else:
  device = "cpu"

print(device)

sampleName="TB_MGOnly_June2025"

savePath="./Projects/ManuscriptData/IntegrateThymus/data/H5AD/BT_Interaction"
adata = sc.read_h5ad(f"{savePath}/seuForCC_MG.h5ad")
adata.obs["CTCC"].cat.categories

adata.obs['cell_type']=adata.obs['CTCC'].astype('category')
adata.obs['cell_type'].cat.categories

adata1=adata[adata.obs.cell_type.isin(['Naiv_Cir_1','Naiv_Cir_2', 'Naiv_Cir_3', 'Naiv_Cir_4', 'Naiv_Cir_5',
adata1.obs['cell_type'].cat.categories
adata1.obs['cell_type'].value_counts()

adata1.obs['Disease_Fig'].cat.categories
adata1.obs['patient'].cat.categories

adata1.obs['PtDis']=adata1.obs['patient'].astype('str')+'_'+adata1.obs['Disease2'].astype('str')
adata1.obs['PtDis']=adata1.obs['PtDis'].astype('category')
adata1.obs['PtDis'].cat.categories
adata1.obs['Study'].cat.categories


#################################################################################
target0=adata1[adata1.obs.Study.isin(["NU"])]

cCount=target0.obs['cell_type'].value_counts()
sel1=cCount[cCount>200]
target=target0[target0.obs.cell_type.isin(sel1.index)]
#################################################################################
sample_key = 'PtDis'
condition_key = 'Disease2'
groupby = 'cell_type'


li.mt.rank_aggregate.by_sample(
  target,
  groupby=groupby,
  resource_name='consensus', 
  sample_key=sample_key, 
  use_raw=False,
  verbose=True, 
  n_perms=500, 
  return_all_lrs=True, 
)
target.uns["liana_res"].sort_values("magnitude_rank").head(10)
target.uns["liana_res"].sort_values("magnitude_rank").tail(10)
target.write_h5ad(f"{savePath}/{sampleName}_addLiana.h5ad")

tensor = li.multi.to_tensor_c2c(target,
                                sample_key=sample_key,
                                score_key='magnitude_rank', 
                                how='outer_cells' 
)

c2c.io.export_variable_with_pickle(tensor, f"{savePath}/{sampleName}_tensor_pre.pkl")

context_dict = target.obs[[sample_key, condition_key]].drop_duplicates()
context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
context_dict = defaultdict(lambda: 'Unknown', context_dict)

tensor_meta = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                  metadata_dicts=[context_dict, None, None, None],
                                                  fill_with_order_elements=True
)
tensor_meta

tensor = c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                                    tensor_meta,
                                                    copy_tensor=True, 
                                                    rank=None, 
                                                    tf_optimization='regular', 
                                                    random_state=0, 
                                                    device=device, 
                                                    elbow_metric='error', 
                                                    smooth_elbow=False, 
                                                    upper_rank=30, 
                                                    tf_init='random', 
                                                    tf_svd='numpy_svd', 
                                                    cmaps=None, 
                                                    sample_col='Element', 
                                                    group_col='Category', 
                                                    output_fig=False,  
)

c2c.io.export_variable_with_pickle(tensor, f"{savePath}/{sampleName}_tensor_gpu.pkl")


