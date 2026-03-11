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
import matplotlib.pyplot as plt
#plt.rcParams["font.family"] = "Arial"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams["figure.figsize"] = (10, 10)

from copy import copy
reds = copy(plt.cm.Reds)
reds.set_under("lightgray")

###%matplotlib inline

# NOTE: to use CPU instead of GPU, set use_gpu = False
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
FigOut="./Thymus/FigOut/C2C"



target=sc.read_h5ad(f"{savePath}/{sampleName}_addLiana.h5ad")

sample_key = 'PtDis'
condition_key = 'Disease2'
groupby = 'cell_type'

target.obs[condition_key].value_counts()

###Build Metadata
context_dict = target.obs[[sample_key, condition_key]].drop_duplicates()
context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
context_dict = defaultdict(lambda: 'Unknown', context_dict)
context_dict

with open(f"{savePath}/{sampleName}_tensor_pre.pkl", 'rb') as f:
  tensor0 = pickle.load(f)

tensor_meta = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor0,
                                                  metadata_dicts=[context_dict, None, None, None],
                                                  fill_with_order_elements=True)
tensor_meta

with open(f"{savePath}/{sampleName}_tensor_gpu.pkl", 'rb') as f:
  tensor = pickle.load(f)

factors, axes = c2c.plotting.tensor_factors_plot(interaction_tensor=tensor,
                                                 metadata = tensor_meta, # This is the metadata for each dimension
                                                 sample_col='Element',
                                                 group_col='Category',
                                                 meta_cmaps = ['viridis', 'Dark2_r', 'tab20', 'tab20'],
                                                 fontsize=10, # Font size of the figures generated
)

plt.savefig(f"{FigOut}/{sampleName}_tensorBarplot.pdf",bbox_inches='tight')


##################Fig. 4A ###########
factors = tensor.factors
_ = c2c.plotting.context_boxplot(context_loadings=factors['Contexts'],
                                 metadict=context_dict,
                                 nrows=3,
                                 figsize=(18, 16),
                                 statistical_test='Mann-Whitney',#'t-test_ind';'Mann-Whitney'
                                 pval_correction='BH',#'fdr_bh'
                                 cmap='tab10',
                                 verbose=False,
)
plt.savefig(f"{FigOut}/{sampleName}_context_Boxplot.pdf",bbox_inches='tight')


################# Fig. 4B ##########
c2c.plotting.ccc_networks_plot(factors,
                               included_factors=['Factor 1'],
                               network_layout='spring',
                               ccc_threshold=0.0646, #
                               nrows=1,
                               panel_size=(6, 6), #
                               
                               edge_width=5,
                               edge_arrow_size=8,
                               edge_alpha=0.8,
                               edge_color='blue',
                               
                               node_size=8,
                               node_label_size=10,
                               node_label_alpha=0.9,
                               node_color='orange',
                               node_label_offset=(0,0),
                               
                               factor_title_size=8
)

plt.savefig(f"{FigOut}/{sampleName}_network_Factor1_Fig4B.pdf",bbox_inches='tight')

c2c.plotting.ccc_networks_plot(factors,
                               included_factors=['Factor 4'],
                               network_layout='spring',
                               ccc_threshold=0.0646, #
                               nrows=1,
                               panel_size=(6, 6), #
                               
                               edge_width=5,
                               edge_arrow_size=8,
                               edge_alpha=0.8,
                               edge_color='green',
                               
                               node_size=8,
                               node_label_size=10,
                               node_label_alpha=0.9,
                               node_color='orange',
                               node_label_offset=(0,0),
                               
                               factor_title_size=8
)

plt.savefig(f"{FigOut}/{sampleName}_network_Factor4_Fig4B.pdf",bbox_inches='tight')

c2c.plotting.ccc_networks_plot(factors,
                               included_factors=['Factor 8'],
                               network_layout='spring',
                               ccc_threshold=0.049, # 
                               nrows=1,
                               panel_size=(6, 6), # 
                               
                               edge_width=5,
                               edge_arrow_size=8,
                               edge_alpha=0.8,
                               edge_color='purple',
                               
                               node_size=8,
                               node_label_size=10,
                               node_label_alpha=0.9,
                               node_color='orange',
                               node_label_offset=(0,0),
                               
                               factor_title_size=8
)

plt.savefig(f"{FigOut}/{sampleName}_network_Factor8_Fig4B.pdf",bbox_inches='tight')


lr_loadings = factors['Ligand-Receptor Pairs']
lr_loadings.to_csv(f'{FigOut}/TensorC2C_lr_loadings.csv')
lr_loadings.sort_values("Factor 4", ascending=False).head(10)

Send_loadings = factors['Sender Cells']
Send_loadings.to_csv(f'{FigOut}/TensorC2C_SenderCells.csv')
Send_loadings.sort_values("Factor 4", ascending=False).head(50)

Rec_loadings = factors['Receiver Cells']
Rec_loadings.to_csv(f'{FigOut}/TensorC2C_ReceiverCells.csv')
Rec_loadings.sort_values("Factor 4", ascending=False).head(50)




