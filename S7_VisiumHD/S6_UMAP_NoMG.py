import scanpy as sc
from banksy.initialize_banksy import initialize_banksy
from banksy.run_banksy import run_banksy_multiparam
import pandas as pd
import matplotlib.pyplot as plt
from copy import copy
reds = copy(plt.cm.Reds)

SaveFolderPath="./Apr2025_H5ADRDS/CellTypist_H5AD"
FigOut="./output/QC_Fig"

sampleIn=['6-Thymus','17-Thymus', '12-thymoma']
allAda=[sc.read_h5ad(f"{SaveFolderPath}/{i}_celltypist.H5AD") for i in sampleIn]
len(allAda)

for i in range(len(sampleIn)):
  allAda[i].obs['SampleName']=sampleIn[i]
allAda[0].obs.head(3) 


cutoff=0.7
Dat=[allAda[0][allAda[0].obs.ConfScore>=cutoff],
     allAda[1][allAda[1].obs.ConfScore>=cutoff],
     allAda[2][allAda[2].obs.ConfScore>=cutoff]]
len(Dat)


rename_dict = {
  "DPT":"DPT",
  "CD4.Naive":"CD4T",
  "mcTEC":"TEC",
  "MemB.Rib+Mit+":"B",
  "Plasma":"Plasma",
  "Macrophage":"Macrophage",
  "Prolf.T":"Prolif.T",
  "KRT15+.mTEC":"TEC",
  "Fibroblast":"Fibroblast",
  "NaivB":"B",
  "VSMC":"VSMC",
  "CD8.Tem":"CD8T",
  "CD8a/b(entry)":"CD8T",
  "CD8T.Naive":"CD8T",
  "MemB.Rib+":"B",
  "T(agonist)":"T(agonist)",
  "SLPI+.mTEC":"TEC",
  "DNT":"DNT",
  "Activated.En":"En",
  "Vascular":"Vascular",
  "Endo":"En",
  "CD8a/a":"CD8T",
  "CD4.Tcm(Th2)":"CD4T",
  "DC":"DC",
  "Pericyte":"Pericyte",
  "Treg":"CD4T",
  "NK":"NK",
  "MemB":"B",
  "CD4.Tcm(Th0)":"CD4T",
  "CD4.Tcm(Th17)":"CD4T",
  "cTEC":"TEC",
  "CD8.Temra":"CD8T",
  "CD8.Tcm":"CD8T",
  "CRTAM+.gdT":"gdT",
  "gdT":"gdT",
  "MemB.Mit+":"B",
  "CD4.Tem(Th1)":"CD4T",
  "Tfh":"CD4T",
  "MemB.ISG":"B",
  "CD4.Tem(Th1/Th17)":"CD4T",
  "Neutrophils":"Neutrophils",
  "ILC":"ILC",
  "Lymphatic":"En",
  'Switch':'B'
}


Dat[0].obs["NewCT"] = [rename_dict[elem] if elem in rename_dict else elem for elem in Dat[0].obs["PredLabel"]]
Dat[1].obs["NewCT"] = [rename_dict[elem] if elem in rename_dict else elem for elem in Dat[1].obs["PredLabel"]]

pd.crosstab(Dat[0].obs.NewCT, Dat[0].obs.PredLabel)


adata17=Dat[1].copy()

sc.pp.highly_variable_genes(adata17, 
                            n_top_genes=2000, 
                            flavor='seurat')

coord_keys = ('x', 'y', 'spatial_cropped150_buffer') 

adata17.obs['x']=adata17.obsm['spatial_cropped150_buffer'][:,0]
adata17.obs['y']=adata17.obsm['spatial_cropped150_buffer'][:,1]

banksy_dict = initialize_banksy(
  adata17, 
  coord_keys=coord_keys,  
  num_neighbours=15,
  nbr_weight_decay='scaled_gaussian'
)

banksy_dict

from banksy_utils.color_lists import spagcn_color

resolutions = [0.3]  # 
pca_dims = [15]  # 
lambda_list = [0.2]  #
max_m=1


adata17.obs['NewCT'] = adata17.obs['NewCT'].astype('category')

results_df = run_banksy_multiparam(
  adata17,
  banksy_dict,
  lambda_list=lambda_list,# Use mixed for best results; # lambda_ = 0.0 (nonspatial), 0.5 (mixed), 1.0 (spatial only)
  resolutions=resolutions,
  
  color_list = spagcn_color,
  max_m = max_m,
  filepath = FigOut,
  
  key = coord_keys,
  
  pca_dims = [10],
  annotation_key = 'NewCT',
  match_labels = False,
  savefig = False,
  add_nonspatial = False,
  variance_balance = False,
)

banksy_dict


from banksy.embed_banksy import generate_banksy_matrix


banksy_dict, banksy_matrix = generate_banksy_matrix(adata17,
                                                    banksy_dict,
                                                    lambda_list,
                                                    max_m)

banksy_dict


from banksy.main import concatenate_all

banksy_dict["nonspatial"] = {
  0.0: {"adata": concatenate_all([adata17.X], 0, adata=adata17), }
}

print(banksy_dict['nonspatial'][0.0]['adata'])


from banksy_utils.umap_pca import pca_umap

pca_umap(banksy_dict,
         pca_dims = pca_dims,
         add_umap = True,
         plt_remaining_var = False,
)

from banksy.cluster_methods import run_Leiden_partition

seed=5

results_df, max_num_labels = run_Leiden_partition(
  banksy_dict,
  resolutions,
  num_nn = 50, # 
  num_iterations = -1, # 
  partition_seed = seed,
  match_labels = True,
)

c_map ='tab20'
banksy_adata = banksy_dict['scaled_gaussian'][0.2]['adata']

banksy_adata.obsm['X_umap'] = banksy_adata.obsm['reduced_pc_15_umap']

banksy_adata.write(f"{SaveFolderPath}/UMAP_for_thymus17_fromBanksy_hiM.h5ad")


#####################
rename_dict2 = {
  "DPT":"T/NK",
  "CD4.Naive":"T/NK",
  "mcTEC":"Ep",
  "MemB.Rib+Mit+":"B/plasma",
  "Plasma":"B/plasma",
  "Macrophage":"Mye",
  "Prolf.T":"T/NK",
  "KRT15+.mTEC":"Ep",
  "Fibroblast":"Stromal",
  "NaivB":"B/plasma",
  "VSMC":"Stromal",
  "CD8.Tem":"T/NK",
  "CD8a/b(entry)":"T/NK",
  "CD8T.Naive":"T/NK",
  "MemB.Rib+":"B/plasma",
  "T(agonist)":"T/NK",
  "SLPI+.mTEC":"Ep",
  "DNT":"T/NK",
  "Activated.En":"En",
  "Vascular":"En",
  "Endo":"En",
  "CD8a/a":"T/NK",
  "CD4.Tcm(Th2)":"T/NK",
  "DC":"Mye",
  "Pericyte":"Stromal",
  "Treg":"T/NK",
  "NK":"T/NK",
  "MemB":"B/plasma",
  "CD4.Tcm(Th0)":"T/NK",
  "CD4.Tcm(Th17)":"T/NK",
  "cTEC":"Ep",
  "CD8.Temra":"T/NK",
  "CD8.Tcm":"T/NK",
  "CRTAM+.gdT":"T/NK",
  "gdT":"T/NK",
  "MemB.Mit+":"B/plasma",
  "CD4.Tem(Th1)":"T/NK",
  "Tfh":"T/NK",
  "MemB.ISG":"B/plasma",
  "CD4.Tem(Th1/Th17)":"T/NK",
  "Neutrophils":"Mye",
  "ILC":"Mye",
  "Lymphatic":"En",
  'Switch':'B/plasma',
  "CD8.Trm":"T/NK"
}


banksy_adata.obs["NewCT2"] = [rename_dict2[elem] if elem in rename_dict2 else elem for elem in banksy_adata.obs["PredLabel"]]
banksy_adata.obs["NewCT2"]=pd.Categorical(banksy_adata.obs["NewCT2"],
                                          categories=['T/NK','B/plasma','En','Ep','Stromal','Mye'],
                                          ordered=True)

############## Fig. S9 #####

plt.rcParams["figure.figsize"] = (8, 7)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
sc.pl.umap(banksy_adata, 
           color='NewCT2',
           frameon=False,
           add_outline=False,
           show=False,
           size=2, 
           palette="tab10",
           legend_fontsize='small',
           legend_fontweight='light')
plt.savefig(f"{FigOut}/umap_NoMG_label.pdf")


plt.rcParams['figure.dpi']=600
plt.rcParams["figure.figsize"] = (8, 7)
sc.pl.umap(banksy_adata, 
           color='NewCT2',
           frameon=False,
           add_outline=False,
           show=False,
           size=2, 
           palette="tab10",
           legend_fontsize='small',
           legend_fontweight='light')
plt.savefig(f"{FigOut}/umap_NoMG_noLabel.tiff")


geneIn2=['MS4A1','CD70','CD83','ADAM28']

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
sc.pl.umap(banksy_adata, 
           color=geneIn2,
           cmap=reds,
           frameon=False,
           add_outline=False,
           show=False,
           ncols=2,
           size=3,
           vmax=1.5,
           vmin=0.2)
plt.savefig(f"{FigOut}/PaperGene_NoMG.pdf")


plt.rcParams['figure.dpi']=600
sc.pl.umap(banksy_adata, 
           color=geneIn2,
           cmap=reds,
           frameon=False,
           add_outline=False,
           show=False,
           ncols=2,
           size=3,
           vmax=1.5,
           vmin=0.2)
plt.savefig(f"{FigOut}/PaperGene_NoMG.tiff")





