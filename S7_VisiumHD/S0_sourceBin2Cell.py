import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
import sys
import subprocess
import fnmatch
import glob
import bin2cell as b2c
import cv2
import tensorflow


help_menu='''\nUsage: python sourceBin2Cell.py -ID xxxx
**parameters**
'''

args=sys.argv
if '-h' in args or '-help' in args or len(args)==1:
  print (help_menu)
sys.exit(0)

if '-ID' not in args:
  sys.exit('SampleID needed')
else:
  i=sys.argv.index('-ID')
sampleName=sys.argv[i+1]


tempH5adPath='./Spatial/Visium_HD'
os.chdir(f'{tempH5adPath}/Output/')

os.makedirs("./output/B2C/stardist", exist_ok=True)
outPath="./output/B2C/stardist"
folderPath='./Apr2025/'
pathRaw='./VisiumHD-set2'

###############
SaveH5adPath="./Apr2025_H5ADRDS"
###############

path = f"{folderPath}/{sampleName}/outs/binned_outputs/square_002um"
source_image_path=glob.glob(pathname=f'{pathRaw}/Brightfield_images/{sampleName}*')[0]
spaceranger_image_path = f"{folderPath}/{sampleName}/outs/binned_outputs/square_002um/spatial"

adata = b2c.read_visium(path, 
                        source_image_path = source_image_path, 
                        spaceranger_image_path = spaceranger_image_path)
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_counts=1)

mpp = 0.5

b2c.scaled_he_image(adata, 
                    mpp=mpp, 
                    spatial_cropped_key="spatial_cropped150_buffer",
                    img_key="mpp0.5",
                    save_path=f"{outPath}/he_{sampleName}.tiff")
b2c.destripe(adata)


mask = ((adata.obs['array_row'] >= 1000) & 
          (adata.obs['array_row'] <= 1200) & 
          (adata.obs['array_col'] >= 1000) & 
          (adata.obs['array_col'] <= 1200)
)


bdata = adata[mask] 

sc.set_figure_params(figsize=[10,10],dpi=100)
sc.pl.spatial(bdata, 
              color=[None, "n_counts_adjusted"],
              alpha=0.5,
              # cmap='gist_rainbow',
              color_map="OrRd",
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              save=f'{sampleName}_destripe.pdf')

b2c.stardist(image_path=f"{outPath}/he_{sampleName}.tiff", 
             labels_npz_path=f"{outPath}/he_{sampleName}.npz", 
             stardist_model="2D_versatile_he", 
             prob_thresh=0.01)

b2c.insert_labels(adata, 
                  labels_npz_path=f"{outPath}/he_{sampleName}.npz", 
                  basis="spatial", 
                  spatial_key="spatial_cropped150_buffer",
                  mpp=mpp, 
                  labels_key="labels_he")

bdata = adata[mask]
bdata = bdata[bdata.obs['labels_he']>0]
bdata.obs['labels_he'] = bdata.obs['labels_he'].astype(str)
sc.pl.spatial(bdata, 
              color=[None, "labels_he"], 
              color_map="OrRd",
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              save=f'{sampleName}_afterHEmodel.pdf')


crop = b2c.get_crop(bdata, basis="spatial", spatial_key="spatial_cropped150_buffer", mpp=mpp)

rendered = b2c.view_stardist_labels(image_path=f"{outPath}/he_{sampleName}.tiff", 
                                    labels_npz_path=f"{outPath}/he_{sampleName}.npz", 
                                    crop=crop)
plt.imshow(rendered)
plt.savefig(f"./figures/{sampleName}_viewStarLabel.pdf", bbox_inches='tight')
plt.close()

crop2 = b2c.get_crop(bdata, basis="spatial", spatial_key="spatial_cropped150_buffer", mpp=mpp)

rendered2 = b2c.view_labels(image_path=f"{outPath}/he_{sampleName}.tiff", 
                            labels_npz_path=f"{outPath}/he_{sampleName}.npz", 
                            crop=crop2)
plt.imshow(rendered2)
plt.savefig(f"./figures/{sampleName}_viewLabel.pdf", bbox_inches='tight')
plt.close()

b2c.expand_labels(adata, 
                  labels_key='labels_he', 
                  expanded_labels_key="labels_he_expanded")
bdata = adata[mask]
bdata = bdata[bdata.obs['labels_he_expanded']>0]
bdata.obs['labels_he_expanded'] = bdata.obs['labels_he_expanded'].astype(str)

sc.pl.spatial(bdata, 
              color=[None, "labels_he_expanded"], 
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              save=f'{sampleName}_afterExpandLabel.pdf')      

b2c.grid_image(adata, 
               "n_counts_adjusted", 
               mpp=mpp, 
               sigma=5, 
               save_path=f"{outPath}/gex_{sampleName}.tiff")

sc.pl.spatial(bdata, 
              color=["labels_he_expanded"], 
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              save=f'{sampleName}_aftergridImage.pdf')      

b2c.stardist(image_path=f"{outPath}/gex_{sampleName}.tiff", 
             labels_npz_path=f"{outPath}/gex_{sampleName}.npz", 
             stardist_model="2D_versatile_fluo", 
             prob_thresh=0.05, 
             nms_thresh=0.5)

b2c.insert_labels(adata, 
                  labels_npz_path=f"{outPath}/gex_{sampleName}.npz", 
                  basis="array", 
                  mpp=mpp, 
                  labels_key="labels_gex")


bdata = adata[mask]
bdata = bdata[bdata.obs['labels_gex']>0]
bdata.obs['labels_gex'] = bdata.obs['labels_gex'].astype(str)
sc.pl.spatial(bdata, 
              color=[None, "labels_gex"], 
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              save=f'{sampleName}_afterLabelGex.pdf')


crop4 = b2c.get_crop(bdata, basis="array", mpp=mpp)

rendered4 = b2c.view_stardist_labels(image_path=f"{outPath}/gex_{sampleName}.tiff", 
                                     labels_npz_path=f"{outPath}/gex_{sampleName}.npz", 
                                     crop=crop4)
plt.imshow(rendered4)   
plt.savefig(f"./figures/{sampleName}_viewStarLabel_Gex.pdf", bbox_inches='tight')
plt.close()  

crop5 = b2c.get_crop(bdata, basis="array", spatial_key="spatial_cropped150_buffer", mpp=mpp)

rendered5 = b2c.view_labels(image_path=f"{outPath}/gex_{sampleName}.tiff", 
                            labels_npz_path=f"{outPath}/gex_{sampleName}.npz",
                            crop=crop5)
plt.imshow(rendered5)    
plt.savefig(f"./figures/{sampleName}_viewLabel_Gex.pdf", bbox_inches='tight')
plt.close()   

b2c.salvage_secondary_labels(adata, 
                             primary_label="labels_he_expanded", 
                             secondary_label="labels_gex", 
                             labels_key="labels_joint")


bdata = adata[mask]
bdata = bdata[bdata.obs['labels_joint']>0]
bdata.obs['labels_joint'] = bdata.obs['labels_joint'].astype(str)

sc.pl.spatial(bdata, 
              color=[None, "labels_joint_source", "labels_joint"], 
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              save=f'{sampleName}_afterLabelJoint.pdf')

cdata = b2c.bin_to_cell(adata, 
                        labels_key="labels_joint", 
                        spatial_keys=["spatial", "spatial_cropped150_buffer"])


cell_mask = ((cdata.obs['array_row'] >= 1000) & 
               (cdata.obs['array_row'] <= 1200) & 
               (cdata.obs['array_col'] >= 1000) & 
               (cdata.obs['array_col'] <= 1200))

ddata = cdata[cell_mask]
sc.pl.spatial(ddata, 
              color=["bin_count", "labels_joint_source"], 
              img_key="mpp0.5", 
              basis="spatial_cropped150_buffer",
              save=f'{sampleName}_withSingleCell.pdf')


cdata.write_h5ad(f'{tempH5adPath}/H5AD/{sampleName}.b2c.h5ad')
adata.write_h5ad(f'{tempH5adPath}/H5AD/{sampleName}.beforeJoin.h5ad')
cdata.write_h5ad(f'{SaveH5adPath}/B2C_H5AD/{sampleName}.b2c.h5ad')
adata.write_h5ad(f'{SaveH5adPath}/B2C_H5AD/{sampleName}.beforeJoin.h5ad')

