library(Seurat)
library(harmony)
library(cowplot)
library(coin)
library(RColorBrewer)
library(scales)
require(gridExtra)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(msigdbr);msigdbr_species()
library(org.Hs.eg.db)
library(plotrix)
library(tidyr)
library(reticulate)
library(anndata)
library(sceasy)
options(Seurat.object.assay.version = 'v4')
library(Matrix)


saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"

####################################################
#############    for pseudobulk   ##################
####################################################
harT2<-readRDS(file=sprintf("%s/annotated_T_NK_ILC__June2025.RDS",saveFolder))
harT2$ptDis2<-sprintf("%s_%s",harT2$patient,harT2$Disease2)

toH5ad<-function(obj=obj,GSEName=GSEName,saveFolder=saveFolder_temp){
  matrix_X<-t(obj@assays$RNA@counts)
  obs_in<-obj@meta.data
  varN<-data.frame(name=colnames(matrix_X))
  varN$name<-gsub("-",".",varN$name)
  row.names(varN)<-varN$name
  colnames(matrix_X)<-varN$name
  matrixFinal<-as(matrix_X, "CsparseMatrix")
  umap_coordinates <- Embeddings(obj, reduction = "umap")
  pca_coordinates <- Embeddings(obj, reduction = "pca")

  anndata <- AnnData(
    X = matrixFinal,
    obs = obs_in,
    var = varN)

  anndata$obsm[['X_pca']] <- pca_coordinates
  anndata$obsm[["X_umap"]] <- umap_coordinates

  write_h5ad(anndata, file=sprintf("%s/%s.h5ad",saveFolder,GSEName))
}

toH5ad(obj=harT2,
       GSEName='forPseudoBulkR',
       saveFolder=FigOut)

