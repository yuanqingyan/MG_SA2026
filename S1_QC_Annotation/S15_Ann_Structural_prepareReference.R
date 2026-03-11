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
#options(Seurat.object.assay.version = 'v4')
library(Matrix)

##### this is seurat object form Yausmizu dataset, will be used as reference ####
seu<-readRDS(file="./Projects/PubD/Thymoma/SeuPublicGSE/Yasumizu_Seu.RDS")
EPcellLabel<-c("mTEC (I)","mTEC (II)","cTEC","nmTEC")
EPcell<-subset(seu,cell_type__custom %in% EPcellLabel)
Idents(EPcell)<-"cell_type__custom"



################################################################################################
##################               output anndata format       ###################################
################################################################################################

saveFolder_temp<-"./Thymus/h5ad"

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

toH5ad(obj=EPcell,
       GSEName='ThymomaEP',
       saveFolder=saveFolder_temp)



##############################################################################

saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"

#################################################################
har_EP2<-readRDS(file=sprintf("%s/Thymoma.har_EP2_forFig.RDS",saveFolder))

toH5ad(obj=har_EP2,
       GSEName='NUThymusEP',
       saveFolder=saveFolder_temp)
