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


##### this is seurat object form Yausmizu dataset, will be used as reference ####
seu<-readRDS(file="./Projects/PubD/Thymoma/SeuPublicGSE/Yasumizu_Seu.RDS")
DimPlot(seu,group.by="cell_type__custom",label=T,repel=T)+NoLegend()
unique(seu$cell_type__custom)

BcellLabel<-c("Naive B cell", "Memory B cell (II)", "Pre GC B cell", "Plasmablast",
              "Unswitched memory B cell", "Memory B cell (I)","Thymic memory B cell","GC B cell")
Bcell<-subset(seu,cell_type__custom %in% BcellLabel)
DimPlot(Bcell,group.by="cell_type__custom",label=T,repel=T)+NoLegend()

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

toH5ad(obj=Bcell,
       GSEName='ReferenceB',
       saveFolder=saveFolder_temp)



##############################################################################
##############################################################################
##############################################################################

saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"
list.files(saveFolder)

All.B<-readRDS(file=sprintf("%s/temp/AllB_withCluster.RDS",saveFolder))
DimPlot(All.B,group.by='C4',label=T,repel=T)

toH5ad(obj=All.B,
       GSEName='Thymoma_All.B',
       saveFolder=saveFolder_temp)


