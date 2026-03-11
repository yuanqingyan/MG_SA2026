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
list.files(saveFolder)


All.B<-readRDS(file=sprintf("%s/B_withSubtype_Annotated_May2025.RDS",saveFolder))
All.B$ptDis2<-sprintf("%s_%s",All.B$patient,All.B$Disease2);head(All.B)

NU_B0<-subset(All.B, (Study %in% 'NU') & (organ %in% 'Tissue'))
NU_B1<-subset(NU_B0, Bsub %in% c("Naiv_Cir_2", "Naiv_Cir_3","Naiv_Cir_5",
                                 "Naiv_Cir_1","Naiv_Cir_4","MemB_Cir",
                                 "Centroblast_GC","Centrocyte_GC",
                                 "Early_GC","Mature_GC","MemB.Mit+",'MemB.EBI3+'), invert=T)

#### to h5ad ########
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


NU_B1@meta.data<-NU_B1@meta.data[,!colnames(NU_B1@meta.data) %in% c("Info_Sample")]
toH5ad(obj=NU_B1,
       GSEName='forPseudoBulk_Bcell',
       saveFolder=FigOut)


