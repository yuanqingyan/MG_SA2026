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

##### this is seurat object form Yausmizu dataset, will be used as reference ####
seu<-readRDS(file="./Projects/PubD/Thymoma/SeuPublicGSE/Yasumizu_Seu.RDS")
DimPlot(seu,group.by="cell_type__custom",label=T,repel=T)+NoLegend()

TcellLabel<-c("CD4 Temra (Th1)","CD8 Temra","CD4 Tnaive","CD4 Tem (Th1/17)","CD4 Tcm (Th0)",
              "DN T cell", "CD4 Tem (Th1)","CD4 Tcm (Tfh)","CD4 Tcm (Th2)","CD8 Tnaive",
              "CD4 Tcm (Th17)","Naive Treg","T agonist","Thymic CD4 T cell (II)", "Activated Treg",
              "CD8 Tem", "Thymic CD8 T cell", "Thymic CD4 T cell (I)","NKT cell (periphery)","NKT cell (thymus)",
              "gd T cell","DP T cell","NK cell","aa CD8 T cell (I)","aa CD8 T cell (II)",
              "CD8 Trm","Cycling DN/DP T cell")
Tcell<-subset(seu,cell_type__custom %in% TcellLabel)
Idents(Tcell)<-"cell_type__custom"
MarkerDiffT<-FindAllMarkers(Tcell,
                            slot='data', 
                            test.use='wilcox', 
                            only.pos=T)
FigOut="./Thymus/FigOut"
write.csv(MarkerDiffT,file=sprintf("%s/MarkerDiffT.csv",FigOut))


####################

library(reticulate)
library(anndata)
library(sceasy)
options(Seurat.object.assay.version = 'v4')
library(Matrix)


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

toH5ad(obj=Tcell,
       GSEName='ThymomaT',
       saveFolder=saveFolder_temp)



##############################################################################

saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"

#################################################################
har_IM<-readRDS(file=sprintf("%s/har_IM_annotated.RDS",saveFolder))
har_IM<-subset(har_IM,CellType_Round1 %in% "Doublet",invert=T)
DimPlot(har_IM,group.by='CellType_Round1',label=T,repel=T,raster=T)
har_T<-subset(har_IM,CellType_Round1 %in% c("B",'Prolif.B','Plasma','pDC','DC','Neutrophils','Macrophage'), invert=T)
saveRDS(har_T,file=sprintf("%s/T.Pop.RDS",FigOut))


toH5ad(obj=har_T,
       GSEName='NUThymusT',
       saveFolder=saveFolder_temp)



