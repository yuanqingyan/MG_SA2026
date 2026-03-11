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


saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"

har_EP2<-readRDS(file=sprintf("%s/Thymoma.har_EP2_forFig.RDS",saveFolder))

celltypist_pre<-read.csv("./Thymus/h5ad/ThymomaEP_EP_predictFromCellTypist.csv",row.names=1)
sum(row.names(har_EP2@meta.data)==row.names(celltypist_pre))#

har_EP2 <- AddMetaData(
  object = har_EP2,
  metadata = celltypist_pre[,(ncol(celltypist_pre)-7):ncol(celltypist_pre)] 
)
saveRDS(har_EP2, file=sprintf("%s/Thymoma.har_EP2_forFig__addCellTypistPred.RDS",saveFolder))



