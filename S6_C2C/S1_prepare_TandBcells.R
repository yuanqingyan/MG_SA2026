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



All.B<-readRDS(file=sprintf("%s/B_withSubtype_Annotated_May2025.RDS",saveFolder))

All.B$Disease2<-dplyr::case_when(
  All.B$Disease_Fig %in% "Others"~'Others',
  All.B$Disease_Fig %in% c("Thymoma_MG","Thymus_MG") ~"MG",
  All.B$Disease_Fig %in% c("Thymoma_no_MG","Thymus_no_MG")~'NoMG',
  .default="unknnown"
)
saveColB<-c("Study",'Disease','patient','patientDisease','patientDisease2','CellType_Round1','Info_Sorting','Info_Patient',
            'Info_SampleType','Info_Sample','Disease_Fig','organ','Bsub','BigC','Disease2')
metaB<-All.B@meta.data[,saveColB]
metaB$whichTB<-"B"
colnames(metaB)[which(colnames(metaB)=="Bsub")]<-'SubCT'



harT2<-readRDS(file=sprintf("%s/annotated_T_NK_ILC__June2025.RDS",saveFolder))
saveColT<-c("Study",'Disease','patient','patientDisease','patientDisease2','CellType_Round1','Info_Sorting','Info_Patient',
            'Info_SampleType','Info_Sample','Disease_Fig','organ','TSub_F','Disease2')## 
metaT<-harT2@meta.data[,saveColT];
colnames(metaT)[which(colnames(metaT)=="TSub_F")]<-'SubCT'
metaT$BigC<-metaT$SubCT
metaT$whichTB<-"T"
allSelColN<-c("Study",'Disease','patient','patientDisease','patientDisease2','CellType_Round1','Info_Sorting','Info_Patient',
              'Info_SampleType','Info_Sample','Disease_Fig','organ','SubCT','BigC','Disease2','whichTB')
allSelColN %in% colnames(metaT)
allSelColN %in% colnames(metaB)

metaT<-metaT[,allSelColN]
metaB<-metaB[,allSelColN]

sum(colnames(metaT)==colnames(metaB))
metaTB<-rbind(metaB,metaT)
metaTB$CellID<-row.names(metaTB)

######################################################################################################################
seu0<-readRDS(file=sprintf("%s/seuThymoma.RDS",saveFolder))
meta0<-seu0@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA","percent.mt",'GSE','ID','Info',
                         'Population','Population2')]
meta0$CellID<-row.names(meta0)

library(dplyr)
metaMerg<-dplyr::left_join(meta0,metaTB,by='CellID')
row.names(metaMerg)<-metaMerg$CellID

sum(metaMerg$CellID==row.names(seu0@meta.data));dim(seu0)
seu0@meta.data<-metaMerg

######################################################################################################################
seuForCC0<-subset(seu0,whichTB %in% c("T",'B'));dim(seuForCC0)
seuForCC<-subset(seuForCC0, SubCT %in% c("NK",'ILC'),invert=T)

seuForCC$CTCC<-seuForCC$SubCT
seuForCC<-subset(seuForCC,organ %in% c('Tissue'))

seuForCC<-subset(seuForCC, Disease2 %in% c("MG",'NoMG'))

saveRDS(seuForCC,file=sprintf("%s/seuForCC.RDS",saveFolder))

######################################################################################################################
######################################################################################################################
######################################################################################################################
saveFolder_temp<-"./Projects/ManuscriptData/IntegrateThymus/data/H5AD/BT_Interaction"

toH5ad<-function(obj=obj,GSEName=GSEName,saveFolder=saveFolder_temp){
  matrix_X<-t(obj@assays$RNA@counts)
  obs_in<-obj@meta.data
  varN<-data.frame(name=colnames(matrix_X))
  varN$name<-gsub("-",".",varN$name)
  row.names(varN)<-varN$name
  colnames(matrix_X)<-varN$name
  matrixFinal<-as(matrix_X, "CsparseMatrix")
  umap_coordinates <- Embeddings(obj, reduction = "umap")

  anndata <- AnnData(
    X = matrixFinal,
    obs = obs_in,
    var = varN)
  
  anndata$obsm[["X_umap"]] <- umap_coordinates
  
  write_h5ad(anndata, file=sprintf("%s/%s.h5ad",saveFolder,GSEName))
}


toH5ad(obj=seuForCC,
       GSEName='seuForCC_MG',
       saveFolder=saveFolder_temp)

