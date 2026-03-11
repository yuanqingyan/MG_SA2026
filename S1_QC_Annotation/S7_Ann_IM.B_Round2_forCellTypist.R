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

saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"
list.files(saveFolder)

#################################################################
har_IM<-readRDS(file=sprintf("%s/har_IM_annotated.RDS",saveFolder))
unique(har_IM$CellType1)
har_IM<-subset(har_IM,CellType_Round1 %in% "Doublet",invert=T)
har_B<-subset(har_IM,CellType_Round1 %in% c("B",'Prolif.B','Plasma'))
DimPlot(har_B,group.by='CellType_Round1',label=T,repel=T,raster=T)
saveRDS(har_B,file=sprintf("%s/Thymoma.B.RDS",FigOut))

######################################################################################
harFun<-function(obj=obj,projN='Test',nf=3000,pcs=40,kP=60,res=1.2){
  rawCount_temp_In<-obj@assays$RNA@counts
  har_obj<-CreateSeuratObject(counts = rawCount_temp_In, project =projN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nf) %>% #
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = pcs, verbose = FALSE)
  har_obj<-AddMetaData(har_obj,metadata=obj@meta.data[,c(4:(ncol(obj@meta.data)))])
  har_obj<-har_obj %>% RunHarmony("ID", plot_convergence = FALSE)
  har_obj <- har_obj %>%
    RunUMAP(reduction = "harmony", dims = 1:pcs) %>%
    FindNeighbors(reduction = "harmony", dims = 1:pcs, k.param = kP) %>%
    FindClusters(resolution = res) %>%
    identity()
  return(har_obj)
}

#######################################################################################


har_B<-readRDS(file=sprintf("%s/Thymoma.B.RDS",FigOut))
har_B$organ<-ifelse(har_B$ID %in% c("MG03_PB","MG03_PT","MG22_PL"),'Blood','Tissue')

All.B0<-harFun(obj=har_B,projN='AllB',nf=3500,pcs=50,kP=50,res=0.2)
DimPlot(All.B0,label=T)+NoLegend()
DefaultAssay(All.B0)<-"RNA";Idents(All.B0)<-"RNA_snn_res.0.2"
FeaturePlot(All.B0,features=c('CD3E'),col=c('grey','red'))

All.B0<-FindSubCluster(All.B0,cluster='5',graph.name="RNA_snn",subcluster.name='SubC',resolution=0.1)
Idents(All.B0)<-'SubC'
DimPlot(All.B0,label=T,repel=T)+NoLegend()
DefaultAssay(All.B0)<-"RNA"
FeaturePlot(All.B0,features=c('CD3E'),col=c('grey','red'))

####################
###remove 5_0 (T); 5_1 (Macrophage??)
All.B0_clean<-subset(All.B0, SubC %in% c("5_0","5_1"), invert=T)
All.B<-harFun(obj=All.B0_clean,projN='AllB',nf=4500,pcs=50,kP=30,res=0.2)#

DefaultAssay(All.B)<-"RNA";Idents(All.B)<-"RNA_snn_res.0.2"
All.B<-FindSubCluster(All.B,cluster='7',graph.name="RNA_snn",
                      subcluster.name='C7',resolution=0.1)
Idents(All.B)<-'C7'
DimPlot(All.B,label=T)+NoLegend()

DefaultAssay(All.B)<-"RNA";Idents(All.B)<-"C7"
All.B<-FindSubCluster(All.B,cluster='11',graph.name="RNA_snn",subcluster.name='C11',resolution=0.2)
DefaultAssay(All.B)<-"RNA";Idents(All.B)<-'C11'
DimPlot(All.B,label=T)+NoLegend()
All.B<-FindSubCluster(All.B,cluster='0',graph.name="RNA_snn",subcluster.name='C0',resolution=0.2)
DefaultAssay(All.B)<-"RNA";Idents(All.B)<-'C0'
DimPlot(All.B,label=T)+NoLegend()
All.B<-FindSubCluster(All.B,cluster='2',graph.name="RNA_snn",subcluster.name='C2',resolution=0.4)
DefaultAssay(All.B)<-"RNA";Idents(All.B)<-'C2'
DimPlot(All.B,label=T)+NoLegend()

DefaultAssay(All.B)<-"RNA";Idents(All.B)<-'C2'
All.B<-FindSubCluster(All.B,cluster='4',graph.name="RNA_snn",subcluster.name='C4',resolution=0.06)
saveRDS(All.B,file=sprintf("%s/temp/AllB_withCluster.RDS",saveFolder))

