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
seu0<-readRDS(file=sprintf("%s/seuThymoma.RDS",saveFolder))
temp_IM<-subset(seu0,Population2 %in% c("Structual"),invert=T)
saveRDS(temp_IM,file=sprintf("%s/temp_IM.RDS",FigOut))


har_IM<-CreateSeuratObject(counts = rawCount_temp_IM, project ="IM", min.cells = 3) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3500) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 38, verbose = FALSE)
har_IM<-AddMetaData(har_IM,metadata=temp_IM@meta.data[,c(5:ncol(temp_IM@meta.data))])
table(har_IM$ID)
har_IM<-har_IM %>% RunHarmony("ID", plot_convergence = FALSE)
har_IM <- har_IM %>%
  RunUMAP(reduction = "harmony", dims = 1:38) %>%
  FindNeighbors(reduction = "harmony", dims = 1:38, k.param = 60) %>%
  FindClusters(resolution = 0.8) %>%
  identity()
saveRDS(har_IM,file=sprintf("%s/har_IM.RDS",FigOut))


################################################################
har_IM<-readRDS(file=sprintf("%s/har_IM.RDS",FigOut))
DefaultAssay(har_IM)<-'RNA';Idents(har_IM)<-"RNA_snn_res.0.8"
IM_marker<-FindAllMarkers(har_IM, slot='data', test.use='wilcox', only.pos=T)
IM_marker_sort<-IM_marker[order(-IM_marker$avg_log2FC),]
IM_marker %>% group_by(cluster) %>% top_n(5, avg_log2FC) %>% as.data.frame()
head(IM_marker_sort[IM_marker_sort$cluster=="21",],20)
write.csv(IM_marker_sort,file=sprintf("%s/IMMarker.csv",FigOut))

marker18<-FindMarkers(har_IM, slot='data', ident.1="18", test.use='wilcox', only.pos=T)
marker18<-marker18[order(-marker18$avg_log2FC),]
#############
Idents(har_IM)<-"RNA_snn_res.0.8"
har_IM<-FindSubCluster(har_IM, cluster="18", graph.name="RNA_snn", subcluster.name = "C18", resolution = 0.1, algorithm = 1)
Idents(har_IM)<-"C18"
har_IM<-FindSubCluster(har_IM, cluster="17", graph.name="RNA_snn", subcluster.name = "C17", resolution = 0.1, algorithm = 1)
DimPlot(har_IM,group.by="C17",label=T,raster=F)
Idents(har_IM)<-"C17"
har_IM<-FindSubCluster(har_IM, cluster="12", graph.name="RNA_snn", subcluster.name = "C12", resolution = 0.5, algorithm = 1)
Idents(har_IM)<-"C12"
har_IM<-FindSubCluster(har_IM, cluster="11", graph.name="RNA_snn", subcluster.name = "C11", resolution = 0.1, algorithm = 1)
Idents(har_IM)<-"C11"

har_IM$CellType_Round1<-dplyr::case_when( 
  har_IM$C11 %in% c("18_0","18_2") ~"Doublet",
  har_IM$C11 %in% c("18_1") ~"Mast",
  har_IM$C11 %in% c("8",'9') ~"B",
  har_IM$C11 %in% c("22") ~"Prolif.B",
  har_IM$C11 %in% c("16") ~"Plasma",
  har_IM$C11 %in% c("17_1") ~"pDC",
  har_IM$C11 %in% c("17_0") ~"DC",
  har_IM$C11 %in% c("12_5") ~"Neutrophils",
  har_IM$C11 %in% c("12_0","12_1","12_2","12_3","12_4","12_6") ~"Macrophage",
  har_IM$C11 %in% c("15") ~"NK",
  har_IM$C11 %in% c("1",'10','14','20') ~"T",
  har_IM$C11 %in% c("3",'23') ~"T",
  har_IM$C11 %in% c('21') ~'T',
  har_IM$C11 %in% c('0','19','24','25') ~'T',
  har_IM$C11 %in% c("7") ~'T',
  har_IM$C11 %in% c("4") ~'T',##
  har_IM$C11 %in% c("13") ~'T',##
  har_IM$C11 %in% c("11_0","11_1") ~'T',##
  har_IM$C11 %in% c("2") ~'T',
  har_IM$C11 %in% c('5') ~'T',
  har_IM$C11 %in% c('6') ~'T',
  .default='Unknown'
);table(har_IM$CellType_Round1)
saveRDS(har_IM,file=sprintf("%s/har_IM_annotated.RDS",saveFolder))

################################################################################################
################################################################################################
################################################################################################
har_IM<-readRDS(file=sprintf("%s/har_IM_annotated.RDS",saveFolder))
har_IM<-subset(har_IM,CellType_Round1 %in% "Doublet",invert=T)

(DimPlot(har_IM,group.by='CellType_Round1',label=T,repel=T,raster=T))


