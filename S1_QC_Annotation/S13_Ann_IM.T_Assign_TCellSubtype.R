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

har_T<-readRDS(file=sprintf("%s/T.Pop.RDS",FigOut))
metaCellTypist<-read.csv("./Thymus/h5ad/T/outMeta/metadata.csv")
sum(metaCellTypist$X==row.names(har_T@meta.data))
metaToAddIn<-metaCellTypist[,(which(colnames(metaCellTypist)=="ThymomaLabel")):ncol(metaCellTypist)]
har_T<-AddMetaData(object = har_T,metadata = metaToAddIn)

##########
harFun<-function(obj=obj,projN='Test',nf=3000,pcs=40,kP=60,res=1.2){
  rawCount_temp_In<-obj@assays$RNA@counts
  har_obj<-CreateSeuratObject(counts = rawCount_temp_In, project =projN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nf) %>% 
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
#######
harT2<-harFun(obj=subset(har_T,CellType1 %in% 'Mast',invert=T),
              projN='Tcells',
              nf=3000,pcs=50,kP=60,res=1.1)

harT2$Sub1<-dplyr::case_when(
  harT2$RNA_snn_res.1.1 %in% c("0",'21','23','26','24')~'DPT',
  harT2$RNA_snn_res.1.1 %in% c("6",'7','31','29')~'DNT',
  harT2$RNA_snn_res.1.1 %in% c("8",'2','30','10','13','25')~'Prolf.T',
  harT2$RNA_snn_res.1.1 %in% c('20')~'NK',
  harT2$RNA_snn_res.1.1 %in% c('22')~'ILC', ##lack of CD3D
  .default='OtherT')
saveRDS(harT2, file=sprintf("%s/Incase_harT2_1.RDS",FigOut))

otherTID<-row.names(subset(harT2,Sub1 %in% 'OtherT')@meta.data)


###################################################################################
###################################################################################
###################################################################################
har_T$cellID<-row.names(har_T@meta.data)
harT3<-harFun(obj=subset(har_T,cellID %in% otherTID),
              projN='OtherT',
              nf=3000,pcs=50,kP=60,res=1.2)

harT3$Sub1<-dplyr::case_when(
  harT3$RNA_snn_res.1.2 %in% c("22",'27','20','3')~'DPT',
  harT3$RNA_snn_res.1.2 %in% c("23")~'Prolf.T',
  .default='OtherT')
saveRDS(harT3, file=sprintf("%s/Incase_harT3_1.RDS",FigOut))

harT3<-readRDS(file=sprintf("%s/Incase_harT3_1.RDS",FigOut))
DPID<-row.names(subset(harT3,Sub1 %in% c("DPT"))@meta.data)
ProfID<-row.names(subset(harT3,Sub1 %in% c("Prolf.T"))@meta.data)

####################################################
harT2$CellID<-row.names(harT2@meta.data)
harT2$Sub1<-ifelse(harT2$CellID %in% DPID, "DPT",harT2$Sub1)
harT2$Sub1<-ifelse(harT2$CellID %in% ProfID, "Prolf.T",harT2$Sub1)

saveRDS(harT2, file=sprintf("%s/Incase_harT2_2.RDS",FigOut))
####################################################


####################################################
otherTID2<-row.names(subset(harT2,Sub1 %in% 'OtherT')@meta.data)
harT4<-harFun(obj=subset(har_T,cellID %in% otherTID2),
              projN='OtherT',
              nf=3000,pcs=50,kP=60,res=0.9)
saveRDS(harT4, file=sprintf("%s/Incase_harT4_1.RDS",FigOut))

DefaultAssay(harT4)<-"RNA"
Idents(harT4)<-"RNA_snn_res.0.9"
Marker_T4<-FindAllMarkers(harT4, slot='data', test.use='wilcox', only.pos=T)
Marker_T4 %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% as.data.frame()
Marker_T4Sort<-Marker_T4[order(-Marker_T4$avg_log2FC),]

Idents(harT4)<-"RNA_snn_res.0.9"
harT4<-FindSubCluster(harT4, cluster='3',graph.name='RNA_snn',subcluster.name = "C3",resolution = 0.3,algorithm = 1)
Idents(harT4)<-"C3"
harT4<-FindSubCluster(harT4, cluster='4',graph.name='RNA_snn',subcluster.name = "C4",resolution = 0.3,algorithm = 1)
Idents(harT4)<-"C4"
harT4<-FindSubCluster(harT4, cluster='2',graph.name='RNA_snn',subcluster.name = "C2",resolution = 0.3,algorithm = 1)
Idents(harT4)<-"C2"
harT4<-FindSubCluster(harT4, cluster='1',graph.name='RNA_snn',subcluster.name = "C1",resolution = 0.3,algorithm = 1)
DimPlot(harT4, group.by='C1',label=T, repel=T,label.size=4,raster=T)+NoLegend()
Idents(harT4)<-"C1"
harT4<-FindSubCluster(harT4, cluster='9',graph.name='RNA_snn',
                      subcluster.name = "C9",resolution = 0.2,algorithm = 1);table(harT4$C9)
Idents(harT4)<-"C9"
harT4<-FindSubCluster(harT4, cluster='10',graph.name='RNA_snn',
                      subcluster.name = "C10",resolution = 0.7,algorithm = 1)
Idents(harT4)<-"C10"
harT4<-FindSubCluster(harT4, cluster='12',graph.name='RNA_snn',
                      subcluster.name = "C12",resolution = 0.2,algorithm = 1)
Idents(harT4)<-'C12'
harT4<-FindSubCluster(harT4, cluster='6',graph.name='RNA_snn',
                      subcluster.name = "C6",resolution = 0.1,algorithm = 1)
Idents(harT4)<-"C6"
#############################################
DefaultAssay(harT4)<-"RNA";Idents(harT4)<-"C6"
Marker_FT4<-FindAllMarkers(harT4, slot='data', test.use='wilcox', only.pos=T)
Marker_FT4 %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% as.data.frame()
Marker_FT4Sort<-Marker_FT4[order(-Marker_FT4$avg_log2FC),]
write.csv(Marker_FT4Sort, file=sprintf("%s/Marker_T4Sort_FinalMarker.csv",FigOut))
#####################

harT4$Sub_Tcell<-dplyr::case_when(
  harT4$C6 %in% c("0")~'CD8a/b(entry)', 
  harT4$C6 %in% c("5")~'Treg', 
  harT4$C6 %in% c("14")~'Tfh', 
  harT4$C6 %in% c("13")~'gdT', 
  harT4$C6 %in% c("3_0")~'CD8T.Naive', 
  harT4$C6 %in% c("3_1",'3_2','15','16','6_1','12_0','12_1')~'CD4T.Naive',
  harT4$C6 %in% c("4_0",'4_2')~'CD4T.Tcm(Th0)', 
  harT4$C6 %in% c("2_0")~'CD4T.Tcm(Th17)', 
  harT4$C6 %in% c("2_1")~'CD4T.Tcm(Th2)', 
  harT4$C6 %in% c("2_2")~'CD4T.Tem(Th1/Th17)', 
  harT4$C6 %in% c("1_3")~'CD4T.Tem(Th1)', 
  harT4$C6 %in% c("9_2",'9_3')~'CD4T.Temra', 
  harT4$C6 %in% c('10_0','10_1','10_3','10_4')~'T(agonist)', 
  harT4$C6 %in% c("11")~'CRTAM+.gdT', 
  harT4$C6 %in% c("10_2",'8')~'CD8a/a', 
  harT4$C6 %in% c("17",'9_1','9_0')~'CD8.Temra', 
  harT4$C6 %in% c("1_2")~'CD8.Trm', 
  harT4$C6 %in% c("1_0",'1_1','12_2')~'CD8.Tem', 
  harT4$C6 %in% c("7",'6_0')~'DPT', #
  harT4$C6 %in% c("4_1")~'CD8.Tcm', 
  .default='Unknown');table(harT4$Sub_Tcell)


####################################################
#############  for meta data      ##################
####################################################
harT4$Info_SampleType<-ifelse(harT4$ID %in% 'MG03_PT',"Thymoma_MG",harT4$Info_SampleType);
harT4$organ<-ifelse(harT4$ID %in% c("MG03_PB","MG03_PT","MG22_PL"),'Blood','Tissue')
saveRDS(harT4, file=sprintf("%s/annotated_harT4__June2025.RDS",saveFolder))
saveRDS(harT4, file=sprintf("%s/harT4_withSubTcell.RDS",FigOut))
saveRDS(harT4@meta.data, file=sprintf("%s/harT4_withSubTcell_meta.RDS",FigOut))

####################################################
T4_meta<-readRDS(file=sprintf("%s/harT4_withSubTcell_meta.RDS",FigOut))


harT2<-readRDS(file=sprintf("%s/Incase_harT2_2.RDS",FigOut))
T4_meta<-T4_meta[,c("cellID",'Sub_Tcell','C6')]
colnames(T4_meta)[1]<-'CellID'
temp_T2<-harT2@meta.data
temp_merg<-dplyr::left_join(temp_T2, T4_meta,by='CellID');head(temp_merg)
temp_merg$TSub_F<-ifelse(temp_merg$Sub1 %in% c("OtherT"),temp_merg$Sub_Tcell,temp_merg$Sub1)
row.names(temp_merg)<-temp_merg$CellID


harT2$TSub_F<-temp_merg$TSub_F
harT2$C6_fromT4<-temp_merg$C6
harT2$TSub_F<-gsub('CD4T','CD4',harT2$TSub_F)


####################################################
#############  for meta data      ##################
####################################################
harT2$Info_SampleType<-ifelse(harT2$ID %in% 'MG03_PT',"Thymoma_MG",harT2$Info_SampleType);
harT2$organ<-ifelse(harT2$ID %in% c("MG03_PB","MG03_PT","MG22_PL"),'Blood','Tissue')

harT2$Disease2<-dplyr::case_when(
  harT2$Disease_Fig %in% "Others"~'Others',
  harT2$Disease_Fig %in% c("Thymoma_MG","Thymus_MG") ~"MG",
  harT2$Disease_Fig %in% c("Thymoma_no_MG","Thymus_no_MG")~'NoMG',
  .default="unknnown"
);
####################################################
saveRDS(harT2, file=sprintf("%s/annotated_T_NK_ILC__June2025.RDS",saveFolder))
####################################################
