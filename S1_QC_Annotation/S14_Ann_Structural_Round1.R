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

SaveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"
list.files(SaveFolder)


#################################################################
seu0<-readRDS(file=sprintf("%s/seuThymoma.RDS",saveFolder))
seu0$Study<-dplyr::case_when(
  seu0$GSE %in% c('ThymusMay242023','Thymus_June2023_GEMonly','Thymus_Sep2023_GEMonly',
                 'Thymus_Dec2023_GEMonly','Thymus_Oct2023_GEMonly')~'NU',
  seu0$GSE %in% c('Xin')~'Xin',
  seu0$GSE %in% c('Yasumizu')~'Yasumizu',
  .default="unknown"
)
temp_Stru<-subset(seu0,Population2 %in% c("Structual"))

rawCount_temp_Stru0<-temp_Stru@assays$RNA@counts
har_Stru0<-CreateSeuratObject(counts = rawCount_temp_Stru0, project ="Structural", min.cells = 3) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE)
har_Stru0<-AddMetaData(har_Stru0,metadata=temp_Stru@meta.data[,c(5:ncol(temp_Stru@meta.data))])
har_Stru0<-har_Stru0 %>% RunHarmony("ID", plot_convergence = FALSE)
har_Stru0 <- har_Stru0 %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50, k.param = 50) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

################remove doublet and reclustering ###########################
temp_StruF<-subset(har_Stru0, RNA_snn_res.0.5 %in% '14',invert=T)
rawCount_temp_Stru<-temp_StruF@assays$RNA@counts
har_Stru<-CreateSeuratObject(counts = rawCount_temp_Stru, project ="Structural", min.cells = 3) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 6000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 40, verbose = FALSE)
har_Stru<-AddMetaData(har_Stru,metadata=temp_StruF@meta.data[,c(4:(ncol(temp_StruF@meta.data)-2))])
har_Stru<-har_Stru %>% RunHarmony("ID", plot_convergence = FALSE)
har_Stru <- har_Stru %>%
  RunUMAP(reduction = "harmony", dims = 1:40) %>%
  FindNeighbors(reduction = "harmony", dims = 1:40, k.param =100) %>%
  FindClusters(resolution = 1.2) %>%
  identity()

har_Stru$MainCell<-dplyr::case_when(
  har_Stru$RNA_snn_res.1.2 %in% c("6","16") ~ "Endothelial",
  har_Stru$RNA_snn_res.1.2 %in% c('7',"12",'21')~"Stromal",
  .default='Epithelial'
)

har_Stru$CellType1<-dplyr::case_when(
  har_Stru$RNA_snn_res.1.2 %in% c('4')~"cTEC", #
  har_Stru$RNA_snn_res.1.2 %in% c('20')~'Myoid', 
  har_Stru$RNA_snn_res.1.2 %in% c('11')~'Neuro',
  har_Stru$RNA_snn_res.1.2 %in% c('18')~'Ciliated',
  har_Stru$RNA_snn_res.1.2 %in% c("15")~"Ionocyte",
  har_Stru$RNA_snn_res.1.2 %in% c("10") ~ "Tuft", 
  har_Stru$RNA_snn_res.1.2 %in% c('1','5','8','13')~'mcTEC',
  har_Stru$RNA_snn_res.1.2 %in% c('0','2','3','19','9','17','14')~'mTEC',
  har_Stru$RNA_snn_res.1.2 %in% c("6","16")~"Endo", 
  har_Stru$RNA_snn_res.1.2 %in% c("12")~"VSMC",  
  har_Stru$RNA_snn_res.1.2 %in% c("7",'21')~"Fibro", 
  .default='unknown'
)

saveRDS(har_Stru,file=sprintf("%s/har_Stru.RDS",FigOut))


######################################################################################
######################################################################################
######################################################################################
harFun<-function(obj=obj,projN='Test',nf=3000,pcs=40,kP=60,res=1.2){
  rawCount_temp_In<-obj@assays$RNA@counts
  har_obj<-CreateSeuratObject(counts = rawCount_temp_In, project =projN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nf) %>% ##2500(38;60)---1742;171;349;49/159--3500 (45/174)
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

generateSubCl<-function(obj=seuST,proN="AdvF",nFeature=3000,npcs=35,res=0.9,kparam=50){
  temp_CT<-obj
  rawCount_temp<-temp_CT@assays$RNA@counts
  har_temp<-CreateSeuratObject(counts = rawCount_temp, project =proN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nFeature) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE)
  har_temp<-AddMetaData(har_temp,metadata=temp_CT@meta.data[,c(5:ncol(temp_CT@meta.data))])
  har_temp<-har_temp %>% RunHarmony("ID", plot_convergence = FALSE)
  har_temp <- har_temp %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs) %>%
    FindNeighbors(reduction = "harmony", dims = 1:npcs, k.param = kparam) %>%
    FindClusters(resolution = res) %>%
    identity()
  return(har_temp)
}

generateSubClFromHar<-function(obj=seuST,proN="AdvF",nFeature=3000,npcs=35,res=0.9,kparam=50){
  temp_CT<-obj
  rawCount_temp<-temp_CT@assays$RNA$counts
  har_temp<-CreateSeuratObject(counts = rawCount_temp, project =proN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nFeature) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE)
  har_temp<-AddMetaData(har_temp,metadata=temp_CT@meta.data[,c(5:ncol(temp_CT@meta.data))])
  har_temp<-har_temp %>% RunHarmony("ID", plot_convergence = FALSE)
  har_temp <- har_temp %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs) %>%
    FindNeighbors(reduction = "harmony", dims = 1:npcs, k.param = kparam) %>%
    FindClusters(resolution = res) %>%
    identity()
  return(har_temp)
}

######################################################################################
har_EP<-subset(har_Stru,CellType1 %in% c("Endo",'VSMC','Fibro'),invert=T)
har_EnSt<-subset(har_Stru,CellType1 %in% c("Endo",'VSMC','Fibro'),invert=F)

har_EP2<-generateSubCl(obj=har_EP,
                       proN="EP",
                       nFeature=5000,#
                       npcs=35,
                       res=0.1,
                       kparam=80)

Idents(har_EP2)<-"RNA_snn_res.0.1"
har_EP2<-FindSubCluster(har_EP2,cluster='2',graph.name="RNA_snn", subcluster.name='C2',resolution=0.1)

############################################################
har_EnSt2<-generateSubCl(obj=har_EnSt,
                         proN="EnSt",
                         nFeature=3500,#
                         npcs=50,
                         res=0.2,
                         kparam=50)

Idents(har_EnSt2)<-"RNA_snn_res.0.2"
har_EnSt2<-FindSubCluster(har_EnSt2,cluster='6',graph.name="RNA_snn", subcluster.name='C6',resolution=0.1)
Idents(har_EnSt2)<-"C6"

har_EnSt3<-generateSubCl(obj=subset(har_EnSt2,C6 %in% '6_1',invert=T),
                         proN="EnSt",
                         nFeature=1000,#
                         npcs=30,
                         res=0.2,
                         kparam=20)#

DefaultAssay(har_EnSt3)<-"RNA";Idents(har_EnSt3)<-"RNA_snn_res.0.2"
har_EnSt3<-FindSubCluster(har_EnSt3,cluster='0',graph.name="RNA_snn", subcluster.name='C0',resolution=0.15)
Idents(har_EnSt3)<-"C0"
DefaultAssay(har_EnSt3)<-"RNA";Idents(har_EnSt3)<-"C0"
har_EnSt3<-FindSubCluster(har_EnSt3,cluster='2',graph.name="RNA_snn", subcluster.name='C2',resolution=0.15)
Idents(har_EnSt3)<-"C2"

har_EnSt3$SubCT<-case_when(
  har_EnSt3$C2 %in% c('6')~"Lymphatic", 
  har_EnSt3$C2 %in% c('0_1')~'Vascular',
  har_EnSt3$C2 %in% c('0_2')~'Cap1', 
  har_EnSt3$C2 %in% c('0_0')~'Endo', 
  har_EnSt3$C2 %in% c('3')~'Activated.En', 
  har_EnSt3$C2 %in% c('5')~'Cap2',
  har_EnSt3$C2 %in% c('9')~'Doublet',
  har_EnSt3$C2 %in% c('2_1')~'Pericyte',
  har_EnSt3$C2 %in% c('2_0')~'VSMC',
  har_EnSt3$C2 %in% c('1','4','7','8')~'Fibroblast', 
  .default='unknown'
)
saveRDS(har_EnSt3,file=sprintf("%s/har_EnSt.RDS",SaveFolder))

EnSt_genes <- c("PECAM1","VWF",
                'MMRN1','CCL21', 
                "DKK2",'SEMA3G','GJA5',
                "IL7R",#Cap1
                "HPGD","IL1RL1",
                'ACKR1',"SELE",  
                "LUM","DCN",
                "COX4I2","HIGD1B", 
                "MYH11","MUSTN1"
)
har_EnSt3_clean<-subset(har_EnSt3,SubCT %in% "Doublet",invert=T)
har_EnSt3_clean$SubCT<-factor(har_EnSt3_clean$SubCT,
                              levels=rev(c("Endo","Lymphatic","Vascular","Cap1","Cap2", "Activated.En", 
                                           "Fibroblast", "Pericyte","VSMC")))

##########################################################################################
###################################     Figure S2A        ################################
##########################################################################################

Idents(har_EnSt3_clean)<-"SubCT";DefaultAssay(har_EnSt3_clean)<-"RNA"
dotEn<-DotPlot(object =har_EnSt3_clean, features = EnSt_genes,cols = c("lightgrey", "red"),scale.by='radius',col.max=1)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf(sprintf("%s/dotplot_EN.pdf",FigOut),width=7.87*1.5,height=7.87/1.5, family="ArialMT")
print(dotEn)
dev.off()


pdf(sprintf("%s/Dim_ENSt.pdf",FigOut),width=7.87*1.1,height=7.87, family="ArialMT")
print(DimPlot(har_EnSt3_clean,group.by='SubCT',label=T,repel=T))
dev.off()
tiff(sprintf("%s/Dim_ENSt.tiff",FigOut),width=7.87*1.1,height=7.87, res=600, units="in",compression='lzw')
print(DimPlot(har_EnSt3_clean,group.by='SubCT',label=F,repel=T,raster=F)+NoLegend())
dev.off()





############################################################
##       annotation of epithelial cells 
############################################################
DefaultAssay(har_EP2)<-"RNA";Idents(har_EP2)<-"RNA_snn_res.0.1"
har_EP2<-FindSubCluster(har_EP2,cluster='2',graph.name="RNA_snn",subcluster.name='C2',resolution=0.2,algorithm=3)
Idents(har_EP2)<-'C2'
har_EP2<-FindSubCluster(har_EP2,cluster='6',graph.name="RNA_snn",subcluster.name='C6',resolution=0.05,algorithm=3)
Idents(har_EP2)<-'C6'
har_EP2$SubCT<-case_when(
  har_EP2$C6 %in% c('3')~"cTEC", #
  har_EP2$C6 %in% c('9')~'Myoid', 
  har_EP2$C6 %in% c('6_0')~'Neuro',
  har_EP2$C6 %in% c('6_1')~'Ciliated',
  har_EP2$C6 %in% c("8")~"Ionocyte",
  har_EP2$C6 %in% c("5") ~ "Tuft",
  har_EP2$C6 %in% c('1')~'mcTEC',
  har_EP2$C6 %in% c('0')~'KRT15+.mTEC',
  har_EP2$C6 %in% c('4','7')~'SLPI+.mTEC',
  har_EP2$C6 %in% c('2_0')~'SLC39A2+.Tumor',
  har_EP2$C6 %in% c('2_1')~'CHGA+.Tumor',
  .default='unknown'
)

EP_genes <- c("EPCAM",'LY75','PSMB11',"CCL25", 
              "CCL21",  #mcTEC
              "KRT15","KRT14", 'MEG3',
              "SLPI","CALML5",'LCN2', 
              
              "MYH3","MYOG","ACTC1","MYLPF",
              "NEUROD1","GNG8","NHLH1","STMN2",
              "GFI1","LHX3","FOXJ1","ATOH1", 
              "FOXI1", "ASCL3", "CLCNKB", 
              "GNB3", "GNAT3", "PLCB2", "OVOL3",  
              
              'SLC39A2','UGT2A1',"CEL", 
              "CHGA","TMEM132A","RAMP1"
)
har_EP2$SubCT<-factor(har_EP2$SubCT,
                      levels=rev(c("cTEC","mcTEC", "KRT15+.mTEC", 'SLPI+.mTEC',
                                   "Myoid","Neuro","Ciliated","Ionocyte","Tuft",
                                   "SLC39A2+.Tumor","CHGA+.Tumor")))

##########################################################################################
###################################     Figure S2A         ################################
##########################################################################################
Idents(har_EP2)<-"SubCT"
DefaultAssay(har_EP2)<-"RNA"
dot<-DotPlot(object =har_EP2, features = EP_genes,cols = c("lightgrey", "red"),scale.by='radius',col.max=1)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf(sprintf("%s/dotplot_EP.pdf",FigOut),width=7.87*1.5,height=7.87/1.5, family="ArialMT")
print(dot)
dev.off()


pdf(sprintf("%s/Dim_EP.pdf",FigOut),width=7.87*1.1,height=7.87, family="ArialMT")
print(DimPlot(har_EP2,group.by='SubCT',label=T,repel=T))
dev.off()
tiff(sprintf("%s/Dim_EP.tiff",FigOut),width=7.87*1.1,height=7.87, res=600, units="in",compression='lzw')
print(DimPlot(har_EP2,group.by='SubCT',label=F,repel=T)+NoLegend())
dev.off()

saveRDS(har_EP2,file=sprintf("%s/Thymoma.har_EP2_forFig.RDS",SaveFolder))



