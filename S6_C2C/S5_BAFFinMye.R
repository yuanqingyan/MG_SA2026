
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
har_IM<-readRDS(file=sprintf("%s/har_IM_annotated.RDS",saveFolder))
DefaultAssay(har_IM)<-'RNA'
Idents(har_IM)<-"CellType_Round1"


IM_MG<-subset(har_IM, Disease2 %in% c("MG",'NoMG'))

FeaturePlot(IM_MG, feature=c('TNFSF13B'),col=c('grey','red'),split.by="Disease2",raster=F)

Mye<-subset(IM_MG, CellType_Round1 %in% c("DC",'pDC','Macrophage','Neutroplils'))
Mye$Disease2<-factor(Mye$Disease2,levels=c("MG",'NoMG'))
(Baff_Mye<-VlnPlot(Mye,feature=c('TNFSF13B'),raster=F,pt.size=0.05,split.by="Disease2"))

##### Fig. S8D ########
pdf(sprintf("%s/BAFF_Mye.pdf",FigOut),width=7.87*1.2,height=7.87/1.2, family="ArialMT")
print(Baff_Mye)
dev.off()

tiff(sprintf("%s/BAFF_Mye.tiff",FigOut),width=7.87*1.2,height=7.87/1.2, family="ArialMT",res=600, units="in",compression='lzw')
print(Baff_Mye)
dev.off()

