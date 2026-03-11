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


All.B<-readRDS(file=sprintf("%s/temp/AllB_withCluster.RDS",saveFolder))
metaPred<-read.csv("./Thymus/h5ad/ThymomaB_May2025_Bcells_predictFromCellTypise.csv",row.names=1)
sum(row.names(metaPred)==row.names(All.B@meta.data))

All.B$predicted_labels<-metaPred$ThymomaLabel
All.B$conf_score<-metaPred$ThymomaConfScore
All.B$Majority_voting<-metaPred$ThymomaMajority_voting

All.B$AtlasLabel<-metaPred$AtlasLabel
All.B$AtlasMajority_voting<-metaPred$AtlasMajority_voting
All.B$AtlasConfScore<-metaPred$AtlasConfScore

All.B$ImHLabel<-metaPred$ImHLabel
All.B$ImHMajority_voting<-metaPred$ImHMajority_voting
All.B$ImHConfScore<-metaPred$ImHConfScore

All.B$ImLLabel<-metaPred$ImLLabel
All.B$ImLMajority_voting<-metaPred$ImLMajority_voting
All.B$ImLConfScore<-metaPred$ImLConfScore

sum(All.B$conf_score>=0.7);sum(All.B$conf_score<0.7)
All.B$predictedLab2<-ifelse(All.B$conf_score>=0.7,All.B$predicted_labels,'Unknown')
All.B$preConf<-ifelse(All.B$conf_score>=0.7,'Good','Poor')
saveRDS(All.B,file=sprintf("%s/temp/AllB_withPrediction.RDS",saveFolder))

##########################################################################################
All.B<-readRDS(file=sprintf("%s/temp/AllB_withPrediction.RDS",saveFolder))
DefaultAssay(All.B)<-"RNA";Idents(All.B)<-'C4'

All.B$Bsub<-dplyr::case_when(
  All.B$C4 %in% c("7_0")~"Centroblast_GC", 
  All.B$C4 %in% c('7_1')~"Centrocyte_GC", 
  All.B$C4 %in% c("11_0")~"Mature_GC", 
  All.B$C4 %in% c("11_1")~"Early_GC",
  All.B$C4 %in% c("4_0","4_1")~"Plasma",
  All.B$C4 %in% c("9")~"Naiv_Cir_1",
  All.B$C4 %in% c("8")~"Naiv_Cir_2",
  All.B$C4 %in% c("13")~"Naiv_Cir_3",
  All.B$C4 %in% c("12")~"Naiv_Cir_4",
  All.B$C4 %in% c("6")~"Naiv_Cir_5",
  All.B$C4 %in% c("1",'0_2')~"NaivB",
  All.B$C4 %in% c("2_0",'2_3')~"MemB_Cir",
  All.B$C4 %in% c("2_1",'2_2')~"MemB",
  All.B$C4 %in% c("5")~"Un_switch",
  All.B$C4 %in% c("3")~"Switch",
  All.B$C4 %in% c("10")~"MemB.ISG",
  All.B$C4 %in% c("0_0")~"MemB.Rib+Mit+",
  All.B$C4 %in% c("0_1")~"MemB.Rib+",
  All.B$C4 %in% c("0_3")~"MemB.Mit+",
  All.B$C4 %in% c("2_4")~"MemB.EBI3+",
  .default='unknown');table(All.B$Bsub)


All.B$BigC<-dplyr::case_when(
  All.B$C4 %in% c("7_0",'7_1','11_0','11_1')~"GC.B",
  All.B$C4 %in% c("4_1",'4_0')~"ASC",
  .default="no.GC.B");table(All.B$BigC)
##############################################################################
saveRDS(All.B, file=sprintf("%s/B_withSubtype_Annotated_May2025.RDS",saveFolder))
##############################################################################


