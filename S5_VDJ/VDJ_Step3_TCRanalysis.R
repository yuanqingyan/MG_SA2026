library(djvdj)
library(Seurat)
library(dplyr)
library(ggplot2)
##############
saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"
list.files(saveFolder)

######################################################################
######################################################################
######################################################################
seuT_NU<-readRDS(file=sprintf("%s/seuT_NU_readVDJ_fromDJVDJ.RDS",saveFolder));head(seuT_NU);dim(seuT_NU)##88,649cells

so_vdj <- seuT_NU |> filter_vdj(n_chains == 2)
seuClean <- so_vdj

seuClean$Disease2<-dplyr::case_when(
  seuClean$Disease_Fig %in% c("Thymoma_MG","Thymus_MG") ~"MG",
  seuClean$Disease_Fig %in% c("Thymoma_no_MG","Thymus_no_MG")~'NoMG',
  .default="Others"
)
seuClean$Disease2<-factor(seuClean$Disease2,levels=rev(c("MG","NoMG","Others")))
keepCells<-(seuClean@meta.data[!is.na(seuClean@meta.data$clonotype_id_F),'mergID'])
seuClean<-subset(seuClean, mergID %in% keepCells)
seuClean$Disease2<-as.character(seuClean$Disease2)
seuClean<-subset(seuClean, Disease2 %in% c("MG",'NoMG'))

seuMG<-subset(seuClean, Disease2 %in% c("MG",'NoMG'))
seuMG$TSub_F<-as.character(seuMG$TSub_F)
uniTsubF<-unique(seuMG$TSub_F)
seuMG<-subset(seuMG, TSub_F %in% uniTsubF[!uniTsubF %in% c("NK",'ILC')])
seuMG_meta<-seuMG@meta.data[,which(colnames(seuMG@meta.data)=='Disease2'):ncol(seuMG@meta.data)]

#############################################################################
########             Clonotype Frequencies                      #############
#############################################################################
###Calculating clonotype frequencies
so_vdj_F <- seuMG |> calc_frequency(data_col = "clonotype_id_F")
so_vdj2 <- seuMG |>calc_frequency(data_col = "clonotype_id_F", cluster_col = "Disease2")
so_vdj2@meta.data[so_vdj2$clonotype_id_F_pct>0.5,]
###################################################

(p1_MG<- (seuMG |>  plot_clone_frequency(data_col = "clonotype_id_F", 
                                       plot_colors = "#3182bd",
                                       cluster_col  = "orig.ident",
                                       panel_scales = "free_x",
                                       units='frequency')))

########## figure S7A ##########
pdf(sprintf("%s/barplot_clonotypeFreq_Tcell.pdf",FigOut),width=6,height=6, family="ArialMT")
print(p1_MG)
dev.off()


(cellPer<-(seuMG |> plot_clone_frequency(data_col = "clonotype_id_F", 
                               cluster_col  = "Disease2", 
                               panel_scales = "free_x",
                               plot_lvls=c("MG",'NoMG'),
                               units='percent',#'percent'
                               alpha=1,
                               plot_colors=c('tomato','light blue'))))

########## figure S7B ##########
pdf(sprintf("%s/barplot_clonotypePercentage_T_byMG.pdf",FigOut),width=6, height=6/1.5, family="ArialMT")
print(cellPer)
dev.off()

(cellty<-(seuMG |> plot_clone_frequency(data_col = "clonotype_id_F",
                                      cluster_col  = "TSub_F", 
                                      units='frequency',#frequency
                                      panel_scales = "free_x")))

########## figure S7C ##########
pdf(sprintf("%s/barplot_clonotypePercentage_T_byCellType.pdf",FigOut),width=12, height=12, family="ArialMT")
print(cellty)
dev.off()


seuMG$dis2PT<-sprintf("%s_%s",seuMG$Disease2,seuMG$Info_Patient)
#############################################################################
########             Repertoire Diversity                       #############
#############################################################################
#Calculating diversity
seuMG2 <- seuMG |>
  calc_diversity(
    data_col    = "clonotype_id_F",
    cluster_col = "Disease2",
    method      = abdiv::shannon,
    downsample  = TRUE,
    n_boots     = 500)

div_fns <- list(
  "simpson" = abdiv::simpson,
  "shannon" = abdiv::shannon,
  "pielou evenness" = abdiv::pielou_e)

seuMG2 |>
  plot_diversity(
    data_col    = "clonotype_id_F",
    cluster_col = "Disease2",
    method      = div_fns,
    plot_lvls=c("MG",'NoMG'),
    plot_colors=c('tomato','light blue'),
    alpha=1
  )+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


seuMG2$orig.ident<-seuMG2$Disease2
(diversity<-seuMG2 |>
    plot_diversity(
      data_col    = "clonotype_id_F",
      cluster_col = "patient",
      group_col   = "orig.ident",
      method      = div_fns,
      downsample  = FALSE,
      n_boots     = 500,
      plot_lvls=c("MG",'NoMG'),
      plot_colors=c('tomato','light blue'),
      alpha=1
    ))


############ Figure 3C #########
pdf(sprintf("%s/Thymus_DiversityAnalysis__Tcells.pdf",FigOut),width=9,height=6, family="ArialMT")
print(diversity)
dev.off()


seuMG2$PTCellCtype<-paste(seuMG2$patient, seuMG2$TSub_F,sep="_")
CTlevel<-c("DPT",'DNT','Prolf.T','CD8a/b(entry)','CD8a/a',
           'CD8T.Naive','CD8.Tcm','CD8.Tem','CD8.Trm','CD8.Temra',
           'T(agonist)',
           'CD4.Naive','CD4.Tcm(Th0)','CD4.Tcm(Th2)','CD4.Tcm(Th17)',
           'CD4.Tem(Th1/Th17)','CD4.Tem(Th1)','CD4.Temra',
           'Tfh','Treg',
           'gdT','CRTAM+.gdT','NK','ILC');length(CTlevel)
celltypeColor<-c("#11FFFF", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                 "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
                 "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                 "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#6A3D9A")
colorCT<-data.frame(CTlevel=CTlevel,celltypeColor=celltypeColor);head(colorCT)
colorCTin<-colorCT[colorCT$CTlevel %in% seuMG2$TSub_F,]

seuMG2$orig.ident<-seuMG2$TSub_F
seuMG2$orig.ident<-factor(seuMG2$orig.ident,levels=colorCTin$CTlevel)




(diversity2<-seuMG2 |>
    plot_diversity(
      data_col    = "clonotype_id_F",
      cluster_col = "PTCellCtype",
      group_col   = "orig.ident",
      method      = div_fns,
      downsample  = FALSE,
      n_boots     = 500,
      plot_lvls=colorCTin$CTlevel,
      plot_colors=colorCTin$celltypeColor,
      alpha=1
    ))

###### Figure 3D ##############
pdf(sprintf("%s/Thymus_DiversityAnalysis__TsubCT.pdf",FigOut),width=18,height=6, family="ArialMT")
print(diversity2)
dev.off()






