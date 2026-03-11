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
seuB_NU<-readRDS(file=sprintf("%s/seuB_NU_readVDJ_fromDJVDJ.RDS",saveFolder))
FeaturePlot(seuB_NU,features=c("n_chains"),col=c('grey','red'))
seuB_NU |> slot("meta.data") |> filter(n_chains > 1) |> head(5)
table(seuB_NU$n_chains)

# To only include clonotypes with exactly 2 chains----remove clonetype, but not cells
so_vdj <- seuB_NU |> filter_vdj(n_chains == 2);dim(so_vdj)
table(so_vdj$Bsub);dim(so_vdj);table(so_vdj$n_chains)

seuClean <- so_vdj;unique(seuClean$Disease_Fig)

seuClean$Disease2<-dplyr::case_when(
  seuClean$Disease_Fig %in% c("Thymoma_MG","Thymus_MG") ~"MG",
  seuClean$Disease_Fig %in% c("Thymoma_no_MG","Thymus_no_MG")~'NoMG',
  .default="Others"
)
seuClean$Disease2<-factor(seuClean$Disease2,levels=rev(c("MG","NoMG","Others")))

keepCells<-(seuClean@meta.data[!is.na(seuClean@meta.data$clonotype_id_F),'mergID'])
seuClean<-subset(seuClean, mergID %in% keepCells)

seuMG<-subset(seuClean, Disease2 %in% c("MG",'NoMG'))

#############################################################################
########             Clonotype Frequencies                      #############
#############################################################################
###Calculating clonotype frequencies
so_vdj_F <- seuClean |> calc_frequency(data_col = "clonotype_id_F")
so_vdj2 <- seuClean |>calc_frequency(data_col = "clonotype_id_F", cluster_col = "Disease2")
so_vdj2@meta.data[so_vdj2$clonotype_id_F_pct>5,]
(p1_MG<-seuMG |>  plot_clone_frequency(data_col = "clonotype_id_F", plot_colors = "#3182bd",units='frequency'))

cellPer<-(seuMG |> plot_clone_frequency(data_col = "clonotype_id_F", 
                               cluster_col  = "Disease2", 
                               panel_scales = "free_x",
                               plot_lvls=c("MG",'NoMG'),
                               units='frequency',#'percent'
                               alpha=1,
                               plot_colors=c('tomato','light blue')))

#### figure 2D #### 
pdf(sprintf("%s/barplot_clonotypeFreq_byMG.pdf",FigOut),width=6, height=6/1.5, family="ArialMT")
print(cellPer)
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


#### figure 2E #### 
pdf(sprintf("%s/Thymus_DiversityAnalysis.pdf",FigOut),width=9,height=6, family="ArialMT")
print(diversity)
dev.off()
#####################################################





