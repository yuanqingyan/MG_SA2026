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

harT2<-readRDS(file=sprintf("%s/annotated_T_NK_ILC__June2025.RDS",saveFolder))
CTlevel<-c("DPT",'DNT','Prolf.T','CD8a/b(entry)','CD8a/a',
           'CD8T.Naive','CD8.Tcm','CD8.Tem','CD8.Trm','CD8.Temra',
           'T(agonist)',
           'CD4.Naive','CD4.Tcm(Th0)','CD4.Tcm(Th2)','CD4.Tcm(Th17)',
           'CD4.Tem(Th1/Th17)','CD4.Tem(Th1)','CD4.Temra',
           'Tfh','Treg',
           'gdT','CRTAM+.gdT','NK','ILC')
celltypeColor<-c("#11FFFF", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                 "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
                 "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                 "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#6A3D9A")
harT2$TSub_F<-factor(harT2$TSub_F,levels=CTlevel)



geneMarker<-c('CD4','CD8A','CD8B','MKI67','TOP2A',
              "TOX2",'SATB1','CCR9', 
              'GNG4','PDCD1', 
              'CCR7', 
              'CD44',
              'EOMES','GZMK',
              'RGS1','CCL4L2',#
              'FCGR3A','NKG7',
              'NR4A1','NFKBID',
              'SELL',
              'FOSB','FOS',
              'GATA3','CCR4',
              'RORC','TNFRSF4',
              'KLRG1',
              'GZMH','CST7',
              'TBX21','CX3CR1',
              'CXCR5','ICOS',
              'CTLA4','FOXP3','IL2RA',
              'TRDC','TRGC1',
              'IKZF2','KLRD1','ITGAD',
              'TRAC',
              'CD3D','AREG','TLE1','IL4I1'
)

Idents(harT2)<-"TSub_F"
(dotT<-DotPlot(harT2,
               features=geneMarker,
               cols = c("lightgrey", "red"),
               scale=T)+
    theme(axis.text.x = element_text(angle = 90, hjust=1)))#+ coord_flip()


############ Fig. S6C ##########
pdf(sprintf("%s/Dotplot_AllT.pdf",FigOut),width=7.87*2.2,height=7.87/1.1, family="ArialMT")
print(dotT)
dev.off()


############ Fig. S6A ##########
pdf(sprintf("%s/T_sub.pdf",FigOut),width=7.87*1.5,height=7.87, family="ArialMT")
print(DimPlot(harT2,group.by='TSub_F',label=T,repel=T,cols=celltypeColor,raster=F))
dev.off()
tiff(sprintf("%s/T_sub.tiff",FigOut),width=7.87*1.05,height=7.87, family="ArialMT",res=600, units="in",compression='lzw')
print(DimPlot(harT2,group.by='TSub_F',label=F,raster=F,cols=celltypeColor) & NoLegend() & NoAxes ())
dev.off()


MG_T2<-subset(harT2, Disease2 %in% c("MG",'NoMG') & 
                organ %in% c("Tissue") & 
                Study %in% 'NU')

MG_T2$TSub_F<-factor(MG_T2$TSub_F,levels=CTlevel)
dimT<-DimPlot(MG_T2, group.by="TSub_F",split.by="Disease2",
        label=T,repel=T,raster=F,
        cols=celltypeColor)+NoLegend()

###### Fig. 3A #########
pdf(sprintf("%s/T_sub_MGnoMG.pdf",FigOut),width=7.87*1.5,height=7.87, family="ArialMT")
print(dimT)
dev.off()

MG_T2$ptDis2<-paste(MG_T2$patient,MG_T2$Disease2,sep='.')
(tab_MG_T2<-table(MG_T2$TSub_F,MG_T2$ptDis2))
propMG<-as.data.frame(prop.table(tab_MG_T2,margin=1))
propMG$Disease<-sapply(strsplit(as.character(propMG$Var2),split="\\."),function(x) x[2])

MG_T2$tempU<-sprintf("%s__%s__%s",MG_T2$patient,MG_T2$Disease2,MG_T2$Disease_Fig)
MG_T2<-subset(MG_T2, TSub_F %in% c("NK",'ILC'),invert=T)
MG_T2$TSub_F<-as.character(MG_T2$TSub_F)
saveRDS(MG_T2, file=sprintf("%s/MG_T2.forProportionalAnalysis.RDS", FigOut))

#########################
MG_T2<-readRDS(file=sprintf("%s/MG_T2.forProportionalAnalysis.RDS", FigOut))
(T_tab<-(table(MG_T2$tempU,MG_T2$TSub_F)))
(colSu<-colSums(T_tab))
(roSum<-rowSums(T_tab))
(T_prop<-prop.table(T_tab,margin=1));sum(T_prop[1,])

plotDat<-t(T_prop)
dist.my <- function(x) as.dist(1-cor(t(x)))
hclust.my <- function(x) hclust(x, method="complete")

pPlot<-data.frame(ID=colnames(plotDat))
pPlot$PT<-sapply(strsplit(colnames(plotDat),split="\\__"),function(x) x[1])
pPlot$Disease<-sapply(strsplit(colnames(plotDat),split="\\__"),function(x) x[2])
pPlot$Type<-sapply(strsplit(sapply(strsplit(colnames(plotDat),split="\\__"),function(x) x[3]),split="\\_"),function(y) y[1])

pPlot<-pPlot[order(pPlot$Disease, pPlot$Type),];pPlot
proScale<-t(scale(t(plotDat[,pPlot$ID])));head(proScale)

library(ComplexHeatmap)
ann_colors = list(
  Disease=c(MG="#ff6347",NoMG="#add8e6"),
  Type=c(Thymus='orange',Thymoma='lightblue'))

top_annotation <- HeatmapAnnotation(
  Disease = factor(pPlot$Disease),
  Type = factor(pPlot$Type),
  col = ann_colors
)

mini<-min(proScale);max<-max(proScale)
negLen<-abs(mini)*100
posLen<-abs(max)*100
mycols<-c(colorRampPalette(colors = c("darkblue","blue","white"))(negLen),
          colorRampPalette(colors = c("white","red","darkred"))(posLen))

heatProp<-ComplexHeatmap::Heatmap(
  proScale,
  col = mycols,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  top_annotation= top_annotation,  # 
  column_names_gp = gpar(fontsize = 10),  # 
  show_column_names = TRUE,
  row_order = NULL,
  column_order = NULL
)

############ Fig. 3B ##########
pdf(sprintf("%s/hc_prop_T.pdf",FigOut),width=10,height=10)
print(heatProp)
dev.off()


##### dendrogram over the columns ######
heatProp_2<-ComplexHeatmap::Heatmap(
  proScale,
  col = mycols,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation= top_annotation,  # 
  column_names_gp = gpar(fontsize = 10),  
  show_column_names = TRUE,
  row_order = NULL,
  column_order = NULL
)

############ Fig. S6D ##########
pdf(sprintf("%s/hc_prop_T_2.pdf",FigOut),width=10,height=10)
print(heatProp_2)
dev.off()

#####################################################################
################# stacked barplot    ################################
#####################################################################
long_T_prop <- data.frame(T_prop);head(long_T_prop)
long_T_prop$Var1<-sapply(strsplit(as.character(long_T_prop$Var1),split="\\__"),function(x) sprintf("%s__%s",x[1],x[3]))
long_T_prop$Var1<-factor(long_T_prop$Var1,
                         levels=c("PT12__Thymoma_MG","PT13__Thymoma_MG","PT9__Thymoma_MG",
                                  "PT11__Thymus_MG","PT12__Thymus_MG","PT15__Thymus_MG","PT2__Thymus_MG","PT5__Thymus_MG", 
                                  "PT6__Thymus_MG", "PT7__Thymus_MG","PT9__Thymus_MG", 
                                  
                                  "PT1__Thymoma_no_MG", "PT10__Thymoma_no_MG", "PT4__Thymoma_no_MG",
                                  "PT16__Thymus_no_MG", "PT3__Thymus_no_MG", "PT4__Thymus_no_MG", "PT8__Thymus_no_MG"))

ctOrder<-c("DPT",'CD8a/a','Prolf.T', 'gdT', 'DNT','CD8.Trm','CD8.Tem','CD4.Temra','CD4.Tem(Th1/Th17)',#
           'CRTAM+.gdT','CD8.Temra','CD4.Tcm(Th17)','CD4.Tcm(Th2)','CD4.Tcm(Th0)','Treg',
           'Tfh','CD8.Tcm','CD4.Tem(Th1)','CD8T.Naive','CD4.Naive','T(agonist)','CD8a/b(entry)')
long_T_prop$Var2<-factor(long_T_prop$Var2,
                         levels=ctOrder)

barplt<-ggplot(long_T_prop, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent)+
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  labs(fill = "Cell Type", x = "", y = "Percentage (%)")

######### Fig. S6E ########
pdf(sprintf("%s/T_composition_barplot.pdf",FigOut),width=8,height=6)
print(barplt)
dev.off()


##########################################################################################################################################
##############################################  miloR to test Differential abundance  ####################################################
##########################################################################################################################################
library(rlang)
library(htmltools)
library(ape)
library(stringr)
library(miloR)
library(scater)
library(scran)
library(patchwork)

FigOut="./Thymus/FigOut"

MG_T2<-readRDS(file=sprintf("%s/MG_T2.forProportionalAnalysis.RDS", FigOut))
CTlevel<-c("DPT",'DNT','Prolf.T','CD8a/b(entry)','CD8a/a',
           'CD8T.Naive','CD8.Tcm','CD8.Tem','CD8.Trm','CD8.Temra',
           'T(agonist)',
           'CD4.Naive','CD4.Tcm(Th0)','CD4.Tcm(Th2)','CD4.Tcm(Th17)',
           'CD4.Tem(Th1/Th17)','CD4.Tem(Th1)','CD4.Temra',
           'Tfh','Treg',
           'gdT','CRTAM+.gdT','NK','ILC');length(CTlevel)
allCTinheat<-CTlevel[CTlevel %in% unique(MG_T2$TSub_F)]
seuIn<-MG_T2
seuIn$Disease2<-as.character(seuIn$Disease2)
seuIn$TSub_F<-as.character(seuIn$TSub_F)
##### downsample DNT; CD4.Naive, CD8a/b(entry), DPT and prolf.T to 4k cells
ncellKeep=4000
metaSeu<-seuIn@meta.data
CT_largN<-c("CD8a/b(entry)","DNT","DPT","Prolf.T",'CD4.Naive')
otherCellID<-row.names(metaSeu)[!metaSeu$TSub_F %in% CT_largN]

set.seed(5)
dn_cellID<-unlist(lapply(CT_largN, function(x){
  temp<-metaSeu[metaSeu$TSub_F %in% x,]
  keepCell<-sample(row.names(temp), ncellKeep)
  return(keepCell)
}))
seuIn_F<-subset(seuIn, CellID %in% c(otherCellID, dn_cellID))

kValue=20; dValue=30
obj_sce<-as.SingleCellExperiment(seuIn_F)
obj_milo <- Milo(obj_sce)
obj_milo <- buildGraph(obj_milo,  
                       k = kValue,  
                       d = dValue)
set.seed(123)
obj_milo <- makeNhoods(obj_milo, 
                       prop = 0.10, 
                       k = kValue, 
                       d=dValue, 
                       refined = TRUE)
plotNhoodSizeHist(obj_milo) ### check neighbourhood size
obj_milo <- countCells(obj_milo, 
                       meta.data = data.frame(colData(obj_milo)), 
                       samples="tempU")
obj_milo <- calcNhoodDistance(obj_milo, d=dValue) ### this step takes forever
saveRDS(obj_milo, file=sprintf("%s/obj_milo_dnSample_tobeused_forTSubtypes.RDS",saveFolder))

obj_milo<-readRDS(file=sprintf("%s/obj_milo_dnSample_tobeused_forTSubtypes.RDS",saveFolder))
#########
milo_design <- data.frame(colData(obj_milo))[,c("tempU", "Disease2",'Disease_Fig')]
(milo_design <- distinct(milo_design))
milo_design$Disease_Fig<-sapply(strsplit(milo_design$Disease_Fig,split="\\_"),function(x) x[1])
(milo_design <- distinct(milo_design))
colnames(milo_design)[3]<-'DiseaseType'
rownames(milo_design) <- milo_design$tempU

(milo_design <- milo_design[colnames(nhoodCounts(obj_milo)), , drop=FALSE])
milo_design$Disease2<-factor(milo_design$Disease2,levels=c("MG",'NoMG'))
table(milo_design$Disease2,milo_design$DiseaseType)


Fda_results <- testNhoods(obj_milo, 
                          design = ~ Disease2+DiseaseType, 
                          design.df = milo_design)
Fda_results %>%
  arrange(SpatialFDR) %>%
  head() 

obj_milo <- buildNhoodGraph(obj_milo)
plotUMAP(obj_milo) + 
  plotNhoodGraphDA(obj_milo, 
                   Fda_results, 
                   alpha=0.05) +
  plot_layout(guides="collect")

####################################
da_results0 <- annotateNhoods(obj_milo, 
                              Fda_results, 
                              coldata_col = "TSub_F");head(da_results0)
da_results0$CT <- ifelse(da_results0$TSub_F_fraction < 0.7, "Mixed", da_results0$TSub_F)

###CD4.Temra has 22 cell count; too small and not in the list. manually add this in to make consistence
cd4temra_row<-data.frame(list(0,0,0,1,1,3756,1,'CD4.Temra',1,'CD4.Temra'))
colnames(cd4temra_row)<-colnames(da_results0);cd4temra_row

da_results<-rbind(da_results0,cd4temra_row);tail(da_results)

ctOrder<-c('Mixed',"DPT",'CD8a/a','Prolf.T', 'gdT', 'DNT','CD8.Trm','CD8.Tem','CD4.Temra','CD4.Tem(Th1/Th17)',#
           'CRTAM+.gdT','CD8.Temra','CD4.Tcm(Th17)','CD4.Tcm(Th2)','CD4.Tcm(Th0)','Treg',
           'Tfh','CD8.Tcm','CD4.Tem(Th1)','CD8T.Naive','CD4.Naive','T(agonist)','CD8a/b(entry)')
da_results$CT<-factor(da_results$CT,
                      levels=rev(ctOrder))

plotDAbeeswarm(da_results, group.by = "CT")

#### Fig. 3B #######
pdf(sprintf("%s/Tsub_DA_miloR.pdf",FigOut),width=7.87,height=7.87, family="ArialMT")
print(plotDAbeeswarm(da_results, group.by = "CT"))
dev.off()


#################################################
harT4<-readRDS(file=sprintf("%s/annotated_harT4__June2025.RDS",saveFolder))
harT4$Sub_Tcell<-gsub('CD4T','CD4',harT4$Sub_Tcell)
CTlevel4<-c("DPT",'CD8a/b(entry)','CD8a/a',
            'CD8T.Naive','CD8.Tcm','CD8.Tem','CD8.Trm','CD8.Temra',
            'T(agonist)',
            'CD4.Naive','CD4.Tcm(Th0)','CD4.Tcm(Th2)','CD4.Tcm(Th17)',
            'CD4.Tem(Th1/Th17)','CD4.Tem(Th1)','CD4.Temra',
            'Tfh','Treg',
            'gdT','CRTAM+.gdT');length(CTlevel4)
celltypeColor4<-c("#11FFFF",  "#4DAF4A", "#984EA3", "#FF7F00",
                  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
                  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                  "#1B9E77", "#D95F02", "#7570B3", "#E7298A")
harT4$Sub_Tcell<-factor(harT4$Sub_Tcell,levels=CTlevel4)

############## Fig. S6B ########
pdf(sprintf("%s/T_sub_harT4.pdf",FigOut),width=7.87*1.5,height=7.87, family="ArialMT")
print(DimPlot(harT4,group.by='Sub_Tcell',label=T,repel=T,cols=celltypeColor4,raster=F))
dev.off()
tiff(sprintf("%s/T_sub_harT4.tiff",FigOut),width=7.87*1.05,height=7.87, family="ArialMT",res=600, units="in",compression='lzw')
print(DimPlot(harT4,group.by='Sub_Tcell',label=F,raster=F,cols=celltypeColor4) & NoLegend() & NoAxes ())
dev.off()
