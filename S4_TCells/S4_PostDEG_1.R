
library(edgeR)
library(dplyr)
library(ggplot2)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)
library(msigdbr)
library(fgsea)
library(data.table)
library(sva)
library(UpSetR)
require(gridExtra)

FigOut="./manu/output"
dataFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData/PseudoBulk/SubTcell"

tfh<-readRDS(file=sprintf("%s/pseudoBulk_tfh_DEG.RDS",dataFolder))$Tfh;names(tfh)

MGvsNoMG<-tfh$MGvsNoMG
cpm<-tfh$cpm
pT<-tfh$pT
sigGene<-row.names(MGvsNoMG[MGvsNoMG$fdr<0.05,])
pltGe<-data.frame(pT, t(cpm[sigGene,]));head(pltGe)
MGvsNoMG$col<-ifelse(MGvsNoMG$fdr<0.05,'red','black');sum(MGvsNoMG$col=='red')
geneLab=MGvsNoMG[sigGene,]

##########heatmap#####
pltGe$Dis<-factor(pltGe$Dis,levels=c('MG','NoMG'))
pltGe<-pltGe[order(pltGe$Dis),];head(pltGe)
plotDat_sig<-pltGe[,8:ncol(pltGe)];head(plotDat_sig)

ann_colors = list(
  Disease=c(MG="tomato",NoMG="lightblue"))

top_annotation <- HeatmapAnnotation(
  Disease = factor(pltGe$Dis),
  col = ann_colors
)

library(ComplexHeatmap)
proScale<-t(scale(plotDat_sig));head(proScale)
(mini<-min(proScale));(max<-max(proScale))
negLen<-abs(mini)*100
posLen<-abs(max)*100
mycols<-c(colorRampPalette(colors = c("#1b8a5a","white"))(negLen),
          colorRampPalette(colors = c("white","#ee3e32"))(posLen))

(heatTFG<-ComplexHeatmap::Heatmap(
  proScale,
  col = mycols,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  top_annotation= top_annotation,  
  column_names_gp = gpar(fontsize = 10),  
  show_column_names = TRUE,
  row_order = NULL,
  column_order = NULL
))

################### Fig. 3E ###################
pdf(sprintf("%s/heat_tfh_DEG.pdf",FigOut),width=10,height=10)
print(heatTFG)
dev.off()


################################################
library(ggplot2)
library(RColorBrewer)
require(gridExtra)
head(pltGe[,1:10])

bp_Plt<-function(GeneID="CD1C",colo=c('tomato',"lightblue")){
  df.plot<-data.frame(exp=as.numeric(pltGe[,GeneID]),group=pltGe$Dis)
  bp_col<-colo
  
  set.seed(123)
  p <- ggplot(df.plot, aes(x=group, y=exp,col=group)) + 
    geom_boxplot()+
    geom_jitter(shape=1, position=position_jitter(0.3),size=3)+ #,colour = df.plot$jitCol
    scale_color_manual(values=bp_col)+
    labs(y = "CPM")+
    labs(x = "")+
    ggtitle(sprintf("%s",GeneID))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.5,size=12,colour = "black"),
          axis.title=element_text(size=12,colour = "black"),
          axis.text.y = element_text(size=12,colour = "black"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none") 
  return(p)
}

print(bp_Plt(GeneID="CHGB",colo=c('tomato','lightblue')))

################### Fig. 3E ###################
width=3;height=3
lapply(c("CHGB","TNFSF13B"),function(x){
  pdf(sprintf("%s/tfh_BoxPlt_%s.pdf", FigOut,x), width=width,height=height)
  p<-bp_Plt(GeneID=x,colo=c('tomato','lightblue'))
  print(p)
  dev.off()
})


bp<-tfh$bp

pltGe2<-data.frame(pT, t(cpm[c('TNFRSF13C','TNFRSF13B','TNFRSF17'),]));head(pltGe2)
pltGe2$Dis<-factor(pltGe2$Dis,levels=c('MG','NoMG'))
pltGe2<-pltGe2[order(pltGe2$Dis),];head(pltGe2)

bp_Plt2<-function(GeneID="CD1C",colo=c('tomato',"lightblue")){
  df.plot<-data.frame(exp=as.numeric(pltGe2[,GeneID]),group=pltGe2$Dis)
  bp_col<-colo
  
  set.seed(123)
  p <- ggplot(df.plot, aes(x=group, y=exp,col=group)) + 
    geom_boxplot()+
    geom_jitter(shape=1, position=position_jitter(0.3),size=3)+ #,colour = df.plot$jitCol
    scale_color_manual(values=bp_col)+
    labs(y = "CPM")+
    labs(x = "")+
    ggtitle(sprintf("%s",GeneID))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.5,size=12,colour = "black"),
          axis.title=element_text(size=12,colour = "black"),
          axis.text.y = element_text(size=12,colour = "black"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none") 
  return(p)
}


MGvsNoMG[c('TNFRSF13C','TNFRSF13B','TNFRSF17'),]

width=3;height=3
lapply(c('TNFRSF13C','TNFRSF13B','TNFRSF17'),function(x){
  pdf(sprintf("%s/tfh_BoxPlt_%s.pdf", FigOut,x), width=width,height=height)
  p<-bp_Plt2(GeneID=x,colo=c('tomato','lightblue'))
  print(p)
  dev.off()
})

#################################################################################
################## single cell validation       #################################
#################################################################################

library(Seurat)
saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"

harT2<-readRDS(file=sprintf("%s/annotated_T_NK_ILC__June2025.RDS",saveFolder))
harT2$Disease_Fig<-harT2$Disease_Fig
table(harT2$Disease_Fig,harT2$Disease2)

harT2$ptDis2<-sprintf("%s_%s",harT2$patient,harT2$Disease_Fig)
Idents(harT2)<-'TSub_F'
harT2$DisStudy<-sprintf("%s_%s",harT2$Study,harT2$Disease_Fig)
table(harT2$DisStudy)
harT2_1<-subset(harT2, DisStudy %in% c("NU_Others","Xin_Carcinoma"),invert=T)

unique(harT2_1$DisStudy)
harT2_1$DisStudy<-factor(harT2$DisStudy,levels=c("NU_Thymus_MG", "NU_Thymus_no_MG", 
                                                 "NU_Thymoma_MG","NU_Thymoma_no_MG","Yasumizu_Thymoma_MG"))
(f_chgb<-FeaturePlot(harT2_1, 
                     features='CHGB',
                     col=c('grey','red'),
                     order = TRUE,
                     raster=F,
                     ncol=3,
                     split.by='DisStudy') & NoLegend() & NoAxes())

################### Fig. 3F ###################
pdf(sprintf("%s/FP_tfh_chgb.pdf", FigOut), width=15,height=4)
print(f_chgb)
dev.off()
tiff(sprintf("%s/FP_tfh_chgb.tiff",FigOut),width=15,height=4, family="ArialMT",res=600, units="in",compression='lzw')
print(f_chgb)
dev.off()


(f_chgb2<-FeaturePlot(harT2_1, 
                      features='CHGB',
                      col=c('grey','red'),
                      order = TRUE,
                      raster=F,
                      ncol=3,
                      split.by='DisStudy') & 
    theme(legend.position = "right"))


pdf(sprintf("%s/FP_tfh_chgb_version2.pdf", FigOut), width=15,height=4)
print(f_chgb2)
dev.off()

forDT<-subset(harT2_1, TSub_F %in% c("Tfh"));Idents(forDT)<-"DisStudy"
library(ggplot2)
dplot<-DotPlot(forDT,features=c('TNFSF13B',"CHGB"), scale=T, cols=c('grey','red'),
               scale.by='size',col.max=1,scale.max=15,dot.scale=10)+ 
  coord_flip()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

################### Fig. 3F ###################
pdf(sprintf("%s/DotPlot_tfh_chgb.pdf", FigOut), width=8,height=3.5)
print(dplot)
dev.off()


