
library(ggplot2)
library(RColorBrewer)
require(gridExtra)

FigOut="./Thymus/FigOut/C2C"
saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"


genesForPlt=c('CD74','CXCR4','LTB','CD40','CD44','TGFB1','ITGB1','ADAM10','NOTCH2',
              'TNFSF13B','TNFRSF17','ADAM28','ITGA4','CD70','CD27')


dataFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData/PseudoBulk/SubTcell"
allT_DEG<-readRDS(file=sprintf("%s/pseudoBulk_allT_DEG_NoPathway.RDS",dataFolder))
allB_DEG<-readRDS(file="./Projects/ManuscriptData/IntegrateThymus/data/ModelData/PseudoBulk/SubBcell/pseudoBulk_allB_DEG.RDS")



getGeDat<-function(ctDat=tfh){
  MGvsNoMG<-ctDat$MGvsNoMG
  cpm<-ctDat$cpm
  pT<-ctDat$pT
  sigGene<-row.names(MGvsNoMG)
  pltGe<-data.frame(pT, t(cpm[sigGene,]))
  pltGe$Dis<-factor(pltGe$Dis,levels=c('MG','NoMG'))
  return(pltGe)
}

bp_Plt<-function(GeneID="CD1C",pltGe=regT, colo=c('tomato',"lightblue")){
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



B_Rec<-c("CD74",'CD40','CD44','ITGB1','TNFRSF17')
B_Lig<-c("ADAM28")

T_Lig<-c("LTB",'ADAM10')

tomCol<-c('tomato','lightblue')
greCol<-c("deeppink",'green2')

(b_CD74<-bp_Plt(GeneID="CD74",
                pltGe=getGeDat(ctDat=allB_DEG[["Switch"]]),
                colo=c('tomato','lightblue')))
(b_CD74_Un<-bp_Plt(GeneID="CD74",
                   pltGe=getGeDat(ctDat=allB_DEG[["Un_switch"]]),
                   colo=c('tomato','lightblue')))
(b_CD40<-bp_Plt(GeneID="CD40",
                pltGe=getGeDat(ctDat=allB_DEG[["Switch"]]),
                colo=c('tomato','lightblue')))

(b_ITGB<-bp_Plt(GeneID="ITGB1",
                pltGe=getGeDat(ctDat=allB_DEG[["Switch"]]),
                colo=c('tomato','lightblue')))

(b_TSF17<-bp_Plt(GeneID="TNFRSF17",
                 pltGe=getGeDat(ctDat=allB_DEG[["Switch"]]),
                 colo=c('tomato','lightblue')))

(b_adam28<-bp_Plt(GeneID="ADAM28",
                  pltGe=getGeDat(ctDat=allB_DEG[["Switch"]]),
                  colo=c('tomato','lightblue')))

Treg_LTB<-bp_Plt(GeneID="LTB",
                       pltGe=getGeDat(ctDat=allT_DEG[["Treg"]]),
                       colo=greCol)

CD417_LTB<-bp_Plt(GeneID="LTB",
                        pltGe=getGeDat(ctDat=allT_DEG[["CD4.Tcm(Th17)"]]),
                        colo=greCol)
CD41_LTB<-bp_Plt(GeneID="LTB",
                       pltGe=getGeDat(ctDat=allT_DEG[["CD4.Tem(Th1)"]]),
                       colo=greCol)

Tfh_cd40l<-bp_Plt(GeneID="CD40LG",
                        pltGe=getGeDat(ctDat=allT_DEG[["Tfh"]]),
                        colo=greCol)



########## Fig. 4D ###############

NewAllPlt<-list(b_CD74_Un,b_CD74, ### cd 74 in un_switch and switch
                b_TSF17, #TNFRSF17
                b_ITGB,
                b_adam28,
                b_CD40,Treg_LTB,CD417_LTB,CD41_LTB,Tfh_cd40l ## CD40
);length(NewAllPlt)

pdf(sprintf("%s/selectedGeneForC2C_NewAllPlt.pdf",FigOut),width=7.87*1.5,height=7.87*0.75*1.5, family="ArialMT")
print(do.call(grid.arrange, c(NewAllPlt, ncol=5,nrow=2)))
dev.off()



############ Fig. S8C ########

geneInSupFig<-c("HLA.DRA",'CTSS')
(b_HLA<-bp_Plt(GeneID=geneInSupFig[1],
               pltGe=getGeDat(ctDat=allB_DEG[["Switch"]]),
               colo=c('tomato','lightblue')))
(b_CTSS<-bp_Plt(GeneID=geneInSupFig[2],
                pltGe=getGeDat(ctDat=allB_DEG[["Switch"]]),
                colo=c('tomato','lightblue')))

HLA_CTSS_Plt<-list(b_HLA, b_CTSS);length(HLA_CTSS_Plt)


pdf(sprintf("%s/HLA_CTSS_boxplot.pdf",FigOut),width=7.87,height=7.87*0.5, family="ArialMT")
print(do.call(grid.arrange, c(HLA_CTSS_Plt, ncol=2,nrow=1)))
dev.off()




