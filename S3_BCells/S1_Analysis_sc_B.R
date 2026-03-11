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


All.B<-readRDS(file=sprintf("%s/B_withSubtype_Annotated_May2025.RDS",saveFolder))

################### Fig S4A #######
pdf(sprintf("%s/B_bigC.pdf",FigOut),width=7.87*1.05,height=7.87, family="ArialMT")
print(DimPlot(All.B,group.by='BigC',label=T))
dev.off()

tiff(sprintf("%s/B_bigC.tiff",FigOut),width=7.87,height=7.87, family="ArialMT",res=600, units="in",compression='lzw')
print(DimPlot(All.B,group.by='BigC',label=F) & NoLegend() & NoAxes ())
dev.off()

Idents(All.B)<-'BigC'
DotPlot(All.B,features=c("CCR7","CD38",'PRDM1'))
pdf(sprintf("%s/B_bigC_dotplot.pdf",FigOut),width=7.87/2*1.5,height=7.87/2, family="ArialMT")
print(DotPlot(All.B,features=c("CCR7","CD38",'PRDM1'))+theme(axis.text.x = element_text(angle = 90)))
dev.off()


#############################################################################
#############################################################################
#############################################################################
Idents(All.B)<-"Bsub"
BsubLevel<-c("Centroblast_GC","Centrocyte_GC","Early_GC","Mature_GC",
             "NaivB", "Naiv_Cir_1", "Naiv_Cir_2", "Naiv_Cir_3", "Naiv_Cir_4", "Naiv_Cir_5",
             "Un_switch","Switch",
             "MemB_Cir","MemB",
             "MemB.ISG",
             "MemB.Rib+Mit+","MemB.Rib+","MemB.Mit+","MemB.EBI3+",
             "Plasma")

celltypeColor<-c("#8A9D12", "#00FF00", "#0000FF", "#FF0000", "#11FFFF", "#FF00FF", "#FFA500", "#800080", "#FFC0CB", "#A52A2A",
                 "#B0F111", "#308080", "#FFD700", "#4682B4", "#DC143C", "#00CED1", "#ADFF2F", "#FF1493", "#2E8B57", "#7B68EE")

NewBgeneDot<-c("AICDA", 
               'MKI67', 
               "IGLL1","IGLC1", 
               'IGHD','IGHM', 
               'TCL1A','CD72',"FCRL1",
               "IGHV3.53",'IGHV3.66', 
               "IGHV3.20","IGHV3.43",
               "IGHV1.2",'IGHV1.14','IGHV1.46',
               "IGHV1.69",'IGHV1.69D',
               "IGHV3.48",'IGHV3.11','IGHV3.21', 
               'CD69','CD83', 
               'CD27', 
               'IFIT3','IFI44L', 'ISG15',
               'percent.mt','percent.ribo', 
               'RPS10','RPL27A','MT.ND1','MT.ND2',
               'EBI3','CCL22','BATF',
               'CD38','IGKC'
);

All.B$Bsub<-factor(All.B$Bsub,levels=BsubLevel)
print(DimPlot(All.B,group.by='Bsub',label=T,repel=T,cols=celltypeColor))

######## Fig S4A #########
pdf(sprintf("%s/B_Bsub.pdf",FigOut),width=7.87*1.05,height=7.87, family="ArialMT")
print(DimPlot(All.B,group.by='Bsub',label=T,repel=T,cols=celltypeColor))
dev.off()
tiff(sprintf("%s/B_Bsub.tiff",FigOut),width=7.87*1.05,height=7.87, family="ArialMT",res=600, units="in",compression='lzw')
print(DimPlot(All.B,group.by='Bsub',label=F,raster=F,cols=celltypeColor) & NoLegend() & NoAxes ())
dev.off()


######## Fig S4A #########
pdf(sprintf("%s/B_Bsub2.pdf",FigOut),width=7.87*1.05,height=7.87, family="ArialMT")
print(DimPlot(All.B,group.by='organ',label=F,raster=F,cols=celltypeColor))
dev.off()
tiff(sprintf("%s/B_Bsub2.tiff",FigOut),width=7.87*1.05,height=7.87, family="ArialMT",res=600, units="in",compression='lzw')
print(DimPlot(All.B,group.by='organ',label=F,raster=F,cols=celltypeColor))
dev.off()
tiff(sprintf("%s/B_Bsub_byOrgan.tiff",FigOut),width=7.87*1.05*2,height=7.87, family="ArialMT",res=600, units="in",compression='lzw')
print(DimPlot(All.B,group.by='Bsub',label=F,split.by="organ",raster=F,cols=celltypeColor) & NoLegend() & NoAxes ())
dev.off()


Idents(All.B)<-"Bsub"
(dotB_all<-DotPlot(All.B,features=(NewBgeneDot),cols = c("lightgrey", "red"),scale=T)+
    theme(axis.text.x = element_text(angle = 90, hjust=1)))#+ coord_flip()

########### Fig. S4B ########
pdf(sprintf("%s/Dotplot_AllB.pdf",FigOut),width=7.87*2.2,height=7.87/1.1, family="ArialMT")
print(dotB_all)
dev.off()


########### Fig. S4C ########
test<-All.B
test$highlight <- ifelse(!test$Bsub %in% c('MemB.Rib+Mit+'), "Other",as.character(test$Bsub));unique(test$highlight)
colPlt<-FeatureScatter(test, feature1 = "percent.mt",feature2 = "percent.ribo",group.by = "highlight") + 
  scale_color_manual(values = c("MemB.Rib+Mit+" = celltypeColor[16], "Other" = "gray")) +theme_minimal()

pdf(sprintf("%s/scatter_hiMithiRip.pdf",FigOut),width=7.87*1.2,height=7.87, family="ArialMT")
print(colPlt)
dev.off()


########### Fig. S4B ########
all.BDis<-subset(All.B,organ=='Tissue')
all.BDis$Disease2<-dplyr::case_when(
  all.BDis$Disease_Fig %in% "Others"~'Others',
  all.BDis$Disease_Fig %in% c("Thymoma_MG","Thymus_MG") ~"MG",
  all.BDis$Disease_Fig %in% c("Thymoma_no_MG","Thymus_no_MG")~'NoMG',
  .default="unknnown"
)

all.BDis$Disease2<-factor(all.BDis$Disease2,levels=rev(c("Others","NoMG","MG")))

##########################################################################################
all.BDis_f<-subset(all.BDis,Disease2 %in% c("MG",'NoMG') & Study %in% 'NU')

########### Fig. 2A ########
tiff(sprintf("%s/B_Bsub_byMG_NUStudy.tiff",FigOut),width=7.87*1.05*2,height=7.87, family="ArialMT",res=600, units="in",compression='lzw')
print(DimPlot(all.BDis_f,group.by='Bsub',split.by="Disease2",raster=F,cols=celltypeColor) & NoLegend() & NoAxes ())
dev.off()

all.BDis_f$tempU<-paste(all.BDis_f$Disease2,all.BDis_f$Info_Patient,sep="___")
(B_tab<-(table(all.BDis_f$tempU,all.BDis_f$Bsub)))
(colSu<-colSums(B_tab))
(roSum<-rowSums(B_tab))
cutoff<-200
ptCutoff<-50
B_tab_F.PT<-B_tab[roSum>=ptCutoff,]

(B_prop_H<-prop.table(B_tab_F.PT,margin=1));sum(B_prop_H[1,])
(B_prop_f<-B_prop_H[,colSu>=cutoff])

plotDat<-t(B_prop_f);head(plotDat)
color.map<-rep(c("#ff6347","#add8e6"),c(8,4))
mini<-min(plotDat);max<-max(plotDat)
len_bk=100
bk<-seq(mini,max,by=(max-mini)/len_bk)
mycols<-c(colorRampPalette(colors = c("darkblue","blue","white"))(length(bk)),
          colorRampPalette(colors = c("white","red","darkred"))(length(bk)))
dist.my <- function(x) as.dist(1-cor(t(x)))
hclust.my <- function(x) hclust(x, method="complete")

library(ggplot2)
library(gplots)
#

########### Fig. 2C ########
pdf(sprintf("%s/hc_prop_B.pdf",FigOut),width=10,height=10)
sidebarcolors <- color.map
heatmap.2(as.matrix(plotDat),
          Rowv=TRUE,
          Colv=FALSE,
          cexRow=1.2,
          cexCol=1.2,
          distfun=dist.my,
          hclustfun = hclust.my,
          ColSideColors=sidebarcolors,
          dendrogram=c("row"),
          col=mycols,
          key=TRUE,
          keysize=0.6,
          symkey=TRUE,
          scale="row",
          density.info="none",
          trace="none",
          margins=c(6,7),
          main="")
dev.off()


########### Fig. S4E ########
pdf(sprintf("%s/hc_prop_B_columnWithDen.pdf",FigOut),width=10,height=10)
sidebarcolors <- color.map
heatmap.2(as.matrix(plotDat),
          Rowv=TRUE,
          Colv=TRUE,
          cexRow=1.2,
          cexCol=1.2,
          distfun=dist.my,
          hclustfun = hclust.my,
          ColSideColors=sidebarcolors,
          dendrogram=c("both"),
          col=mycols,
          key=TRUE,
          keysize=0.6,
          symkey=TRUE,
          scale='row',
          density.info="none",
          trace="none",
          margins=c(6,7),
          main="")
dev.off()


#####################################################################
################# stacked barplot    ################################
#####################################################################
long_B_prop_f <- data.frame(B_prop_f);head(long_B_prop_f)
long_B_prop_f$Var1<-sapply(strsplit(as.character(long_B_prop_f$Var1),split="___",fixed=T),function(x) x[2])
long_B_prop_f$Var1<-factor(long_B_prop_f$Var1,
                           levels=c("PT11", "PT12", "PT15", "PT2", "PT5", "PT6", "PT7",
                                    "PT9", "PT16", "PT3", "PT4", "PT8"))
allCTinheat<-c("Switch",'Un_switch','MemB','MemB_Cir','Plasma','MemB.Rib+Mit+',
               'NaivB','MemB.ISG','MemB.Rib+')
long_B_prop_f$Var2<-factor(long_B_prop_f$Var2,
                           levels=allCTinheat)

barplt<-ggplot(long_B_prop_f, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent)+
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  labs(fill = "Cell Type", x = "", y = "Percentage (%)")

########### Fig. S4F ########
pdf(sprintf("%s/B_composition_barplot.pdf",FigOut),width=8,height=6)
print(barplt)
dev.off()


##########################################################################################################################################
##############################################  miloR to test Differential abundance  ####################################################
##########################################################################################################################################
library(rlang)
library(dplyr)
library(ggplot2)
#library(glmGamPoi)
library(htmltools)
library(ape)
library(Seurat)
library(Matrix)
library(stringr)
library(miloR)
library(scater)
library(scran)
library(patchwork)

FigOut="./Thymus/FigOut"
allCTinheat<-c("Switch",'Un_switch','MemB','MemB_Cir','Plasma','MemB.Rib+Mit+',
               'NaivB','MemB.ISG','MemB.Rib+')
seuIn<-subset(all.BDis_f, Bsub %in% allCTinheat)
seuIn$Disease2<-as.character(seuIn$Disease2)
seuIn$Bsub<-as.character(seuIn$Bsub)

kValue=20; dValue=30
obj_sce<-as.SingleCellExperiment(seuIn)
obj_milo <- Milo(obj_sce)

#### construct KNN graph ##
obj_milo <- buildGraph(obj_milo,  
                       k = kValue,  
                       d = dValue)

#### Defining representative neighbourhoods
set.seed(5)
obj_milo <- makeNhoods(obj_milo, 
                       prop = 0.12, 
                       k = kValue, 
                       d=dValue, 
                       refined = TRUE)
plotNhoodSizeHist(obj_milo) ### check neighbourhood size

###Counting cells in neighbourhoods
obj_milo <- countCells(obj_milo, 
                       meta.data = data.frame(colData(obj_milo)), 
                       samples="tempU")

###Differential abundance testing
milo_design <- data.frame(colData(obj_milo))[,c("tempU", "Disease2")]
(milo_design <- distinct(milo_design))
rownames(milo_design) <- milo_design$tempU

(milo_design <- milo_design[colnames(nhoodCounts(obj_milo)), , drop=FALSE])
milo_design$Disease2<-factor(milo_design$Disease2,levels=c("MG",'NoMG'))

obj_milo <- calcNhoodDistance(obj_milo, d=dValue)
Fda_results <- testNhoods(obj_milo, 
                          design = ~ Disease2, 
                          design.df = milo_design)
Fda_results %>%
  arrange(SpatialFDR) %>%
  head() 

### visualize neighbourhoods displaying DA
obj_milo <- buildNhoodGraph(obj_milo)
plotUMAP(obj_milo) + 
  plotNhoodGraphDA(obj_milo, 
                   Fda_results, 
                   alpha=0.05) +
  plot_layout(guides="collect")

####################################
da_results <- annotateNhoods(obj_milo, 
                             Fda_results, 
                             coldata_col = "Bsub");head(da_results)
da_results$CT <- ifelse(da_results$Bsub_fraction < 0.7, "Mixed", da_results$Bsub)

ctOrder<-c('Mixed',"Switch",'Un_switch','MemB','MemB_Cir','Plasma','MemB.Rib+Mit+',
           'NaivB','MemB.ISG','MemB.Rib+')
da_results$CT<-factor(da_results$CT,
                      levels=rev(ctOrder))

plotDAbeeswarm(da_results, group.by = "CT")

########### Fig. 2C ########
pdf(sprintf("%s/Bsub_DA_miloR.pdf",FigOut),width=7.87,height=7.87, family="ArialMT")
print(plotDAbeeswarm(da_results, group.by = "CT"))
dev.off()





(FP_M<-FeaturePlot(all.BDis_f, 
                   features = c("CD83",'CD70'),
                   col=c('grey','red'),
                   raster=F,
                   min.cutoff = 'q5',
                   max.cutoff='q95',
                   order = TRUE,
                   pt.size=0.05) & NoLegend() & NoAxes())
########### Fig. 2B ########
tiff(sprintf("%s/B_CD83_CD70.tiff",FigOut),width=7.87*2/1.2/1.2,height=7.87/1.2/1.2, family="ArialMT",res=600, units="in",compression='lzw')
print(FP_M)
dev.off()

(FP_IGMIGD<-FeaturePlot(all.BDis_f, 
                        features = c('IGHD',"IGHA1"),
                        col=c('grey','red'),
                        raster=F,
                        min.cutoff = 'q5',
                        max.cutoff='q95',
                        order = TRUE,
                        pt.size=0.05) & NoLegend() & NoAxes())

########### Fig. S4D ########
tiff(sprintf("%s/B_IGHDIGHA1_byMG_NUStudy.tiff",FigOut),width=7.87*2/1.2/1.2,height=7.87/1.2/1.2, family="ArialMT",res=600, units="in",compression='lzw')
print(FP_IGMIGD)
dev.off()



(FP<-FeaturePlot(all.BDis_f, 
                 features = c("CD1C"),
                 col=c('grey','red'),
                 split.by='Disease2',
                 raster=F,
                 order = TRUE,
                 pt.size = 0.5) & NoLegend() & NoAxes())

########### Fig. 3F ########
tiff(sprintf("%s/B_CD1C_byMG_NUStudy.tiff",FigOut),width=7.87*2/1.2/1.2,height=7.87/1.2/1.2, family="ArialMT",res=600, units="in",compression='lzw')
print(FP)
dev.off()



test2<-all.BDis_f
test2$grp<-sprintf("%s(%s)",test2$Bsub,test2$Disease2)

(dt<-DotPlot(
  subset(test2, Bsub %in% c("Switch", "Un_switch")), 
  features = c("CD1C"), 
  group.by = "grp", 
  dot.scale = 8, 
  cols = c("lightgrey", "red"), 
  scale = FALSE  # Turn off per-gene scaling
))

########### Fig. S4G ########
pdf(sprintf("%s/DotPlot__Cd1c_swith.pdf",FigOut),width=7.87,height=7.87/2, family="ArialMT")
print(dt)
dev.off()


# 
####################################################################################
##########################    pathway analysis       ###############################
####################################################################################
library(msigdbr);library(fgsea);library(data.table)
BP<-as.data.frame(msigdbr(species = "Homo sapiens", category ="C5",subcategory = "GO:BP"))[,c("gs_name","human_gene_symbol")];head(BP)
BP<-BP[!BP$gs_name %in% c("GOBP_FLAVONE_METABOLIC_PROCESS",
                          'GOBP_FORMATION_OF_QUADRUPLE_SL_U4_U5_U6_SNRNP',
                          'GOBP_NEGATIVE_REGULATION_OF_INTERLEUKIN_21_PRODUCTION',
                          "GOBP_LEUKOTRIENE_B4_CATABOLIC_PROCESS",
                          'GOBP_NEGATIVE_REGULATION_OF_MYOBLAST_PROLIFERATION'),]
BP_list<-lapply(split(BP,f=BP$gs_name),function(x) x$human_gene_symbol);length(BP_list)#7654

all.BDis_f_BP<-all.BDis_f
for(xpath in 1:length(BP_list)){print(xpath)
  all.BDis_f_BP<-AddModuleScore(
    object=all.BDis_f_BP,
    features=BP_list[xpath],
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    k = FALSE,
    assay = NULL,
    name =names(BP_list)[xpath],
    seed = 1,
    search = FALSE,
    slot = "data")}
saveRDS(all.BDis_f_BP,file=sprintf("%s/temp/all.BDis_f_BP__addRec.RDS",saveFolder))
all.BDis_f_BP<-readRDS(file=sprintf("%s/temp/all.BDis_f_BP__addRec.RDS",saveFolder))

#########################################################################
deg_vsMG<-function(df,istart=30,lmmodel=TRUE){
  result_p<-lapply(istart:ncol(df),function(iPath){print(iPath)
    pathw<-colnames(df)[iPath]
    temp<-as.data.frame(df[,c(which(colnames(df)=="Disease2"),iPath)])
    colnames(temp)[2]<-"Score"
    if(sum(temp$Score==0)==nrow(temp)){
      outdf=NULL
    }else{
      if(lmmodel==TRUE){
        lm = lm(Score ~ Disease2, data = temp)
        outdf<-data.frame(coef=as.data.frame(summary(lm)$coefficients)['Disease2NoMG',"Estimate"], 
                          p=as.data.frame(summary(lm)$coefficients)['Disease2NoMG',"Pr(>|t|)"])
        outdf$comp<-'NoMGvsMG'
        outdf$Pathway<-pathw
      }else{
        p<-wilcox.test(Score~Disease2,data=temp)$p.value
        temp$Disease2<-as.character(temp$Disease2)
        median_MG<-median(temp[temp$Disease2=="MG",'Score'])
        median_NoMG<-median(temp[temp$Disease2=="NoMG",'Score'])
        diff<-median_MG-median_NoMG
        outdf<-data.frame(p=p,median_MG=median_MG,median_NoMG=median_NoMG,diff=diff)
      }
    }
    return(outdf)
  });names(result_p)<-colnames(df)[istart:ncol(df)]
  df_result_p<-do.call('rbind',result_p)
  df_result_p$fdr<-p.adjust(df_result_p$p,"BH");head(df_result_p)
  return(df_result_p)
}

meta_BP<-all.BDis_f_BP@meta.data;colnames(meta_BP)[1:50]
uniCT<-c("Switch",'Un_switch')
list_difBP<-lapply(1:length(uniCT),function(x){
  ct<-uniCT[x]
  temp<-meta_BP[meta_BP$Bsub %in% ct,]
  tempOut<-deg_vsMG(df=temp,istart=which(colnames(meta_BP)=="GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS1"),lmmodel=TRUE)
  return(tempOut)
});names(list_difBP)<-uniCT

df_allBP<-do.call('rbind',lapply(1:length(list_difBP),function(x){
  temp<-list_difBP[[x]]
  temp$CellType<-names(list_difBP)[x]
  return(temp)}))
df_allBP<-df_allBP[order(df_allBP$CellType,df_allBP$fdr),]
write.csv(df_allBP,file=sprintf("%s/df_allBP.csv",FigOut))


list_difBP2<-lapply(1:length(uniCT),function(x){
  ct<-uniCT[x]
  temp<-meta_BP[meta_BP$Bsub %in% ct,]
  tempOut<-deg_vsMG(df=temp,istart=which(colnames(meta_BP)=="GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS1"),lmmodel=FALSE)
  return(tempOut)
});names(list_difBP2)<-uniCT

df_allBP2<-do.call('rbind',lapply(1:length(list_difBP2),function(x){
  temp<-list_difBP2[[x]]
  temp$CellType<-names(list_difBP2)[x]
  return(temp)}))
df_allBP2<-df_allBP2[order(df_allBP2$CellType,df_allBP2$fdr),];head(df_allBP2)
write.csv(df_allBP2,file=sprintf("%s/df_allBP2.csv",FigOut))

twoPath<-c("GOBP_NEGATIVE_REGULATION_OF_REGULATORY_T_CELL_DIFFERENTIATION1",
           "GOBP_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_STIMULUS1")
df_allBP[twoPath,]
df_allBP2[twoPath,]


pltBP<-subset(all.BDis_f_BP, Bsub %in% c("Switch",'Un_switch'))
pltBP$Bsub<-factor(pltBP$Bsub,levels=c("Un_switch",'Switch'))
pltBP$Disease2<-factor(pltBP$Disease2,levels=c('MG',"NoMG"))
v1<-VlnPlot(pltBP,
            group.by='Bsub',
            features=c("GOBP_NEGATIVE_REGULATION_OF_REGULATORY_T_CELL_DIFFERENTIATION1"),
            split.by="Disease2",
            col=c('tomato','light blue'),
            pt.size=0)+
  VlnPlot(pltBP,
          group.by='Bsub',
          features=c("GOBP_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_STIMULUS1"),
          split.by="Disease2",
          col=c('tomato','light blue'),
          pt.size=0)

pdf(sprintf("%s/VlnPlot__sig_BP_twoCT.pdf",FigOut),width=7.87/2,height=7.87*1.2, family="ArialMT")
print(v1)
dev.off()



metaBP=pltBP@meta.data

(pBP1<-ggplot(metaBP, 
              aes(x=Bsub, 
                  y=GOBP_NEGATIVE_REGULATION_OF_REGULATORY_T_CELL_DIFFERENTIATION1, 
                  fill = Disease2)) +
    geom_boxplot(width = .5, 
                 alpha = 1, 
                 median.size = 1.2,
                 show.legend = TRUE) +
    scale_fill_manual(values = c('tomato','light blue'),
                      name = "Disease") +
    ggtitle("Negative regulation of regulatory T cell differentiation")+
    labs(y = "Enrichment score")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1,size=12,colour = "black"),
          axis.title=element_text(size=12,colour = "black"),
          axis.text.y = element_text(size=12,colour = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")))


(pBP2<-ggplot(metaBP, 
              aes(x=Bsub, 
                  y=GOBP_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_STIMULUS1, 
                  fill = Disease2)) +
    geom_boxplot(width = .5, 
                 alpha = 1, 
                 median.size = 1.2,
                 show.legend = TRUE) +
    scale_fill_manual(values = c('tomato','light blue'),
                      name = "Disease") +
    ggtitle("Negative regulation of cellular response to TGF-beta stimulus")+
    labs(y = "Enrichment score")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1,size=12,colour = "black"),
          axis.title=element_text(size=12,colour = "black"),
          axis.text.y = element_text(size=12,colour = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")))

########### Fig. 2G ########
pdf(sprintf("%s/boxPlot__sig_BP_twoCT.pdf",FigOut),width=7.87*2,height=7.87/1.2, family="ArialMT")
print(pBP1+pBP2)
dev.off()


