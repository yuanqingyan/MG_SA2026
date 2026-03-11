
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

FigOut="./manu/output"
dataFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData/PseudoBulk"
allCT<-list.files(dataFolder,pattern="sampleMeta_");allCT
CT<-gsub("sampleMeta_","",gsub(".csv","",allCT));CT 

################################################################################################################
###########            meta with addition information                                                    #######
################################################################################################################
meta_addi<-read.csv("./Projects/ManuscriptData/IntegrateThymus/data/RDS/thymus_meta.csv")
meta_addi$PtDis<-paste(meta_addi$Info_SampleType,meta_addi$Info_Patient,sep="_")
meta_addi_B<-meta_addi[meta_addi$SubCT %in% c("B","Prolif.B"),c(2:6,8)];dim(meta_addi_B)

#############
cellType<-"Bcell"
meta0<-read.csv(sprintf("%s/sampleMeta_%s.csv",dataFolder,cellType));head(meta0)
dat<-read.delim(sprintf("%s/sampleValues_%s.txt",dataFolder,cellType),header=T,row.names=1,sep="\t",stringsAsFactors = F);dim(dat)
data.frame(meta=meta0$PatientDis,datName=colnames(dat))
row.names(meta0)<-colnames(dat)
meta0$ID<-colnames(dat)

pdata<-meta0[,2:4];pdata
pdata$PT<-sapply(strsplit(pdata$ID,split="\\_"),function(x) x[length(x)]);pdata
pdata$Organ<-sapply(strsplit(pdata$ID,split="\\_"),function(x) x[1]);pdata
pdata$Dis<-sapply(strsplit(pdata$ID,split="\\_PT"),function(x) x[1]);pdata
pdata$chem<-ifelse(pdata$PT %in% c("PT1",'PT2','PT3','PT4'),'V3','V5');table(pdata$chem,pdata$Dis)

pdata$sex<-ifelse(pdata$PT %in% c("PT3",'PT10','PT12','PT13','PT16'),'M','F')
pdata$age<-ifelse(pdata$PT %in% c('PT10','PT14','PT16'),'Old','Young')## 60 is the cutoff

###############
pB<-pdata[pdata$Dis %in% c("Thymus_MG",'Thymus_no_MG','Thymoma_MG','Thymoma_no_MG'),]
pB$Dis2<-pB$Dis
pB$MG<-ifelse(pB$Dis2 %in% c("Thymus_MG",'Thymoma_MG'),'MG','NoMG');table(pB$MG)
table(pB$chem,pB$MG)

datF<-dat[,paste(pB$PatientDis)];dim(datF)

#############################################################################################################
############                        for coding gene                                                ##########
#############################################################################################################
anoFile<-read.delim("./work/reference/RNAseqRef/RefDownload/mart_export_ensemble110_GRCh38.p14.txt",header=T,sep=",",stringsAsFactors=FALSE);head(anoFile)
anoFile$Gene.name<-gsub("-",".",anoFile$Gene.name)
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]

#############################################################################################################
mypalette <- colorRampPalette(brewer.pal(name="Paired", n = 12))(5)
mypalette <- pals::watlington(5)

datF0<-datF[rowSums(datF)>0,];dim(datF0)
datF_b<-datF0[row.names(datF0) %in% codingGene$Gene.name,];dim(datF_b)
edgeInput1<-DGEList(datF_b,group=as.factor(pB$MG))
edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(edgeInput1)*0.2, ]
edgeInput2<- calcNormFactors(edgeInput2)
cpm<-cpm(edgeInput2,prior.count=1,log=TRUE)

pB$MG<-factor(as.character(pB$MG),levels=c("NoMG","MG"))
design.mat <- model.matrix(~pB$MG+pB$chem);design.mat
print(is.fullrank(design.mat))
colnames(design.mat)<-sapply(strsplit(colnames(design.mat),split="\\$"),function(x) x[2])
colnames(design.mat)[1]<-"(Intercept)"
rownames(design.mat) <- colnames(edgeInput2)

###
edgeInput3<- estimateDisp(edgeInput2, design.mat, robust=TRUE) ##similar to estimateCommonDisp and estimateTagwiseDisp.
edgeInput3$common.dispersion;#plotBCV(edgeInput3)
fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)

mypalette <- rev(brewer.pal(length(unique(pB$Dis)),"Paired"))
mypalette[2]<-"purple"

par(mfrow=c(2,2))
plotMDS(cpm, col=mypalette[factor(pB$Dis)],pch=19,cex=2,xlab="Dim1",ylab="Dim2")
plotMDS(cpm, col=mypalette[factor(pB$Dis2)],pch=19,cex=2,xlab="Dim1",ylab="Dim2")
plotMDS(cpm, col=mypalette[factor(pB$Dis)],pch=19,cex=2,xlab="Dim1",ylab="Dim2")
plotMDS(cpm, col=mypalette[factor(pB$MG)],pch=19,cex=2,xlab="Dim1",ylab="Dim2")

plotMDS(cpm, col=mypalette[factor(pB$MG)],pch=19,cex=2,xlab="Dim1",ylab="Dim2")
plotMDS(cpm, col=mypalette[factor(pB$Dis)],pch=19,cex=1,xlab="Dim1",ylab="Dim2",labels=pB$ID)


###########################################################################
### Sensitive analysis:  consider age and sex as convariates  #############
###########################################################################

design.mat2 <- model.matrix(~pB$MG+pB$chem+pB$sex+pB$age);design.mat2
print(is.fullrank(design.mat2))
colnames(design.mat2)<-sapply(strsplit(colnames(design.mat2),split="\\$"),function(x) x[2])
colnames(design.mat2)[1]<-"(Intercept)"
rownames(design.mat2) <- colnames(edgeInput2)

###
edgeInput3_full<- estimateDisp(edgeInput2, design.mat2, robust=TRUE) ##similar to estimateCommonDisp and estimateTagwiseDisp.
edgeInput3_full$common.dispersion;#plotBCV(edgeInput3)
fit_full <- glmQLFit(edgeInput3_full, design.mat2, robust=TRUE)
colnames(fit_full);plotQLDisp(fit_full)

MGvsNoMG_full<-glmQLFTest(fit_full, coef=2)$table;min(MGvsNoMG_full$PValue)
MGvsNoMG_full$fdr<-p.adjust(MGvsNoMG_full$PValue,"BH");min(MGvsNoMG_full$fdr)
sum(MGvsNoMG_full$fdr<0.05)
MGvsNoMG_full[MGvsNoMG_full$fdr<0.05,]
MGvsNoMG_full[c('CD1C','LILRA4'),]


########################################################################################################
####   due to the limited samples size; we do not include the sex and age as covariates   #############
########################################################################################################

MGvsNoMG<-glmQLFTest(fit, coef=2)$table;min(MGvsNoMG$PValue)
MGvsNoMG$fdr<-p.adjust(MGvsNoMG$PValue,"BH");min(MGvsNoMG$fdr)
sum(MGvsNoMG$fdr<0.05)
MGvsNoMG[MGvsNoMG$fdr<0.05,]

par(mfrow=c(1,1))
plot(MGvsNoMG$logFC, -log10(MGvsNoMG$PValue))

#########
sigGene<-row.names(MGvsNoMG[MGvsNoMG$fdr<0.05,])
pltGe<-data.frame(pB, t(cpm[sigGene,]));head(pltGe)
boxplot(LILRA4~MG,data=pltGe)
boxplot(CD1C~MG,data=pltGe)

MGvsNoMG$col<-ifelse(MGvsNoMG$fdr<0.05,'red','black');sum(MGvsNoMG$col=='red')#6
geneLab=MGvsNoMG[sigGene,]


#####################################################################################
#####################################################################################

############## Fig. 1B ########
pdf(sprintf("%s/MA_bulkRNAseq__B_MGvsNoNG.pdf",FigOut), width=8,height=7)
par(mfrow=c(1,1))    
plot(MGvsNoMG$logCPM, MGvsNoMG$logFC,xlab="Log CPM",ylab="Log FC",col=MGvsNoMG$col)
abline(h=0,lty=2,col='grey')
text(geneLab$logCPM,geneLab$logFC,labels=row.names(geneLab),adj=0,col='red')
dev.off()



################################################
library(ggplot2)
library(RColorBrewer)
require(gridExtra)
pltGe<-data.frame(pB, t(cpm));head(pltGe[,1:10])

bp_Plt<-function(GeneID="CD1C",colo=c("light blue",'magenta')){
  df.plot<-data.frame(exp=as.numeric(pltGe[,GeneID]),group=pltGe$MG)
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

print(bp_Plt(GeneID="CD1C",colo=c('light blue','tomato')))


############## Fig. 1C ########

width=3;height=3
lapply(sigGene,function(x){
  pdf(sprintf("%s/MG_BoxPlt_%s.pdf", FigOut,x), width=width,height=height)
  p<-bp_Plt(GeneID=x,colo=c('light blue','tomato'))
  print(p)
  dev.off()
})

#####################################################################################
#####################################################################################

MGvsNoMG$stat_gsea<-zscoreT(sign(MGvsNoMG$logFC)*sqrt(MGvsNoMG$F),df=glmQLFTest(fit, coef=2)$df.total)
MGvsNoMG[is.na(MGvsNoMG$stat_gsea),]

library(msigdbr)
library(fgsea)
library(ggplot2)
library(data.table)
require(gridExtra)

forGSEA<-MGvsNoMG[!is.na(MGvsNoMG$stat_gsea),]
ranks <- forGSEA$stat_gsea;head(ranks)
names(ranks) <- row.names(forGSEA)
ranks<-ranks[order(-ranks)];head(ranks)

stats <- ranks
geneSets <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = "GO:BP")
geneSets <- geneSets[geneSets$gene_symbol %in% names(ranks),]
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)

set.seed(135)
eaRes <- fgsea(pathways = m_list, stats = stats, nperm = 5e4, minSize = 15)
eaRes <-eaRes[order(eaRes$padj),]
sum(eaRes$padj<0.05)
eaRes[eaRes$padj<0.05,]
fwrite(eaRes,file=sprintf("%s/MG_GOBP_GSEA.csv",FigOut))
fwrite(eaRes[eaRes$padj<0.05,],file=sprintf("%s/MG_GOBP_GSEA_Sig.csv",FigOut))


geneSetsReac <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = "CP:REACTOME")
geneSetsReac <- geneSetsReac[geneSetsReac$gene_symbol %in% names(ranks),]
m_listReac <- geneSetsReac %>% split(x = .$gene_symbol, f = .$gs_name)
set.seed(135)
eaResReac <- fgsea(pathways = m_listReac, stats = stats, nperm = 1e5, minSize = 15)
eaResReac <-eaResReac[order(eaResReac$padj),]
sum(eaResReac$padj<0.05)
eaResReac[eaResReac$padj<0.05,]

geneSetsKegg <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = "CP:KEGG")
geneSetsKegg <- geneSetsKegg[geneSetsKegg$gene_symbol %in% names(ranks),]
m_listKegg <- geneSetsKegg %>% split(x = .$gene_symbol, f = .$gs_name)
set.seed(135)
eaResKegg <- fgsea(pathways = m_listKegg, stats = stats, nperm = 1e5, minSize = 15)
eaResKegg <-eaResKegg[order(eaResKegg$padj),]
sum(eaResKegg$padj<0.05)
eaResKegg[eaResKegg$padj<0.05,]


geneSetsH <- msigdbr(species = "Homo sapiens", category = 'H')
geneSetsH <- geneSetsH[geneSetsH$gene_symbol %in% names(ranks),]
m_listH <- geneSetsH %>% split(x = .$gene_symbol, f = .$gs_name)
set.seed(135)
eaResH <- fgsea(pathways = m_listH, stats = stats, nperm = 5e4, minSize = 15)
eaResH <-eaResH[order(eaResH$padj),]
sum(eaResH$padj<0.05)
eaResH[eaResH$padj<0.05,]


################################################
AllOut<-list(MGvsNoMG=MGvsNoMG,
             fit=fit,
             cpm=cpm,
             pB=pB,
             datF0=datF0,
             
             eaRes=eaRes,
             m_list=m_list,
             eaResReac=eaResReac,
             m_listReac=m_listReac,
             eaResKegg=eaResKegg,
             m_listKegg=m_listKegg,
             geneSetsH=geneSetsH,
             m_listH=m_listH,
             stats=stats,
             forGSEA=forGSEA)

saveRDS(AllOut,file=sprintf("%s/Stat_MG_B.PseudoBulk.RDS",dataFolder))
################################################



################################################
selPath<-c("REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
           'REACTOME_HSF1_ACTIVATION',
           'REACTOME_METABOLISM_OF_LIPIDS',
           'REACTOME_FATTY_ACID_METABOLISM',
           'REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE')
eaResReac[eaResReac$pathway %in% selPath,]

selRec<-lapply(1:length(selPath),function(x){
  return(plotEnrichment(m_listReac[[selPath[x]]], stats) +
           labs(title=selPath[x]))
})

############## Fig. 1D ########
pdf(sprintf("%s/bulk_ReacPathway.pdf",FigOut),width=7.87*1.2,height=7.87/1.2, family="ArialMT")
do.call(grid.arrange, c(selRec, ncol=3,nrow=2))
dev.off()



###############################################################
(selGene<-eaResReac[eaResReac$pathway %in% selPath[1],]$leadingEdge[[1]])

inDat<-pltGe[order(as.factor(pltGe$MG)),];head(inDat[,1:10])
plotDat<-inDat[,colnames(inDat) %in% c("MG",selGene)];dim(plotDat)
color.map<-c("green","magenta")[factor(plotDat$MG)]
mini<-min(plotDat[,2:ncol(plotDat)]);max<-max(plotDat[,2:ncol(plotDat)])
len_bk=100
bk<-seq(mini,max,by=(max-mini)/len_bk)
mycols<-c(colorRampPalette(colors = c("darkblue","blue","white"))(length(bk)),
          colorRampPalette(colors = c("white","red","darkred"))(length(bk)))
dist.my <- function(x) as.dist(1-cor(t(x)))
hclust.my <- function(x) hclust(x, method="ward.D2")

library(gplots)
pdf(sprintf("%s/heat_%s.pdf",FigOut,selPath[1]),width=10,height=10)
sidebarcolors <- color.map
heatmap.2(as.matrix(t(plotDat[,2:ncol(plotDat)])),Rowv=TRUE,Colv=FALSE,cexRow=1.2,cexCol=1.2,
          distfun=dist.my,hclustfun = hclust.my,ColSideColors=sidebarcolors,
          dendrogram=c("row"),col=mycols,key=TRUE,keysize=0.6,symkey=TRUE,scale="row",
          density.info="none",trace="none",margins=c(8,14),main="")
dev.off()


