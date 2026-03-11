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

dataFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData/PseudoBulk/SubEP"
list.files(dataFolder)

allCT<-list.files(dataFolder,pattern="sampleMeta_");allCT
CT<-gsub("sampleMeta_","",gsub(".csv","",allCT));CT 

############for coding gene ##########
anoFile<-read.delim("./work/reference/RNAseqRef/RefDownload/mart_export_ensemble110_GRCh38.p14.txt",header=T,sep=",",stringsAsFactors=FALSE)
anoFile[grepl("^MT\\-",anoFile$Gene.name),]
anoFile[grepl("-",anoFile$Gene.name),1:4]
anoFile$Gene.name<-gsub("-",".",anoFile$Gene.name)
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]

geneSets0 <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = "GO:BP")

CT<-"EP"
allT_DEG<-lapply(1:length(CT),function(iC){
  cellType<-CT[iC];print(cellType)
  
  meta0<-read.csv(sprintf("%s/sampleMeta_%s.csv",dataFolder,cellType));head(meta0)
  dat<-read.delim(sprintf("%s/sampleValues_%s.txt",dataFolder,cellType),header=T,row.names=1,sep="\t",stringsAsFactors = F);dim(dat)
  row.names(meta0)<-colnames(dat)
  meta0$ID<-colnames(dat)
  meta0$patient<-sapply(strsplit(meta0$ptDis2,split="\\_"),function(iIn) iIn[1])
  
  pdata<-meta0[,2:ncol(meta0)];pdata
  pdata$Dis<-sapply(strsplit(pdata$ID,split="\\_"),function(x) x[2]);pdata
  pdata$chem<-ifelse(pdata$patient %in% c("PT1",'PT2','PT3','PT4'),'V3','V5');table(pdata$chem,pdata$Dis)
  pdata$Study<-ifelse(pdata$patient %in% paste("PT", 1:16,sep=""),'NU','Public')
  
  pT<-pdata[pdata$Dis %in% c("MG",'NoMG') & pdata$Study %in% 'NU' & pdata$n_cells>=25,]
  if (sum(table(pT$Dis)>2)==2){
    pT$Dis<-factor(pT$Dis,levels=c('NoMG','MG'))
    datF<-dat[,paste(pT$ID)];dim(datF)
    
    datF0<-datF[rowSums(datF)>0,];dim(datF0)
    datF_b<-datF0[row.names(datF0) %in% codingGene$Gene.name,];dim(datF_b)
    edgeInput1<-DGEList(datF_b,group=as.factor(pT$Dis))
    edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(edgeInput1)*0.1,]
    edgeInput2<- calcNormFactors(edgeInput2)
    cpm<-cpm(edgeInput2,prior.count=1,log=TRUE)
    
    design.mat <- model.matrix(~pT$Dis+pT$chem);design.mat#+pT$chem
    print(is.fullrank(design.mat))
    colnames(design.mat)<-sapply(strsplit(colnames(design.mat),split="\\$"),function(x) x[2]) 
    colnames(design.mat)[1]<-"(Intercept)"
    rownames(design.mat) <- colnames(edgeInput2)
    edgeInput3<- estimateDisp(edgeInput2, design.mat, robust=TRUE) 
    fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)
    MGvsNoMG<-glmQLFTest(fit, coef=2)$table;min(MGvsNoMG$PValue)
    MGvsNoMG$fdr<-p.adjust(MGvsNoMG$PValue,"BH");min(MGvsNoMG$fdr)
    
    outList<-list(
      cpm=cpm,
      pT=pT,
      pdata=pdata,
      MGvsNoMG=MGvsNoMG
    )}
  
  return(outList)
});names(allT_DEG)<-CT

lapply(allT_DEG,function(x) x$MGvsNoMG[x$MGvsNoMG$fdr<0.05,])
saveRDS(allT_DEG, file=sprintf("%s/pseudo_DEG_EP.RDS",dataFolder))



MGvsNoMG<-allT_DEG$EP$MGvsNoMG
cpm<-allT_DEG$EP$cpm
pT<-allT_DEG$EP$pT

chrnGene<-c("CHRNA1",'CHRNA2',"CHRNA3",'CHRNA4',"CHRNA5",'CHRNA6',"CHRNA7",'CHRNA9',"CHRNA10",
            'CHRNE','CHRND','CHRNG','CHRNB1','CHRNB2','CHRNB4')

### none of chrnGene is significant
MGvsNoMG[chrnGene,]

pltGene<-chrnGene[chrnGene %in% row.names(MGvsNoMG)];length(pltGene)
colnames(cpm)==pT$ptDis2
df_plot<-data.frame(pT, t(cpm[pltGene,]));head(df_plot)


ipage.plot<-lapply(1:length(pltGene),function(iGene){
  whichGene<-pltGene[iGene]
  df.plot<-data.frame(exp=as.numeric(df_plot[,whichGene]),
                      df_plot[,c("chem","Dis")])
  df.plot$jitCol<-c('purple','cyan')[as.factor(df.plot$chem)]
  bp_col<-c("grey",brewer.pal(n=5,name="Set1"))
  #df.plot$Dis<-as.character(df.plot$Dis)
  
  set.seed(123)
  p <- ggplot(df.plot, aes(x=Dis, y=exp,col=Dis)) + 
    geom_boxplot()+
    geom_jitter(shape=1, position=position_jitter(0.3),size=4,colour = df.plot$jitCol)+
    scale_color_manual(values=bp_col)+
    labs(y = "CPM")+
    labs(x = "")+
    ggtitle(sprintf("%s--FDR: %s; logFC: %s",whichGene, 
                    round(MGvsNoMG[whichGene,'fdr'],3),
                    round(MGvsNoMG[whichGene,'logFC'],3)))+
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
})


pdf(sprintf("%s/boxplot_EP_pseudobulk.pdf",FigOut),width=7.87*2,height=7.87*2, family="ArialMT")
do.call(grid.arrange, c(ipage.plot, ncol=5,nrow=3))
dev.off()

