
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

dataFolder="./Projects/ManuscriptData/IntegrateThymus/data/ModelData/PseudoBulk/SubBcell"


allCT<-list.files(dataFolder,pattern="sampleMeta_");allCT
CT<-gsub("sampleMeta_","",gsub(".csv","",allCT));CT 

############for coding gene ##########
anoFile<-read.delim("./work/reference/RNAseqRef/RefDownload/mart_export_ensemble110_GRCh38.p14.txt",header=T,sep=",",stringsAsFactors=FALSE);head(anoFile)
anoFile[grepl("^MT\\-",anoFile$Gene.name),]
anoFile[grepl("-",anoFile$Gene.name),1:4]
anoFile$Gene.name<-gsub("-",".",anoFile$Gene.name)
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]

geneSets0 <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = "GO:BP")

CT<-c("allB")

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
  
  pT<-pdata[pdata$Dis %in% c("MG",'NoMG') & pdata$Study %in% 'NU' & pdata$n_cells>=30,]
  if (sum(table(pT$Dis)>=2)==2){
    pT$Dis<-factor(pT$Dis,levels=c('NoMG','MG'))
    datF<-dat[,paste(pT$ID)];dim(datF)
    
    datF0<-datF[rowSums(datF)>0,];dim(datF0)
    datF_b<-datF0[row.names(datF0) %in% codingGene$Gene.name,];dim(datF_b)
    edgeInput1<-DGEList(datF_b,group=as.factor(pT$Dis))
    edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(edgeInput1)*0.1,]
    edgeInput2<- calcNormFactors(edgeInput2)
    cpm<-cpm(edgeInput2,prior.count=1,log=TRUE)
    
    if(cellType=='Un_switch'){
      design.mat <- model.matrix(~pT$Dis);design.mat#
    }else{
      design.mat <- model.matrix(~pT$Dis+pT$chem);design.mat#
    }
    
    print(is.fullrank(design.mat))
    colnames(design.mat)<-sapply(strsplit(colnames(design.mat),split="\\$"),function(x) x[2]) 
    colnames(design.mat)[1]<-"(Intercept)"
    rownames(design.mat) <- colnames(edgeInput2)
    edgeInput3<- estimateDisp(edgeInput2, design.mat, robust=TRUE) 
    fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)
    MGvsNoMG<-glmQLFTest(fit, coef=2)$table;min(MGvsNoMG$PValue)
    MGvsNoMG$fdr<-p.adjust(MGvsNoMG$PValue,"BH");min(MGvsNoMG$fdr)

    MGvsNoMG$stat_gsea<-zscoreT(sign(MGvsNoMG$logFC)*sqrt(MGvsNoMG$F),df=glmQLFTest(fit, coef=2)$df.total)
    MGvsNoMG[is.na(MGvsNoMG$stat_gsea),]
    
    forGSEA<-MGvsNoMG[!is.na(MGvsNoMG$stat_gsea),]
    ranks <- forGSEA$stat_gsea;head(ranks)
    names(ranks) <- row.names(forGSEA)
    ranks<-ranks[order(-ranks)];head(ranks)
    stats <- ranks
    
    geneSets <- geneSets0[geneSets0$gene_symbol %in% names(ranks),]
    m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
    
    set.seed(135)
    eaRes <- fgsea(pathways = m_list, stats = stats, nperm = 1e5, minSize = 20)
    eaRes <-eaRes[order(eaRes$padj),]
 
    outList<-list(
      cpm=cpm,
      pT=pT,
      pdata=pdata,
      MGvsNoMG=MGvsNoMG,
      stats=stats,
      m_list=m_list,
      bp=eaRes
    )
    
  }else{
    outList<-list(
      cpm=NULL,
      pT=NULL,
      pdata=NULL,
      MGvsNoMG=NULL,
      stats=NULL,
      m_list=NULL,
      bp=NULL
    )
  }
  return(outList)
});names(allT_DEG)<-CT


lapply(allT_DEG,function(x) x$MGvsNoMG[x$MGvsNoMG$fdr<0.05,])
saveRDS(allT_DEG,file=sprintf("%s/pseudoBulk_allB_DEG.RDS",dataFolder))

