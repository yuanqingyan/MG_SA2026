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

SaveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"
list.files(SaveFolder)

har_EP2<-readRDS(file=sprintf("%s/Thymoma.har_EP2_forFig.RDS",SaveFolder))
nmCell_0<-subset(har_EP2, Study %in% 'Yasumizu')
nmCell_0$yalabel<-sapply(strsplit(nmCell_0$Info,split="\\."),function(x) x[6]);table(nmCell_0$yalabel) ## only 23 nmTEC
nmTEC_ID<-nmCell_0@meta.data[nmCell_0@meta.data$yalabel %in% c('nmTEC'),'CellID'];nmTEC_ID

################ Fig S2.B ###########
(dimNMTEC<-DimPlot(har_EP2, 
                   label=F, 
                   group.by="SubCT", 
                   raster=F,
                   cells.highlight= list(nmTEC_ID),
                   cols.highlight = c("red"), cols= "grey") & NoAxes() & NoLegend())

pdf(sprintf("%s/nmTEC_DimPlot.pdf",FigOut),width=7.87*1.1,height=7.87, family="ArialMT")
print(dimNMTEC)
dev.off()

tiff(sprintf("%s/nmTEC_DimPlot.tiff",FigOut),width=7.87*1.1,height=7.87, family="ArialMT",res=600, units="in",compression='lzw')
print(dimNMTEC)
dev.off()


har_EP2_MG<-subset(har_EP2, Disease2 %in% c("MG",'NoMG'))
DefaultAssay(har_EP2_MG)<-'RNA'
chrnGene<-c("CHRNA1",'CHRNA2',"CHRNA3",'CHRNA4',"CHRNA5",'CHRNA6',"CHRNA7",'CHRNA9',"CHRNA10",
            'CHRNE','CHRND','CHRNG','CHRNB1','CHRNB2','CHRNB4')
Idents(har_EP2_MG)<-'SubCT'

allChr<-lapply(1:length(chrnGene),function(x){
  gen<-chrnGene[x]
  vp<-VlnPlot(har_EP2_MG,features=gen,split.by='Disease2')
  return(vp)
})

############## Figure S3  #############
require(gridExtra)
pdf(sprintf("%s/all_chrnGene.pdf",FigOut),width=7.87*3,height=7.87*3, family="ArialMT")
do.call(grid.arrange, c(allChr, ncol=4,nrow=4))
dev.off()



######################################################
####. for pseudobulk RNAseq

library(anndata)
library(sceasy)
toH5ad<-function(obj=obj,GSEName=GSEName,saveFolder=saveFolder_temp){
  matrix_X<-t(obj@assays$RNA@counts)
  obs_in<-obj@meta.data
  varN<-data.frame(name=colnames(matrix_X))
  varN$name<-gsub("-",".",varN$name)
  row.names(varN)<-varN$name
  colnames(matrix_X)<-varN$name
  matrixFinal<-as(matrix_X, "CsparseMatrix")
  umap_coordinates <- Embeddings(obj, reduction = "umap")
  pca_coordinates <- Embeddings(obj, reduction = "pca")

  anndata <- AnnData(
    X = matrixFinal,
    obs = obs_in,
    var = varN)

  anndata$obsm[['X_pca']] <- pca_coordinates
  anndata$obsm[["X_umap"]] <- umap_coordinates

  write_h5ad(anndata, file=sprintf("%s/%s.h5ad",saveFolder,GSEName))
}

har_EP2_MG<-subset(har_EP2, Disease2 %in% c("MG",'NoMG'))
outHarEP2<-har_EP2_MG
outHarEP2$ptDis2<-sprintf("%s_%s",outHarEP2$patient,outHarEP2$Disease2)
toH5ad(obj=outHarEP2,
       GSEName='forPseudoBulk_EP',
       saveFolder=FigOut)


