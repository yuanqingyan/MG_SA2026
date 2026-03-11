library(dplyr)
library(Seurat)
library(Matrix)
library(singleGEO)
library(SoupX)
library(reticulate)
library(anndata)
library(scDblFinder)


minFeature <- 200  #
maxFeature <- 7500
minCount <- 400    #
maxCount <- 40000
maxMT <- 10        #

GSEID<-"Xin"
###################################
CountFolder="./Projects/PubD/Xin_NC_TEC/download/SRR/countData"
outFolder="./Projects/PubD/Xin_NC_TEC/download/SRR/countData/RDS"
saveFolder<-"./Projects/PubD/Thymoma/SeuPublicGSE"
tmpOutput="./SC/Thymus/public/temp"

(AllSample<-list.files(CountFolder,pattern="HRX"))
pData<-read.delim("./SC/Thymus/PublicData/dnSRA/Xin_NC2022/pdata.txt",header=T,sep="\t")
sum(pData$HRX %in% AllSample)

objName=GSEID
numRef=6
###################################
RawDataList_soupx<-sapply(pData$HRX,function(x) {
  sample_i<-x;print(sample_i)
  SCReadFolder<-sprintf("%s/%s/outs",CountFolder,sample_i)
  sc = load10X(SCReadFolder)
  sc = autoEstCont(sc)
  out_forSeu = adjustCounts(sc)
  
  raw_dat <- CreateSeuratObject(counts = out_forSeu,
                                project =sprintf("%s",sample_i),
                                min.cells = 1,
                                min.features = 1)
  set.seed(12345)
  raw_sce <- scDblFinder(GetAssayData(raw_dat, slot="counts"))
  raw_dat$scDblFinder.class=raw_sce$scDblFinder.class
  raw_dat$scDblFinder.score=raw_sce$scDblFinder.score
  
  raw_dat[['percent.mt']] <- PercentageFeatureSet(raw_dat, pattern = '^MT-') #
  raw_dat[['percent.rps']] <- PercentageFeatureSet(raw_dat, pattern = '^RPS') #
  raw_dat[['percent.rpl']] <- PercentageFeatureSet(raw_dat, pattern = '^RPL') #
  raw_dat$percent.rp <- raw_dat$percent.rps + raw_dat$percent.rpl #
  raw_dat[['patient']] <- sample_i
  raw_dat[['ID']] <- sample_i
  
  tiff(sprintf("%s/%s_FeaturePlot_rawData_soupx.tiff",tmpOutput,sample_i),res=300,width=10,height=7,compression = "lzw",unit="in")
  VlnPlot<-VlnPlot(raw_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(VlnPlot)
  dev.off()
  
  return(raw_dat)
})
if(is.null(names(RawDataList_soupx))){names(RawDataList_soupx)<-pData$HRX}
#########################################################################################################
############################   Filter data                           ####################################
#########################################################################################################
lapply(rawIs,function(x) dim(x))
filter.dat<-lapply(rawIs,function(x){
  temp<-subset(x,
               subset = nFeature_RNA > minFeature &
                 nFeature_RNA < maxFeature &
                 nCount_RNA > minCount &
                 nCount_RNA < maxCount &
                 percent.mt < maxMT &
                 scDblFinder.class %in% "singlet")
  filterDat<-GetAssayData(temp,slot="counts")
  return(filterDat)
})

filter.Seu<-lapply(rawIs,function(x){
  temp<-subset(x,
               subset = nFeature_RNA > minFeature &
                 nFeature_RNA < maxFeature &
                 nCount_RNA > minCount &
                 nCount_RNA < maxCount &
                 percent.mt < maxMT &
                 scDblFinder.class %in% "singlet")
  return(temp)
})

rawSeu<-MakeSeuObj_FromRawRNAData(RawList=filter.dat,
                                  GSE.ID=sprintf("%s",objName),
                                  MtPattern='^MT-',
                                  MinFeature=1,
                                  MaxFeature=750000,
                                  MinCount=1,
                                  MaxCount=4000000,
                                  MaxMT=10,
                                  Norm.method="lognorm",
                                  Feature.selection.method = "vst",
                                  Resolution=2,
                                  Nfeatures = 3000)
seu <- SeuObj_integration(Object.list =rawSeu,
                          Object.list2 = NULL,
                          Frow.which.Norm="lognorm",
                          SampleNameAsReference=NULL,
                          NumberOfSampleForReference=numRef,
                          Nfeatures = 3000,
                          Do.scale = TRUE,
                          Do.center = TRUE,
                          Anchor.reduction = 'rpca',
                          Dims.anchor=1:30,
                          Dims.umap=1:30,
                          Resolution=1.5,
                          Future.globals.maxSize = 96000*1024^2)
saveRDS(seu,file=sprintf("%s/%s_Soupx_beforeAnnotation.rds",outFolder,objName))

seuMeta2<-seu@meta.data
seuMeta2$colID<-row.names(seuMeta2)
seuMeta2$HRX<-seuMeta2$ID
seuMeta_m2<-dplyr::left_join(seuMeta2,pData,by="HRX");dim(seuMeta_m2)
row.names(seuMeta_m2)<-seuMeta_m2$colID
seuMeta_m2<-seuMeta_m2[,!colnames(seuMeta_m2) %in% c("colID")]
seu$Info<-sprintf("%s.%s.%s.%s.%s",
                  gsub("\\.","_",seu$Disease),
                  gsub("\\.","_",seu$Gender),
                  gsub("\\.","_",seu$STAGE),
                  gsub("\\.","_",seu$samName),
                  gsub("\\.","_",seu$Age))
saveRDS(seu,file=sprintf("%s/%s_Seu.RDS",saveFolder,GSEID))


####################################################################################################
###############################        to h5ad                     #################################
####################################################################################################

library(reticulate)
library(anndata)
library(sceasy)
library(Seurat)
options(Seurat.object.assay.version = 'v4')
library(Matrix)


seuPath="./Projects/PubD/Thymoma/SeuPublicGSE"

keepColName<-c("orig.ident","nCount_RNA", "nFeature_RNA", "percent.mt", "GSE","ID","Info")
saveFolder<-"./Projects/ManuscriptData/IntegrateThymus/data/toAnndata"
saveFolder_temp<-"./SC/Thymus/public/temp"

################################################################################################
##################     create some functions                 ###################################
################################################################################################
saveFiles<-function(obj=obj,GSEName=GSEName,saveFolder=saveFolder,keepCol=keepColName){
  matrix_X<-t(obj@assays$RNA@counts)
  obs_in<-obj@meta.data[,keepCol]
  varN<-data.frame(name=colnames(matrix_X))
  varN$name<-gsub("-",".",varN$name)
  colnames(matrix_X)<-varN$name
  saveRDS(matrix_X,file=sprintf("%s/%s.matrix.rds",saveFolder,GSEName))
  saveRDS(obs_in,file=sprintf("%s/%s.obs_in.rds",saveFolder,GSEName))
  saveRDS(varN,file=sprintf("%s/%s.varN.rds",saveFolder,GSEName))
}


# ################################################################################################
# ##################            read data and save to temp folder   ##############################
# ################################################################################################
allGSEName<-c("Xin")

lapply(1:length(allGSEName),function(x){
  GSEName<-allGSEName[x]
  seuTemp<-readRDS(sprintf("%s/%s_Seu.RDS",seuPath,GSEName))
  if(!sum(keepColName %in% colnames(seuTemp@meta.data))==length(keepColName)){
    stop(sprintf("%s has wrong meta colnames",GSEName))}else{
      temp1<-saveFiles(obj=seuTemp,
                       GSEName=GSEName,
                       saveFolder=saveFolder_temp,
                       keepCol=keepColName)}
})


list.files(sprintf("%s",saveFolder_temp))
(FindSaveFileName<-sapply(strsplit(list.files(sprintf("%s",saveFolder_temp),pattern=".varN.rds"),split="\\."),function(x) x[1]))

lapply(1:length(FindSaveFileName),function(x){
  GSEN<-FindSaveFileName[[x]];print(GSEN)
  matrixFile<-readRDS(sprintf("%s/%s.matrix.rds",saveFolder_temp,GSEN))
  obsFile<-readRDS(sprintf("%s/%s.obs_in.rds",saveFolder_temp,GSEN))
  obsFile$batchKey<-paste(GSEN,obsFile$orig.ident,sep="_")
  varFile<-readRDS(sprintf("%s/%s.varN.rds",saveFolder_temp,GSEN))
  
  
  varFinal<-varFile
  matrixFinal<-matrixFile
  matrixFinal<-as(matrixFinal, "CsparseMatrix")
  
  
  row.names(varFinal)<-varFinal$name
  
  anndata <- AnnData(
    X = matrixFinal,
    obs = obsFile,
    var = varFinal)
  
  write_h5ad(anndata, file=sprintf("%s/%s.h5ad",saveFolder,GSEN))
})