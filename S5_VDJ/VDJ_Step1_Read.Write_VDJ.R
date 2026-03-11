######################
library(djvdj)
library(Seurat)
library(dplyr)
library(ggplot2)
##############
saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"
FigOut="./Thymus/FigOut"
####################


path_Jun2023<-"./Projects/SC/Thymus/June2023_nucore/VDJ"
path_Sep2023<-"./Projects/SC/Thymus/Sep2023_nucore/VDJ"
path_Oct2023<-"./Projects/SC/Thymus/Oct2023_nucore/VDJ"
path_Dec2023<-"./Projects/SC/Thymus/Dec2023_nucore/VDJ"

pathList<-list(path_Jun2023=path_Jun2023,
               path_Sep2023=path_Sep2023,
               path_Oct2023=path_Oct2023,
               path_Dec2023=path_Dec2023)


#############################################################################
######     load annotated B cell subtypes of scRNAseq    ####################
#############################################################################
all.B<-readRDS(file=sprintf("%s/B_withSubtype_Annotated_May2025.RDS",saveFolder))
all.B$barcode<-sapply(strsplit(row.names(all.B@meta.data),split="\\_"),function(x) x[1])
seuB_NU0<-subset(all.B,GSE %in% c("Xin",'Yasumizu','ThymusMay242023'),invert=T)

### remove one sample which has no V(D)J sequencing
seuB_NU<-subset(seuB_NU0,ID %in% c("PT9.S13.MetastasisTumor2.Mixed.vdj_b"),invert=T) 

###########
(selPID<-names(table(seuB_NU$ID)[table(seuB_NU$ID)>10]))
BcrDat_djvdj<-unlist(lapply(1:length(pathList),function(iStudy){
  fileCount<-sprintf("%s/pBCR.txt",pathList[[iStudy]])
  if(file.exists(fileCount)){
    pTemp<-read.delim(fileCount,header=F,sep="\t")
  }else{
    pTemp<-read.delim(sprintf("%s/pCount.txt",pathList[[iStudy]]),header=F,sep="\t")
  }
  pTemp<-pTemp[sprintf("VDJ_%s",pTemp$V1) %in% selPID,,drop=F]
  
  readData<-lapply(1:nrow(pTemp),function(iSample){
    sampleName<-pTemp$V1[iSample];print(sampleName)
    tempPath<-sprintf('%s/VDJ_%s/outs/per_sample_outs/VDJ_%s/vdj_b',pathList[[iStudy]],sampleName,sampleName)
    tempSeu<-subset(seuB_NU,ID==sprintf("VDJ_%s",sampleName))
    tempSeu_vdj <- tempSeu |> import_vdj(vdj_dir = file.path(tempPath),include_mutations = TRUE)
    return(tempSeu_vdj)
  }) %>% `names<-`(pTemp$V1)
  return(readData)
}),recursive=FALSE)
names(BcrDat_djvdj)

##############
BcrDat_djvdj_meta<-lapply(BcrDat_djvdj,function(x) return(x@meta.data));names(BcrDat_djvdj_meta)
lapply(BcrDat_djvdj_meta,function(x) sum(colnames(x)==colnames(BcrDat_djvdj_meta[[1]])))
rbindMeta<-do.call('rbind',BcrDat_djvdj_meta);dim(rbindMeta)
rbindMetaMerg<-rbindMeta[,which(colnames(rbindMeta)=='clonotype_id'):ncol(rbindMeta)]
rbindMetaMerg$mergID<-sapply(strsplit(row.names(rbindMetaMerg),split="\\."),function(x) x[2])
nuMeta<-seuB_NU@meta.data
nuMeta$mergID<-row.names(nuMeta)
sum(rbindMetaMerg$mergID %in% nuMeta$mergID)
mergMeta<-dplyr::left_join(nuMeta,rbindMetaMerg,by="mergID")
sum(is.na(mergMeta$clonotype_id))
row.names(mergMeta)<-mergMeta$mergID

seuB_NU@meta.data<-mergMeta

seuB_NU <- define_clonotypes(
  seuB_NU,
  data_cols = "cdr3",         
  clonotype_col = "clonotype_id_F"
)
saveRDS(seuB_NU,file=sprintf("%s/seuB_NU_readVDJ_fromDJVDJ.RDS",saveFolder))



#############################################################################
##################   for T cells      ######################################
#############################################################################
all.T<-readRDS(file=sprintf("%s/annotated_T_NK_ILC__June2025.RDS",saveFolder))
all.T$barcode<-sapply(strsplit(row.names(all.T@meta.data),split="\\_"),function(x) x[1])
seuT_NU0<-subset(all.T,GSE %in% c("Xin",'Yasumizu','ThymusMay242023'),invert=T)
seuT_NU<-subset(seuT_NU0,ID %in% c("PT9.S13.MetastasisTumor2.Mixed.vdj_b"),invert=T)

(selPID<-names(table(seuT_NU$ID)))
TcrDat_djvdj<-unlist(lapply(1:length(pathList),function(iStudy){
  fileCount<-sprintf("%s/pTCR.txt",pathList[[iStudy]])
  if(file.exists(fileCount)){
    pTemp<-read.delim(fileCount,header=F,sep="\t")
  }else{
    pTemp<-read.delim(sprintf("%s/pCount.txt",pathList[[iStudy]]),header=F,sep="\t")
  }
  pTemp<-pTemp[sprintf("VDJ_%s",pTemp$V1) %in% selPID,,drop=F]
  
  readData<-lapply(1:nrow(pTemp),function(iSample){
    sampleName<-pTemp$V1[iSample];print(sampleName)
    tempPath<-sprintf('%s/VDJ_%s/outs/per_sample_outs/VDJ_%s/vdj_t',pathList[[iStudy]],sampleName,sampleName)
    tempSeu<-subset(seuT_NU,ID==sprintf("VDJ_%s",sampleName))
    tempSeu_vdj <- tempSeu |> import_vdj(vdj_dir = file.path(tempPath),include_mutations = TRUE)
    return(tempSeu_vdj)
  }) %>% `names<-`(pTemp$V1)
  return(readData)
}),recursive=FALSE)
names(TcrDat_djvdj)

##############
TcrDat_djvdj_meta<-lapply(TcrDat_djvdj,function(x) return(x@meta.data))
lapply(TcrDat_djvdj_meta,function(x) sum(colnames(x)==colnames(TcrDat_djvdj_meta[[1]])))
rbindMeta<-do.call('rbind',TcrDat_djvdj_meta)
rbindMetaMerg<-rbindMeta[,which(colnames(rbindMeta)=='clonotype_id'):ncol(rbindMeta)]
rbindMetaMerg$mergID<-sapply(strsplit(row.names(rbindMetaMerg),split="\\."),function(x) x[2])
nuMeta<-seuT_NU@meta.data
nuMeta$mergID<-row.names(nuMeta)
sum(rbindMetaMerg$mergID %in% nuMeta$mergID)
mergMeta<-dplyr::left_join(nuMeta,rbindMetaMerg,by="mergID")
sum(is.na(mergMeta$clonotype_id))
row.names(mergMeta)<-mergMeta$mergID

seuT_NU@meta.data<-mergMeta

seuT_NU <- define_clonotypes(
  seuT_NU,
  data_cols = "cdr3",         
  clonotype_col = "clonotype_id_F"
)
saveRDS(seuT_NU,file=sprintf("%s/seuT_NU_readVDJ_fromDJVDJ.RDS",saveFolder))



