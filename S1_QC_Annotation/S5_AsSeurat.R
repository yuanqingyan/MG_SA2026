library(Seurat)
library(cowplot)
saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"


creatSeuObj<-function(pythonDataPath="./temp/scVI",
                      projectName="FromAnndata"){
  dataPath=sprintf("%s/outData",pythonDataPath)
  metaPath=sprintf("%s/outMeta",pythonDataPath)
  NorData<-Seurat::Read10X(data.dir=dataPath)
  metaData<-read.csv(sprintf("%s/metadata.csv",metaPath))
  row.names(metaData)<-metaData[,1]
  embed_umap<-read.csv(sprintf("%s/X_umap.csv",metaPath),row.names=1)
  sObj<-Seurat::CreateSeuratObject(counts = NorData,
                                   meta.data=metaData,
                                   project =projectName,
                                   min.cells = 1,
                                   min.features = 1)
  embedding <- as.matrix(embed_umap)
  colnames(embedding) <- c("umap_1", "umap_2")
  sObj[["umap"]] <- CreateDimReducObject(embedding, key = "umap_",assay="RNA")
  return(sObj)
}

####################################################################################
pythonDataPath="./SC/Thymus/output/scVI/output/Seurat"
seu<-creatSeuObj(pythonDataPath=pythonDataPath,
                 projectName="Thymus_Manu")
saveRDS(seu,file=sprintf("%s/seuThymoma.RDS",saveFolder))
####################################################################################

