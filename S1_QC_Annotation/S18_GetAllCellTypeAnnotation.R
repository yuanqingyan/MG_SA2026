
############################################
har_EnSt3<-readRDS(file="./Projects/ManuscriptData/IntegrateThymus/data/RDS/har_EnSt.RDS");head(har_EnSt3)
har_EP2<-readRDS(file="./Projects/ManuscriptData/IntegrateThymus/data/RDS/Thymoma.har_EP2_forFig.RDS");head(har_EP2)

EnSt3_meta<-har_EnSt3@meta.data[,c("SubCT","Info_Sorting", "Info_Patient", "Info_SampleType","Info_Sample")]
Ep_meta<-har_EP2@meta.data[,c("SubCT","Info_Sorting", "Info_Patient", "Info_SampleType","Info_Sample")]
str_meta<-rbind(Ep_meta,EnSt3_meta);dim(str_meta)


############################################
har_IM<-readRDS(file="./Projects/ManuscriptData/IntegrateThymus/data/RDS/har_IM_annotated.RDS");head(har_IM)
har_IM$SubCT<-har_IM$CellType1

IM_meta<-har_IM@meta.data[,c("SubCT","Info_Sorting", "Info_Patient", "Info_SampleType","Info_Sample")]

thymus_CT<-rbind(str_meta,IM_meta)
thymus_CT$ID<-row.names(thymus_CT)

write.csv(thymus_CT,file="./Projects/ManuscriptData/IntegrateThymus/data/RDS/thymus_meta.csv")

