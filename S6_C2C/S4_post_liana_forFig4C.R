FigOut="./Thymus/FigOut/C2C"
saveFolder="./Projects/ManuscriptData/IntegrateThymus/data/RDS"

loading<-read.csv(sprintf("%s/TensorC2C_lr_loadings.csv",FigOut));head(loading)
sender<-read.csv(sprintf("%s/TensorC2C_SenderCells.csv",FigOut));head(sender)
receiver<-read.csv(sprintf("%s/TensorC2C_ReceiverCells.csv",FigOut));head(receiver)

AllSta<-list(loading=loading, sender=sender, receiver=receiver)
saveRDS(AllSta, file=sprintf("%s/T_B_TensorC2C_stat.RDS",saveFolder))

AllSta<-readRDS(file=sprintf("%s/T_B_TensorC2C_stat.RDS",saveFolder))
loading<-AllSta$loading
sender<-AllSta$sender
receiver<-AllSta$receiver

###################### Fig. 4C #########
loadingPlot<-function(loadF='Factor.1'){
  temp<-loading[order(-loading[,loadF]),]
  pdf(sprintf("%s/TC2C_loading_%s.pdf",FigOut,loadF),width=7.87*1.2,height=7.87/1.5, family="ArialMT")
  plot(1:nrow(temp),temp[,loadF],xlab="rank",ylab='Loading value',main=loadF)
  dev.off()
}

