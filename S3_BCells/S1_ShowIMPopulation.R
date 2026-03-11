
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

har_IM<-readRDS(file=sprintf("%s/har_IM_annotated.RDS",saveFolder))
### correct T cell population based on nnotated_T_NK_ILC__June2025.RDS which is from T cell analysis of S13_Ann_IM.T_Assign_TCellSubtype.R
harT2<-readRDS(file=sprintf("%s/annotated_T_NK_ILC__June2025.RDS",saveFolder))
harT2_meta<-harT2@meta.data;head(harT2_meta)

harIM_Meta<-har_IM@meta.data
harIM_Meta$CellID<-row.names(harIM_Meta)

mergMeta<-dplyr::left_join(harIM_Meta, harT2_meta[,c("CellID",'TSub_F')],by='CellID')
row.names(mergMeta)<-mergMeta$CellID
har_IM@meta.data<-mergMeta

har_IM<-subset(har_IM,CellType_Round1 %in% "Doublet",invert=T)

har_IM$pltCT<-dplyr::case_when(
  har_IM$C11 %in% c("1",'14','10','20')~'Prolif.T',
  har_IM$C11 %in% c("3",'23')~'DNT',
  har_IM$C11 %in% c("21",'24','0','25','19','7')~'DPT',
  har_IM$C11 %in% c("11_0",'11_1')~'gdT',
  har_IM$C11 %in% c("13",'4')~'CD4T',
  har_IM$C11 %in% c("15")~'NK',
  har_IM$C11 %in% c("2",'5','6')~'CD8T',
  .default=as.character(har_IM$CellType_Round1)
)

IM_genes <- c("MS4A2",'CPA3',
              "SCT",
              "CLEC9A",
              "CD14","CD68","FOLR2", 
              "FCGR3B","CSF3R", 
              
              "GNLY","SPON2", 
              "IGKC","MZB1","MS4A1",
              "MKI67","TOP2A", 
              
              "TRAC","TRDC", 
              "CD4", 'CD40LG',"CD8A",'CD8B'
)

har_IM$pltCT<-factor(har_IM$pltCT, 
                     levels=(c("Mast", "pDC","DC", "Macrophage",'Neutrophils',
                               'NK','Plasma','B','Prolif.B',
                               "Prolif.T",'DNT',"DPT",
                               'CD4T',
                               'CD8T','gdT')))

### Fig. 1A
pdf(sprintf("%s/Dim_IM.pdf",FigOut),width=7.87*1.1,height=7.87, family="ArialMT")
print(DimPlot(har_IM,group.by='pltCT',label=T,repel=T,raster=T))
dev.off()
tiff(sprintf("%s/Dim_IM.tiff",FigOut),width=7.87*1.1,height=7.87, res=600, units="in",compression='lzw')
print(DimPlot(har_IM,group.by='pltCT',label=F,repel=T,raster=T) & NoAxes() & NoLegend())
dev.off()


har_IM$pltCT<-factor(har_IM$pltCT, 
                     levels=rev(c("Mast", "pDC","DC", "Macrophage",'Neutrophils',
                                  'NK','Plasma','B','Prolif.B',
                                  "Prolif.T",'DNT',"DPT",
                                  'CD4T',
                                  'CD8T','gdT')))
DefaultAssay(har_IM)<-"RNA";Idents(har_IM)<-"pltCT"
dotIM<-Seurat::DotPlot(object =har_IM, 
                       features = IM_genes,
                       cols = c("lightgrey", "red"),
                       scale.by='radius',col.max=1)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


pdf(sprintf("%s/dotplot_IM.pdf",FigOut),width=7.87*1.5,height=7.87/1.5, family="ArialMT")
print(dotIM)
dev.off()


saveRDS(har_IM@meta.data,file=sprintf("%s/IM_CT_forRef.RDS",FigOut))



