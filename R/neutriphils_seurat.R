rm(list = ls())
library(dplyr, verbose = F)
library(Seurat, verbose = F)
library(patchwork,verbose = F)
library(ggplot2, verbose = F)
library(ggpubr, verbose = F)  ##theme_pubr
library(harmony,verbose = F)
library(stringr)
library(ggplot2)
library("ggpubr")
library("ggrepel")

#source("./PANDO/PCnum.R")

#setwd("./PANDO/harmony")
set.seed(123)

srat_obj<-readRDS("./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.rds")
###-----
Neucells<-subset(srat_obj,idents ="Neutrophils")

DimPlot(Neucells, reduction = "umap",group.by = "celltype",label = T,label.size = 3)

Neucells <- NormalizeData(Neucells, normalization.method = "LogNormalize", scale.factor = 10000)
Neucells <- ScaleData(Neucells, features = rownames(Neucells))
Neucells <- FindVariableFeatures(Neucells, selection.method = "vst", nfeatures = 2000)
Neucells <- RunPCA(Neucells, features = VariableFeatures(object = Neucells),npcs = 50)
Neucells <- RunHarmony(Neucells,group.by.vars = "orig.ident",plot_convergence = F) #

reduce="harmony"
pcs<-run_PCnum(Neucells,reduction = reduce)

Neucells <- FindNeighbors(Neucells, dims = 1:pcs,reduction = reduce)
Neucells <- FindClusters(Neucells, resolution =1)
Neucells <- RunUMAP(Neucells, dims = 1:pcs,reduction = reduce)

####---
new.cluster.ids<-c("S100A4+ Neu",  # 0 
                   "IL1B+ Neu",  # 1 
                   "TNFAIP6+ Neu",  # 2 
                   "IL1B+ Neu",  # 3  
                   "IL1B+ Neu",  # 4  
                   "IL1B+ Neu",  #5   
                   "GBP1+ Neu",     #6
                   "IL1B+ Neu" #7
                   )  

names(new.cluster.ids) <- levels(Neucells)
Neucells <- RenameIdents(Neucells, new.cluster.ids)
levels(Neucells)<-c("TNFAIP6+ Neu","IL1B+ Neu","S100A4+ Neu","GBP1+ Neu")
Neucells$Neucelltype<-Idents(Neucells)


saveRDS(Neucells, file = "./PANDO/harmony/Neucells/Neucells_harmony_res1.annotaion.rds")

Neucells.markers <- FindAllMarkers(Neucells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Neucells.markers,file="./PANDO/harmony/Neucells/Neucells_harmony_res1.marker.csv")

##
mycolor<-c("#58C9B9","#A593E0","#DE6449","#F6B352") #"#30A9DE",'#6778AE',,"#D499B9"

marker_genes<-c("FCGR3B","CSF3R","S100A8","S100A9","MME","S100A12","S100A4","FCN1","CDA","IL1RN","MARCKS","TNFAIP6","TNFAIP3","CCL4","IL1B","MX1","IFIT2","XAF1","GBP1")

pdf("./PANDO/harmony/Neucells/Neucells_markergene_vlnplot.pdf",width = 3,height = 5)
VlnPlot(Neucells, features=marker_genes,group.by = "Neucelltype",split.by = "Neucelltype",cols = mycolor,stack = T,flip = T)+
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.9),
        axis.text = element_text(size = 10),legend.text = element_text(size = 10),
        strip.text = element_text(size = 8))+
  NoLegend()+
  xlab("")
dev.off()

###

pdf("./PANDO/harmony/Neucells/Neucells_umap.pdf",width = 4,height = 2.5)
DimPlot(Neucells, reduction = "umap", group.by = "Neucelltype",split.by = "sample", label = F,label.size = 1,pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),axis.title = element_text(size = 9),
        legend.title = element_text(size=8),legend.text = element_text(size = 8),
        title = element_text(size = 9),strip.text = element_text(size = 8))+
  scale_color_manual(values = mycolor)
dev.off()



