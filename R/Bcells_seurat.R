rm(list = ls())
library(dplyr, verbose = F)
library(Seurat, verbose = F)
library(patchwork,verbose = F)
library(ggplot2, verbose = F)
library(ggpubr, verbose = F)  ##theme_pubr
library(harmony,verbose = F)
library(stringr)
library(ggplot2)
library(ggsci)
library("ggpubr")
library("ggrepel")

srat_obj<-readRDS("./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.rds")
###-----
Bcells<-subset(srat_obj,idents =c("B cells"))

DimPlot(Bcells, reduction = "umap",group.by = "celltype",label = T,label.size = 3)

Bcells <- NormalizeData(Bcells, normalization.method = "LogNormalize", scale.factor = 10000)
Bcells <- FindVariableFeatures(Bcells, selection.method = "vst", nfeatures = 3000)
Bcells <- ScaleData(Bcells, features = rownames(Bcells))
Bcells <- RunPCA(Bcells, features = VariableFeatures(object = Bcells),npcs = 100)
Bcells <- RunHarmony(Bcells,group.by.vars = "orig.ident",plot_convergence = F)

reduce="harmony"
pcs<-run_PCnum(Bcells,reduction = reduce)

Bcells <- FindNeighbors(Bcells, dims = 1:pcs,reduction = reduce)
Bcells <- FindClusters(Bcells, resolution = 0.4)
Bcells <- RunUMAP(Bcells, dims = 1:pcs,reduction = reduce)

####---
new.cluster.ids<-c("TMX4+ Bm",  # 0
                   "FOSB+ Bm",   # 1
                   "B_Follicular",  #2
                   "TMX4+ Bm",  #3
                   "EGR+ Bm",  # 4
                   "RGD+ Bm"  #5
                   )  
names(new.cluster.ids) <- levels(Bcells)
Bcells <- RenameIdents(Bcells, new.cluster.ids)
levels(Bcells)<-c("TMX4+ Bm","FOSB+ Bm","EGR+ Bm","RGD+ Bm","B_Follicular")
Bcells$Bcelltype<-Idents(Bcells)
saveRDS(Bcells, file = "./PANDO/harmony/Bcells/Bcells_harmony_res0.4.annotaion.rds")

Bcells.markers <- FindAllMarkers(Bcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Bcells.markers,file="./PANDO/harmony/Bcells/Bcells_harmony_res0.4.marker.csv")

#######---marks
Bcells<-readRDS("./PANDO/harmony/Bcells/Bcells_harmony_res0.4.annotaion.rds")

genes<-c("CD19","CD79A","MS4A1","CXCR4", "IGHG3", "IGHG1", "IGHD","CD22", "FCER2","IL4R","IER5","FOSB","FOS","EGR3","NR4A3","RGS2","RGS1","DUSP4")
VlnPlot(Bcells, features=genes,group.by = "Bcelltype",split.by = "Bcelltype",cols = mycolor,stack = T,flip = T)+
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.9),
        axis.text = element_text(size = 10),legend.text = element_text(size = 10),
        strip.text = element_text(size = 8))+
  xlab("")


#Bcells<-readRDS("D:/project/scRNAseq/XiaoCaixia_Cellranger_23.1.29/harmony/Bcells/Bcells_harmony_res0.4.annotaion.rds")
data1<-Bcells@meta.data
data<-Embeddings(Bcells,'umap') %>% as.data.frame()
mydata<-merge(data,data1[,c('sample','Bcelltype')],by=0,all.x=T) %>% column_to_rownames("Row.names")
theme_dot<-theme(panel.grid = element_blank(),
                 panel.background = element_blank(),
                 axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 axis.line = element_blank(),
                 panel.border = element_blank(),
                 legend.text = element_text(size=6),
                 legend.title = element_text(size=8),
                 strip.text = element_text(size=9),#分面背景修改的参数
                 #分面label的背景设置是theme里面的strip.background调整，设置一串element_rect给他
                 plot.title = element_text(size=9))

pdf("./PANDO/harmony/Bcells/Bcells_umap.pdf",width = 4,height = 2.4)
DimPlot(Bcells, reduction = "umap", group.by = "Bcelltype",split.by = "sample", label = F,label.size = 1,pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),axis.title = element_text(size = 9),
        legend.title = element_text(size=8),legend.text = element_text(size = 8),
        title = element_text(size = 9),strip.text = element_text(size = 8))+
  scale_color_manual(values = mycolor)+
  ggtitle("MS4A1+ B cell")
dev.off()

##--cell proportion
cell.prop<-as.data.frame(prop.table(table(Bcells@meta.data$sample,Bcells@meta.data$Bcelltype),1))
colnames(cell.prop)<-c("sample","cluster","proportion")

p1<-ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = mycolor) +
  ggtitle("")+
  theme_bw()+
  #  coord_flip()+ 
  theme(axis.title.x  = element_blank(),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.text = element_text(size = 8),legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))

