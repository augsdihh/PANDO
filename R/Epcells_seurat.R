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

set.seed(123)

srat_obj<-readRDS("./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.rds")
###-----
EpCells<-subset(srat_obj,idents =c("Epithelial cells"))

DimPlot(EpCells, reduction = "umap",group.by = "celltype",label = T,label.size = 3)

EpCells <- ScaleData(EpCells, features = rownames(EpCells))
EpCells <- FindVariableFeatures(EpCells, selection.method = "vst", nfeatures = 2000)
EpCells <- RunPCA(EpCells, features = VariableFeatures(object = EpCells),npcs = 50)
EpCells <- RunHarmony(EpCells,group.by.vars = "orig.ident",plot_convergence = F) #

reduce="harmony"
pcs<-run_PCnum(EpCells,reduction = reduce)

EpCells <- FindNeighbors(EpCells, dims = 1:pcs,reduction = reduce)
EpCells <- FindClusters(EpCells, resolution = 0.2)
EpCells <- RunUMAP(EpCells, dims = 1:pcs,reduction = reduce)

####---
new.cluster.ids<-c("MALAT1+ ECs",  # 0
                   "Basal ECs",  # 1
                   "Secretory ECs",  # 2
                   "CDH19+ ECS"   # 3
)  
names(new.cluster.ids) <- levels(EpCells)
EpCells <- RenameIdents(EpCells, new.cluster.ids)
levels(EpCells)<-c("Secretory ECs","MALAT1+ ECs","Basal ECs", "CDH19+ ECS")
EpCells$EPcelltype<-Idents(EpCells)

saveRDS(EpCells, file = "./PANDO/harmony/Epcells/Epcells_harmony_res0.2.annotaion.rds")

Epcells.markers <- FindAllMarkers(EpCells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Epcells.markers,file="./PANDO/harmony/Epcells/EpCells_harmony_res0.2.marker.csv")

#######---Dimplot
EpCells<-readRDS("./PANDO/harmony/Epcells/Epcells_harmony_res0.2.annotaion.rds")
pdf("./PANDO/harmony/Epcells/Epcells_umap.pdf",width = 4,height = 2.5)
mycolor<-c("#A593E0","#30A9DE","#8CD790","#F6B352") 

DimPlot(EpCells, reduction = "umap",group.by ="EPcelltype",split.by = "sample",label = F,label.size = 1,pt.size = 0.1)+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),axis.title = element_text(size = 9),
        legend.title = element_text(size=8),legend.text = element_text(size = 8),
        title = element_text(size = 9),strip.text = element_text(size = 8))+
  scale_color_manual(values = mycolor)
dev.off()

#######---marks genes
Epcells.markers <-read.csv("./PANDO/harmony/Epcells/EpCells_harmony_res0.2.marker.csv")
top10 <- Epcells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf("./PANDO/harmony/Epcells/Epcells_markergene_dotplot.pdf",width = 3.5,height = 7)
DotPlot(EpCells,features = unique(top10$gene),dot.scale = 4,scale = T,group.by = "EPcelltype")+coord_flip()+  #+NoLegend()
  theme(axis.text.x = element_text(size = 8,angle = 90,vjust = 0.6,hjust = 1),
        axis.text.y = element_text(size = 8),axis.title = element_text(size = 9),
        legend.title = element_text(size=9),legend.text = element_text(size = 8))+
  scale_color_gradient(low="grey90",high ='red')
dev.off()  

##-----cell proportion
EpCells<-readRDS("./PANDO/harmony/Epcells/Epcells_harmony_res0.2.annotaion.rds")

cell.prop<-as.data.frame(prop.table(table(EpCells@meta.data$sample,EpCells@meta.data$EPcelltype),1))
colnames(cell.prop)<-c("sample","cluster","proportion")

pdf("./PANDO/harmony/Epcells/Epcells_cell_distribution.pdf",width = 2,height = 2.5)
ggplot(cell.prop,aes(x=sample,proportion,fill=cluster))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = mycolor) +
  ggtitle("")+
  theme_bw()+
  #  coord_flip()+ 
  theme(axis.title.x  = element_blank(),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.text = element_text(size = 8),legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))
dev.off()

#####---genes1 Dotplot
genes1<-c("TNF","IL1B","CXCL8","CXCL1","CXCL2","CXCL3","MMP2","MMP10")

DotPlot(EpCells,features = genes1,group.by = "sample",dot.scale = 4,scale = T) + coord_flip() + theme_bw()+
  theme(axis.text.y= element_text(size = 7), axis.text.x=element_text(size = 8,hjust = 1,vjust=1,angle = 45),
        legend.title = element_text(size=7),legend.text = element_text(size = 6),
        legend.key.height = unit(15, "pt"),legend.key.width = unit(10, "pt"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradient2(low = "#34a186",mid = 'white',high  = '#b5182b')

###---genes2 Dotplot
genes2<-c("CDKN2A","CDKN1B","TP53","RB1","AKT1","AKT2")
DotPlot(EpCells,features = genes2, assay='RNA',group.by = "sample",dot.scale = 4,scale = T) + coord_flip() + theme_bw()+
  theme(axis.text.y= element_text(size = 7), axis.text.x=element_text(size = 8,hjust = 1,vjust=1,angle = 45),
        legend.title = element_text(size=7),legend.text = element_text(size = 6),
        legend.key.height = unit(15, "pt"),legend.key.width = unit(10, "pt"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradient2(low = "#34a186",mid = 'white',high  = '#b5182b')


##-----DEG of subgroup
EpCells<-readRDS("./PANDO/harmony/Epcells/Epcells_harmony_res0.2.annotaion.rds")
sc<- SetIdent(EpCells, value = "sample")
##---
cluster<-c("Secretory ECs","MALAT1+ ECs","Basal ECs", "CDH19+ ECS")
for(i in 1:length(cluster)){
  DGE_inflamvsearly_inflam <- FindMarkers(subset(sc,EPcelltype==cluster[i]),ident.1 = "inflam",ident.2 = "early-inflam",logfc.threshold = 0.01, min.pct = 0.01,min.cells.group = 1)
  DGE_inflamvsearly_inflam$gene<-rownames(DGE_inflamvsearly_inflam)
  logFC=log2(1.2)
  P.Value_t=0.01
  DGE_inflamvsearly_inflam$change = ifelse(DGE_inflamvsearly_inflam$p_val_adj < P.Value_t & DGE_inflamvsearly_inflam$avg_log2FC < (-logFC),"down",
                                           ifelse(DGE_inflamvsearly_inflam$p_val_adj < P.Value_t & DGE_inflamvsearly_inflam$avg_log2FC > logFC,"up","nosig"))
  
  write.csv(DGE_inflamvsearly_inflam,file = paste0("./PANDO/harmony/Epcells/DGE_inflamvsearly_inflam_",cluster[i],".csv"))
}

##
genes1<-c("MUC5AC","MUC5B","TFF1","TFF3","FCGBP","ZG16B")
mydata<-FetchData(EpCells,vars=c(genes1,'sample',"EPcelltype"))
mydata<-reshape2::melt(mydata)
colnames(mydata)<-c("sample","EPcelltype","gene","expression")
library(gghalves)
ggplot(mydata)+
  geom_half_violin(aes(x=gene,y=expression,split=sample,fill=sample),position = "identity",trim = F)+
  coord_cartesian(ylim = c(0,0.3), expand = F,  default = FALSE, clip = "on")
  #ylim(0,1)
  #geom_half_boxplot(data = mydata[mydata$sample=="inflam",],aes(x=group,y=value),width = 0.15,side = "l",outlier.shape = NA)+
  #geom_half_boxplot(data = histone_1[histone_1$sample=="D2_34D",],aes(x=group,y=value),width = 0.15,side = "r",outlier.shape = NA)+
  labs(x = "diffence of DNA methylatin between D2 and AC", y = "Normalized read density")+
  ggpubr::theme_pubr()

p<-VlnPlot(EpCells,features =  genes1, stack = T,split.by = 'sample',split.plot = T,
        sort = TRUE, cols = c("#78C2C4","#C73E3A")) 

p<-VlnPlot(EpCells,features =  genes1,stack=T,pt.size=0,flip = T,add.noise = T,split.by = 'sample',
           #group.by = "celltype",
           cols = c("#78C2C4","#C73E3A"),
           split.plot = T)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')

genes2<-c("PI3","SLPI","BPIFA1")
genes3<-c("C3","SAA1","SAA2","S100A8","S100A9")
genes4<-c("SAA1","S100A8","S100A9","C3","PI3")
genes5<-c("SAA1","SAA2","S100A8","S100A9","C3")




