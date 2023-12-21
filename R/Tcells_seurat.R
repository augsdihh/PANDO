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

source("D:/project/scRNAseq/pipeline/src/01-PCnum.R")
source("D:/project/scRNAseq/pipeline/src/03-sceplot.R")

set.seed(123)

scRNA<-readRDS("./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.rds")
###-----
Tcells<-subset(scRNA,idents ="T cells")

DimPlot(Tcells, reduction = "umap",group.by = "celltype",label = T,label.size = 3)

Tcells <- NormalizeData(Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)
Tcells <- ScaleData(Tcells, features = rownames(Tcells))
Tcells <- RunPCA(Tcells, features = VariableFeatures(object = Tcells),npcs = 50)
Tcells <- RunHarmony(Tcells,group.by.vars = "orig.ident",plot_convergence = F) #

reduce="harmony"
pcs<-run_PCnum(Tcells,reduction = reduce)

Tcells <- FindNeighbors(Tcells, dims = 1:pcs,reduction = reduce)
Tcells <- FindClusters(Tcells, resolution = 0.5)
Tcells <- RunUMAP(Tcells, dims = 1:pcs,reduction = reduce)
####---
new.cluster.ids<-c("CD4+Tcm",  # 0
                   "CD8+Tm",  # 1
                   "CD4+Tem",  # 2
                   "CD8+Tm",  # 3
                   "CD4+Tcm",  # 4
                   "Treg",  #5
                   "CD8+Tm",  #6
                   "NKT",  #7
                   "NK",   #8
                   "CD4+Tex",  #9
                   "Proliferative-T"  #10
                   )  

names(new.cluster.ids) <- levels(Tcells)
Tcells <- RenameIdents(Tcells, new.cluster.ids)
levels(Tcells)<-c("CD8+Tm","CD4+Tcm","CD4+Tem","CD4+Tex","Treg","NKT","NK","Proliferative-T")
Tcells$Tcelltype<-Idents(Tcells)
DimPlot(Tcells, reduction = "umap",label = T,label.size = 3)
DimPlot(Tcells, reduction = "umap",group.by = "Tcelltype",label = T,label.size = 3)

saveRDS(Tcells, file = "./PANDO/harmony/Tcells/Tcells_harmony_res0.5.annotaion.rds")

Tcells.markers <- FindAllMarkers(Tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tcells.markers,file="./PANDO/harmony/Tcells/Tcells_harmony_res0.5.marker.csv")

##
Tcells<-readRDS("./PANDO/harmony/Tcells/Tcells_harmony_res0.5.annotaion.rds")

mycolor<-c("#30A9DE","#58C9B9","#8CD790","#5CAB7D","#A593E0","#DE6449","#D499B9","#F6B352")

pdf("./PANDO/harmony/Tcells/Tcells_umap.pdf",width = 5.5,height = 3)
DimPlot(Tcells, reduction = "umap", group.by = "Tcelltype",split.by = "sample", label = F,label.size = 1,pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),axis.title = element_text(size = 9),
        legend.title = element_text(size=8),legend.text = element_text(size = 8),
        title = element_text(size = 9),strip.text = element_text(size = 8))+
  scale_color_manual(values = mycolor)+
  dev.off()


###---mark genes
marker_genes<-c("CD3E","CD4","CD8A",'IL7R',"CCR7","SELL","CD28","FOXP3","TIGIT","TOX2","PDCD1","TOX","LAG3",
                "IFNG","GNLY","NKG7","GZMB","PRF1","NCAM1","FCGR3A","CDCA5","CDC25A","CCNA2")

pdf("./PANDO/harmony/Tcells/Tcells_makergene_vlnplot.pdf",width = 3,height = 7)
VlnPlot(Tcells, features=marker_genes,group.by = "Tcelltype",split.by = "Tcelltype",cols = mycolor,stack = T,flip = T)+
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.9),
        axis.text = element_text(size = 10),legend.text = element_text(size = 10),
        strip.text = element_text(size = 8))+
  xlab("")
dev.off()

##-----DEG
Tcells<-readRDS("./PANDO/harmony/Tcells/Tcells_harmony_res0.5.annotaion.rds")
sc<- SetIdent(Tcells, value = "sample")

##---
cluster<-c("CD8+Tm", "CD4+Tcm","CD4+Tem", "CD4+Tex","Treg","NKT","NK","Proliferative-T")
for(i in 1:length(cluster)){
  DGE_inflamvsearly_inflam <- FindMarkers(subset(sc,Tcelltype==cluster[i]),ident.1 = "inflam",ident.2 = "early-inflam",logfc.threshold = 0.01, min.pct = 0.01)
  DGE_inflamvsearly_inflam$gene<-rownames(DGE_inflamvsearly_inflam)
  logFC=log2(1.2)
  P.Value_t=0.01
  DGE_inflamvsearly_inflam$change = ifelse(DGE_inflamvsearly_inflam$p_val_adj < P.Value_t & DGE_inflamvsearly_inflam$avg_log2FC < (-logFC),"down",
                                           ifelse(DGE_inflamvsearly_inflam$p_val_adj < P.Value_t & DGE_inflamvsearly_inflam$avg_log2FC > logFC,"up","nosig"))
  
  write.csv(DGE_inflamvsearly_inflam,file = paste0("./PANDO/harmony/Tcells/DGE_inflamvsearly_inflam_",cluster[i],".csv"))
  
}

DGE_inflamvsearly_inflam<-read.csv("./PANDO/harmony/Tcells/DGE_inflamvsearly_inflam_CD4+Tcm.csv")
p<-ggplot(data = DGE_inflamvsearly_inflam, aes(x = avg_log2FC, y = -log10(p_val_adj), colour=change)) +theme_pubr()+
  geom_point(alpha=1, size=1.5) +
  scale_color_manual(values=c("#4cb1c4", "grey","#b5182b")) +
  #xlim(c(-2.5, 2)) + ylim(c(0,100)) +
  geom_vline(xintercept=c(-log2(1.2),log2(1.2)),lty=2,col="grey2",lwd=0.5) +
  geom_hline(yintercept = -log10(0.01),lty=2,col="grey2",lwd=0.5) +
  labs(x="log2(fold change)",y="-log10(adjusted p-value)",title="inflamvsearly_inflam in CD4+Tcm") +
  theme(plot.title = element_text(hjust = 0.5,size = 10), legend.position="none",
        axis.text = element_text(size = 8), axis.title = element_text(size = 10))
p<-p + geom_text_repel(data = subset(DGE_inflamvsearly_inflam, DGE_inflamvsearly_inflam$p_val_adj < 0.01 & abs(DGE_inflamvsearly_inflam$avg_log2FC) >= log2(1.2)),
                       aes(label = gene),size = 2.5,box.padding = unit(0.25, "lines"),point.padding = unit(1e-05, "lines"),
                       label.padding = 0.3,segment.color = "black", show.legend = FALSE )




