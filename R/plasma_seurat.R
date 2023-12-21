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

srat_obj<-readRDS("./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.rds")
###-----
plasmacells<-subset(srat_obj,idents =c("Plasma cells"))

DimPlot(plasmacells, reduction = "umap",group.by = "celltype",label = T,label.size = 3)

plasmacells <- NormalizeData(plasmacells, normalization.method = "LogNormalize", scale.factor = 10000)
plasmacells <- ScaleData(plasmacells, features = rownames(plasmacells))
plasmacells <- FindVariableFeatures(plasmacells, selection.method = "vst", nfeatures = 2000)
plasmacells <- RunPCA(plasmacells, features = VariableFeatures(object = plasmacells),npcs = 50)
plasmacells <- RunHarmony(plasmacells,group.by.vars = "orig.ident",plot_convergence = F) #

reduce="harmony"
pcs<-run_PCnum(plasmacells,reduction = reduce)

plasmacells <- FindNeighbors(plasmacells, dims = 1:pcs,reduction = reduce)
plasmacells <- FindClusters(plasmacells, resolution =0.6)
plasmacells <- RunUMAP(plasmacells, dims = 1:pcs,reduction = reduce)

####---注释
#plasmacells<-readRDS("D:/project/scRNAseq/XiaoCaixia_Cellranger_23.1.29/harmony/plasmacells/plasmacells_harmony_res0.6.rds")
new.cluster.ids<-c("IgA2+ plasma",  # 0
                   "IgG1+ plasma"   # 1
                   )  

names(new.cluster.ids) <- levels(plasmacells)
plasmacells <- RenameIdents(plasmacells, new.cluster.ids)
levels(plasmacells)<-c("IgA2+ plasma","IgG1+ plasma")
plasmacells$plasmacell<-Idents(plasmacells)

saveRDS(plasmacells, file = "./PANDO/harmony/plasmacells/plasmacells_harmony_res0.6.annotaion.rds")

plasmacells.markers <- FindAllMarkers(plasmacells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(plasmacells.markers,file="./PANDO/harmony/plasmacells/plasmacells_harmony_res0.6.marker.csv")

####---
plasmacells<-readRDS("./PANDO/harmony/plasmacells/plasmacells_harmony_res0.6.annotaion.rds")
mycolor<-c("#58C9B9","#A593E0") #"#30A9DE",'#6778AE',,"#D499B9"

pdf("./PANDO/harmony/plasmacells/plasmacells_umap.pdf",width = 4,height = 2)
DimPlot(plasmacells, reduction = "umap", group.by = "plasmacell",split.by = "sample", label = F,label.size = 1,pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),axis.title = element_text(size = 9),
        legend.title = element_text(size=8),legend.text = element_text(size = 8),
        title = element_text(size = 9),strip.text = element_text(size = 8)) +
  scale_color_manual(values = mycolor )
dev.off()

##--marker genes
marker_genes<-c("JUN", "MALAT1","HSP90AA1","HSPA1B","HSPA1A", "NFKBIA", "KLF6","FOS" , "AHNAK", "IGHM","IGHA2", 
                "IGKC", "IGHG1","MSN")

pdf("D:/project/scRNAseq/XiaoCaixia_Cellranger_23.1.29/20230814/plasmacells_makergene_dotplot.pdf",width = 2.5,height = 3)
DotPlot(plasmacells,features =marker_genes,dot.scale = 4,scale = T,group.by = "plasmacell")+coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(),axis.text.y= element_text(size = 7), axis.text.x=element_text(size = 8,hjust = 1,vjust=1,angle = 45),
        legend.title = element_text(size=7),legend.text = element_text(size = 6),
        legend.key.height = unit(12, "pt"),legend.key.width = unit(10, "pt"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradient2(high = "#34a186",mid = 'white',low = '#9FA3A8')
dev.off()


## cell proportion
cell.prop<-as.data.frame(prop.table(table(plasmacells@meta.data$sample,plasmacells@meta.data$plasmacell),1))
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



