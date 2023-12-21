rm(list = ls())
library(dplyr, verbose = F)
library(Seurat, verbose = F)
library(patchwork,verbose = F)
library(ggplot2, verbose = F)
library(ggpubr, verbose = F) 
library(ggthemes, verbose = F)
library(harmony,verbose = F)
library(stringr)
library(ggplot2)

#setwd("D:/project/scRNAseq/XiaoCaixia_Cellranger_23.1.29/")
set.seed(123)
## ======== 1. load the dataset and create seurat object ====
dir = list.files("./PANDO/sample/")
sample <- sapply(dir,function(x) strsplit(x,"[.]")[[1]][1])

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HHBG1" ,"HBG2","HBM","HBQ1","HBZ") 

scRNAlist<-list()
for(i in 1:length(dir))
{
  counts <- Read10X(data.dir = paste0("./PANDO/sample/",dir[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts,project = sample[i],min.cells=3)
  
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample[i])
  
  doublet<-read.table(paste0("./PANDO/",sample[i],"_doubletfinder_result.txt"),header = T)
  scRNAlist[[i]][["DoubletFinder"]] <- doublet$DoubletFinder
  
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  
  scRNAlist[[i]][["percent.rp"]] <- PercentageFeatureSet(scRNAlist[[i]],pattern = "^RP[SL]")
  
  HB.genes <- CaseMatch(HB.genes,rownames(scRNAlist[[i]]))
  scRNAlist[[i]][["percent.HB"]] <- PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
  
}

names(scRNAlist) <- sample
scRNA <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])

scRNA$sample<-ifelse(scRNA$orig.ident %in% c("P1","P5","S1"),"late-inflam","early-inflam")
scRNA$orig.ident<-factor(scRNA$orig.ident,levels = c("P3","P4","P1","P5","S1"))

scRNA = NormalizeData(scRNA)
scRNA = CellCycleScoring(scRNA,s.features=cc.genes.updated.2019$s.genes,g2m.features=cc.genes.updated.2019$g2m.genes)

scRNA <- subset(scRNA,subset = (DoubletFinder=="Singlet" & percent.mt<40 & percent.HB<10 
                                      & nFeature_RNA > 200  & nFeature_RNA < 4000 ))

mt.genes<-grep(pattern = "^MT-",rownames(scRNA),value=T) 
Rp.genes<-grep(pattern = "^RP[SL]",rownames(scRNA),value=T) # ^Rp[sl][[:digit:]]?
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HHBG1" ,"HBG2","HBM","HBQ1","HBZ") 
phase.genes = unlist(cc.genes.updated.2019)

scRNA<-scRNA[!rownames(scRNA) %in% c(mt.genes,HB.genes,Rp.genes,phase.genes),] #
dim(scRNA)
#[1] 22242 25791

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),npcs = 50)
scRNA <- RunHarmony(scRNA,group.by.vars = "orig.ident",plot_convergence = T)

hms<-38
scRNA <- FindNeighbors(scRNA, dims = 1:hms,reduction = "harmony")

scRNA <- FindClusters(scRNA, resolution = 0.2)
scRNA <- RunUMAP(scRNA, dims = 1:hms,reduction = "harmony")

#
new.cluster.ids<-c("T cells",  # 0
                   "Fibroblasts",   # 1
                   "T cells",  #2
                   "B cells",  #3
                   "Neutrophils",  # 4
                   "Endothelial cells",  #5
                   "Smooth muscle cells",#6
                   "Monocytes", #7
                   "Plasma cells", #8
                   "Epithelial cells", #9
                   "Mast cells", #10
                   "Plasmacytoid dendritic cells" #11
                   )  
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
levels(scRNA)<-c("T cells","B cells","Plasma cells","Monocytes","Neutrophils","Plasmacytoid dendritic cells",
                 "Mast cells","Smooth muscle cells","Fibroblasts","Endothelial cells","Epithelial cells")
scRNA$celltype<-Idents(scRNA)
saveRDS(scRNA, file = "./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.rds")

scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(scRNA.markers,file="./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.marker.csv")


###
mycol<-c("#58C9B9","#30A9DE","#476D87","#A593E0","#F6B352","#D499B9","#E3E36A","#DE6449","#5CAB7D","#D09E88","#7C7877")
scRNA<-readRDS("./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.rds")

p<-DimPlot(scRNA, reduction = "umap", group.by = "celltype", label = T,label.size = 4,pt.size = 0.1)
p+scale_color_manual(values = mycol )+theme_few()+NoLegend()+ggtitle("")

p<-DimPlot(scRNA, reduction = "umap", group.by = "celltype",split.by = "sample", label = F,label.size = 1,pt.size = 0.1) +
  scale_color_manual(values = mycol )
p+theme_few()+NoLegend()+ggtitle("")

p<-DimPlot(scRNA, reduction = "umap", group.by = "celltype",split.by = "orig.ident", label = F,label.size = 1,pt.size = 0.1) +
  scale_color_manual(values = mycol )
p+theme_void()+NoLegend()+ggtitle("")


##---density
mydata<-FetchData(scRNA,vars = c("sample","celltype","UMAP_1","UMAP_2"))
ggplot(mydata,aes(UMAP_1,UMAP_2)) +
  stat_density_2d(aes(fill = ..ndensity..),geom="raster",contour = FALSE,h=c(1,1))+
  scale_fill_gradientn(values=c(0,0.2,0.5,0.8,1),colours = c("white","#AFEEEE","#FFFFCC","#FD8D3C","#E31A1C"))+
  facet_wrap(~sample) + theme_void() + NoLegend()+
  geom_segment(aes(x=min(mydata$UMAP_1)-1,y=min(mydata$UMAP_2)-1,xend=min(mydata$UMAP_1)+8,yend=min(mydata$UMAP_2)-1),arrow = arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x=min(mydata$UMAP_1)-1,y=min(mydata$UMAP_2)-1,xend=min(mydata$UMAP_1)-1,yend=min(mydata$UMAP_1)+8),arrow = arrow(length = unit(0.1,"cm")))+
  annotate("text",x=min(mydata$UMAP_1)+3,y=min(mydata$UMAP_2)-2,label="UMAP_1",color="grey20",size=3,fontface="bold")+
  annotate("text",x=min(mydata$UMAP_1)-2.5,y=min(mydata$UMAP_2)+3,label="UMAP_2",color="grey20",size=3,fontface="bold",angle=90)


##---propotion of cell number
cell.prop<-as.data.frame(prop.table(table(scRNA@meta.data$sample,scRNA@meta.data$celltype),1))
colnames(cell.prop)<-c("sample","cluster","proportion")
p1<-ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values = mycol) +
  ggtitle("")+
  theme_bw()+
  theme(axis.title = element_text(size = 12),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.text = element_text(size = 11),legend.text = element_text(size = 10))+
  guides(fill=guide_legend(title=NULL))

#
prop<-prop.table(table(scRNA@meta.data$sample,scRNA@meta.data$celltype),1)
#
fold.prop<-as.data.frame(log2(prop[2,]/prop[1,]))
colnames(fold.prop)<-"fold"
fold.prop$cluster<-rownames(fold.prop)
fold.prop$cluster<-factor(fold.prop$cluster,levels = rev(fold.prop$cluster))
#fold.prop<-fold.prop[order(fold.prop$fold,decreasing = T),]
#fold.prop$cluster<-factor(fold.prop$cluster,levels = fold.prop$cluster)
p2<-ggplot(fold.prop,aes(cluster,fold,fill=cluster))+
  coord_flip()+ 
  geom_bar(stat="identity")+
  scale_fill_manual(values = rev(mycol)) +scale_x_discrete(position="top")+
  ggtitle("")+
  theme_bw()+
  theme(axis.text = element_text(size = 11,family = "B"),axis.title = element_text(size = 12))+
  guides(fill="none") #+labs(y="log2 fold change (inflam/early_inflam)")

##
ggplot(cell.prop[cell.prop$cluster %in% c("T cells","B cells"),],aes(cluster,proportion,fill=sample))+
  geom_bar(stat="identity",position=position_dodge2(padding=0.3))+
  geom_text(aes(label=round(proportion*100,digits = 2)),size=4,vjust = -0.3, position = position_dodge(.9))+
  scale_fill_manual(values = mycol) +
  ggtitle("")+
  theme_pubr()+
  theme(axis.title = element_text(size = 12),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.text = element_text(size = 11),legend.text = element_text(size = 10),legend.position ="right")+
  guides(fill=guide_legend(title=NULL))+ylab("percent of cell number(%)")

##
ggplot(cell.prop[cell.prop$cluster %in% c("Mast cells","Smooth muscle cells","Fibroblasts","Endothelial cells","Epithelial cells"),],aes(cluster,proportion,fill=sample))+
  geom_bar(stat="identity",position=position_dodge2(padding=0.3))+
  geom_text(aes(label=round(proportion*100,digits = 2)),size=4,vjust = -0.2, position = position_dodge(.9))+
  scale_fill_manual(values = mycol) +
  ggtitle("")+
  theme_pubr()+
  theme(axis.title = element_text(size = 12),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.text = element_text(size = 11),legend.text = element_text(size = 10),legend.position ="right")+
  guides(fill=guide_legend(title=NULL))+ylab("percent of cell number(%)")

##--mark genes
Tcell<-c("CD3D","CD3E") #,"PTPRC"
Bcell<-c("CD19","MS4A1","CD79A")
plasma<-c("IGHG1","CD79A") #"MZB1","SDC1",
Mon<-c("CSF1R","LYZ")
Np<-c("CSF1R","LYZ","CSF3R","FCGR3B")
Mast<-c("CPA3","KIT")
smc<-c("MYH11", "ACTA2")
Fb<-c("FGF7","COL1A1", "PDGFRA") #"MME",
Endc<-c("PECAM1", "ACKR1")
Epc<-c("KRT5", "KRT8","EPCAM")
pDC<-c("GZMB", "IL3RA")

genes_to_check = unique(c(Tcell,Bcell,plasma,Mon,Np,pDC,Mast,smc,Fb,Endc,Epc))

pdf("./PANDO/Total_cells_markgene_dotplot.pdf",width = 4,height = 5)
DotPlot(scRNA,features = genes_to_check, assay='RNA',dot.scale = 4,scale = T) + coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(),axis.text.y= element_text(size = 7), axis.text.x=element_text(size = 8,hjust = 1,vjust=1,angle = 45),
        legend.title = element_text(size=7),legend.text = element_text(size = 6))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradient2(high = "#34a186",mid = 'white',low = '#b5182b')
dev.off()

###
