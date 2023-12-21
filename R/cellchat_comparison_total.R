library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
library(ggplot2)
library(Seurat)

options(stringsAsFactors = FALSE) 
options(futrue.globlas.Maxsize=2*1024**3)

setwd("./PANDO/harmony/CellChat/")

cellchat_inflam<-readRDS("./cellchat_late_inflam.rds")
cellchat_Einflam<-readRDS("./cellchat_early_inflam.rds")

object.list <- list(inflam = cellchat_inflam,early_inflam = cellchat_Einflam)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# To better control the node size and edge weights of the inferred networks across different datasets, 
# we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets. 
weight.max <- getMaxWeight(object.list, attribute = c("idents","count")) #count
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("number of interactions - ", names(object.list)[i]))
}

## Identify the upgulated and down-regulated signaling ligand-receptor pairs
# Identify dysfunctional signaling by using differential expression analysis
# perform differential expression analysis
pos.dataset = "late-inflam"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name,
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
write.csv(net,"./cellchat_netMappingDEG_inflam_vs_early_inflam.csv",row.names = F)

# extract the ligand-receptor pairs with upregulated ligands in KO
net.up <- subsetCommunication(cellchat, net = net, datasets = "inflam",ligand.logFC = 0.1, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in WT, i.e.,downregulated in KO
net.down <- subsetCommunication(cellchat, net = net, datasets = "early_inflam",ligand.logFC = -0.1, receptor.logFC = -0.1)

#Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
par(mfrow = c(1,1), xpd=TRUE)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("up-regulated signaling in ", names(object.list)[1]))
netVisual_chord_gene(object.list[[2]], slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("down-regulated signaling in ", names(object.list)[1]))

saveRDS(cellchat, file = "./CellChat/cellchat_inflam_vs_earyly_inflam.rds")


