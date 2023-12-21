library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
library(ggplot2)
library(Seurat)

options(stringsAsFactors = FALSE) 
options(futrue.globlas.Maxsize=2*1024**3)
suppressWarnings(suppressMessages(future::plan("multiprocess", workers = 8)))

srat_obj<-readRDS("./PANDO/harmony/scRNA_res0.2.annotaion.levels.merge.rds")
srat_sub<-subset(srat_obj,idents =c("T cells","Epithelial cells","Fibroblasts","B cells" ,"Neutrophils","Monocytes"))

E_inflam<-subset(srat_sub,subset=sample=="early-inflam")
L_inflam<-subset(srat_sub,subset=sample=="late-inflam")

###-----early-inflam
cellchat <- createCellChat(object = E_inflam)

#--database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat@DB <- CellChatDB.use

#
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = T)  
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

saveRDS(cellchat, file = "./PANDO/CellChat/cellchat_early_inflam.rds")

###-----late-inflam
cellchat <- createCellChat(object = L_inflam)

cellchat@DB <- CellChatDB.use

#
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = T)  
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

saveRDS(cellchat, file = "./PANDO/CellChat/cellchat_late_inflam.rds")

