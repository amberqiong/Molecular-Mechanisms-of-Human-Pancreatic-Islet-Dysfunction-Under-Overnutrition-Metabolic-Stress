library(CellChat)
library(patchwork)
library(Seurat)

lipoglucotoxicity=readRDS("/gpfs/home/lg23w/lipoglucotoxicity/lipoglucotoxicity.integrated.120823.rds")
Idents(lipoglucotoxicity)<- "cell.type.final"
lipo.subset=subset(lipoglucotoxicity,
                   idents=c("acinar","alpha","beta","delta","ductal","epsilon",
                            "endothelial","fibroblast","immune","pp"))
rm(lipoglucotoxicity)

set.seed(123)
options(future.seed=TRUE)

lipo.data.input=lipo.subset[["RNA"]]$data
metadata=lipo.subset@meta.data
cell.use.PA=rownames(metadata)[metadata$condition=="PA"]

lipo.data.input_1=lipo.data.input[,cell.use.PA]
metadata_1=metadata[cell.use.PA,]
metadata_cellchat.PA=data.frame(labels=metadata_1$cell.type.final,row.names = rownames(metadata_1))


PA_cellchat= createCellChat(object=lipo.data.input_1,meta=metadata_cellchat.PA,group.by = "labels")



cell.use.Ctrl=rownames(metadata)[metadata$condition=="Ctrl"]
lipo.data.input_2=lipo.data.input[,cell.use.Ctrl]
metadata_2=metadata[cell.use.Ctrl,]
metadata_cellchat.Ctrl=data.frame(labels=metadata_2$cell.type.final,row.names = rownames(metadata_2))

Ctrl_cellchat= createCellChat(object=lipo.data.input_2,meta=metadata_cellchat.Ctrl,group.by = "labels")

cell.type.order=c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune")
Ctrl_cellchat=setIdent(Ctrl_cellchat,ident.use = "labels",levels =cell.type.order)
PA_cellchat=setIdent(PA_cellchat,ident.use = "labels",levels =cell.type.order)

rm(lipo.subset)
CellchatDB=CellChatDB.human

CellchatDB.use <- CellchatDB
Ctrl_cellchat@DB <- CellchatDB.use

Ctrl_cellchat<-subsetData(Ctrl_cellchat)
future::plan("multisession",workers=6)
Ctrl_cellchat <- identifyOverExpressedGenes(Ctrl_cellchat)
Ctrl_cellchat <- identifyOverExpressedInteractions(Ctrl_cellchat)
Ctrl_cellchat <-computeCommunProb(Ctrl_cellchat,type="triMean")
Ctrl_cellchat <-computeCommunProbPathway(Ctrl_cellchat)
Ctrl_cellchat <-aggregateNet(Ctrl_cellchat)
Ctrl_cellchat <- netAnalysis_computeCentrality(Ctrl_cellchat, slot.name = "netP")

#####PA
PA_cellchat@DB <- CellchatDB.use

PA_cellchat<-subsetData(PA_cellchat)
future::plan("multisession",workers=6)
PA_cellchat <- identifyOverExpressedGenes(PA_cellchat)
PA_cellchat <- identifyOverExpressedInteractions(PA_cellchat)

PA_cellchat <-computeCommunProb(PA_cellchat,type="triMean")
PA_cellchat <-computeCommunProbPathway(PA_cellchat)
PA_cellchat <-aggregateNet(PA_cellchat)
PA_cellchat <- netAnalysis_computeCentrality(PA_cellchat, slot.name = "netP")

object.list <- list(Ctrl=Ctrl_cellchat,PA=PA_cellchat)
merged.cellchat <- mergeCellChat(object.list,add.names = names(object.list))

save(object.list, file = "merged_cellchat_withoutdoubletsNunknown_object.list.RData")
save(merged.cellchat, file="merged_cellchat_withoutdoubletsNunknown.RData")

saveRDS(Ctrl_cellchat,file="Ctrl_cellchat_reorder.rds")
saveRDS(PA_cellchat,file="PA_cellchat_reorder.rds")
