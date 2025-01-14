mat <- sweep(outgoing_combined_reorder, 1L, apply(outgoing_combined_reorder, 1, max), '/', check.margin = FALSE)
color.heatmap.use =colorRampPalette((brewer.pal(n = 9, name = "BuGn")))(100)

reordered_mat =mat[,c("alpha_ctrl","alpha_PA",
                      "beta_ctrl","beta_PA",
                      "delta_ctrl","delta_PA",
                      "pp_ctrl","pp_PA")]

Heatmap(reordered_mat,cluster_columns = FALSE,col = color.heatmap.use,
        name = "Relative strength")


indeg_control <- lapply(object.list$Ctrl@netP$centr, function(signaling) {
  unlist(signaling$indeg)
})

indeg_control_df <- do.call(rbind,indeg_control)

indeg_PA <- lapply(object.list$PA@netP$centr, function(signaling) {
  unlist(signaling$indeg)
})

indeg_PA_df <- do.call(rbind,indeg_PA)

colnames(indeg_control_df)=c("alpha_ctrl","beta_ctrl",
                             "delta_ctrl","pp_ctrl")
colnames(indeg_PA_df)=c("alpha_PA","beta_PA",
                             "delta_PA","pp_PA")

indeg_combined_df <- merge(indeg_control_df, indeg_PA_df, by = "row.names", all = TRUE)

indeg_combined_df[is.na(indeg_combined_df)] <- 0

rownames(indeg_combined_df)=indeg_combined_df$Row.names

indeg_combined_df=indeg_combined_df[,-1]

indeg_combined_df_reorder=indeg_combined_df[,c("alpha_ctrl","alpha_PA",
                                                       "beta_ctrl","beta_PA",
                                                       "delta_ctrl","delta_PA",
                                                       "pp_ctrl","pp_PA")]

indeg_mat <- sweep(indeg_combined_df_reorder, 1L, apply(indeg_combined_df_reorder, 1, max), '/', check.margin = FALSE)

color.heatmap.use =colorRampPalette((brewer.pal(n = 9, name = "GnBu")))(100)

Heatmap(indeg_mat,cluster_columns = FALSE,col = color.heatmap.use,
        name = "Relative strength")

all_control_df=outdeg_ctrl_df+indeg_ctrl_df
all_PA_df=outdeg_PA_df+indeg_PA_df

all_combined_df <- merge(all_control_df, all_PA_df, by = "row.names", all = TRUE)

all_combined_df[is.na(all_combined_df)] <- 0

rownames(all_combined_df)=all_combined_df$Row.names

all_combined_df=all_combined_df[,-1]

all_combined_df_reorder=all_combined_df[,c("alpha_ctrl","alpha_PA",
                                               "beta_ctrl","beta_PA",
                                               "delta_ctrl","delta_PA",
                                               "pp_ctrl","pp_PA")]

all_mat <- sweep(all_combined_df_reorder, 1L, apply(all_combined_df_reorder, 1, max), '/', check.margin = FALSE)

color.heatmap.use =colorRampPalette((brewer.pal(n = 9, name = "OrRd")))(100)

Heatmap(all_mat,cluster_columns = FALSE,col = color.heatmap.use,
        name = "Relative strength")

gg1 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)

gg1 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", sources.use = c("alpha","beta","delta","pp"), targets.use = c("alpha","beta","delta","pp"), stacked = T, do.stat = TRUE)

gg2 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "alpha", targets.use =c("alpha","beta","delta","pp"), stacked = T, do.stat = TRUE)
gg3 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "beta", targets.use =c("alpha","beta","delta","pp"), stacked = T, do.stat = TRUE)
gg4 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "delta", targets.use =c("alpha","beta","delta","pp"), stacked = T, do.stat = TRUE)
gg5 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "pp", targets.use =c("alpha","beta","delta","pp"), stacked = T, do.stat = TRUE)

library(gridExtra)

grid.arrange(gg2, gg3, gg4, gg5, ncol = 4)

gg12 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", sources.use = "alpha", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE)
gg13 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", sources.use = "beta", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE)
gg14 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", sources.use = "delta", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE)
gg15 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", sources.use = "pp", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE)
gg16 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", sources.use = "epsilon", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE)

library(gridExtra)

grid.arrange(gg12, gg13, gg14, gg15,gg16, nrow=1)

pathways.show <- c("EPHA","Netrin","SEMA6","TGFb","CLDN","NRXN","APP","GRN","WNT","JAM","PTPR","ADGRL","Cholesterol","DESMOSOME",
                   "CADM","PTPRM","LAMININ","GIPR","AGT")



pathways.show <- c("AGT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("TGFb") 
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(object.list$Ctrl, signaling = pathways.show, layout = "circle")

nPatterns = 3
lipo.object.list.all$Ctrl <- identifyCommunicationPatterns(lipo.object.list.all$Ctrl, pattern = "outgoing", k = nPatterns)

netAnalysis_river(object.list$Ctrl, pattern = "incoming")

netAnalysis_dot(object.list$Ctrl, pattern = "incoming")

selectK(object.list$PA, pattern = "outgoing")
lipo.object.list.all$PA <- identifyCommunicationPatterns(lipo.object.list.all$PA, pattern = "outgoing", k = nPatterns)
netAnalysis_river(object.list$PA, pattern = "outgoing")
netAnalysis_dot(object.list$PA, pattern = "outgoing")


lipo.object.list.all$PA <- identifyCommunicationPatterns(lipo.object.list.all$PA, pattern = "incoming", k = nPatterns)
lipo.object.list.all$Ctrl <- identifyCommunicationPatterns(lipo.object.list.all$Ctrl, pattern = "incoming", k = nPatterns)

unique_in_PA <- setdiff(colnames(lipo.object.list.all$PA@netP$pattern$outgoing$data),
                        colnames(lipo.object.list.all$Ctrl@netP$pattern$outgoing$data))

# Find unique column names in the second vector compared to the first vector
unique_in_Ctrl <- setdiff(colnames(lipo.object.list.all$Ctrl@netP$pattern$outgoing$data),
                          colnames(lipo.object.list.all$PA@netP$pattern$outgoing$data))

unique_in_PA_in <- setdiff(colnames(lipo.object.list.all$PA@netP$pattern$incoming$data),
                        colnames(lipo.object.list.all$Ctrl@netP$pattern$incoming$data))

# Find unique column names in the second vector compared to the first vector
unique_in_Ctrl_in <- setdiff(colnames(lipo.object.list.all$Ctrl@netP$pattern$incoming$data),
                          colnames(lipo.object.list.all$PA@netP$pattern$incoming$data))


signaling_pattern_list=c("CADM","OCLN","MPZ","PDGF","Netrin","UNC5","SEMA6","EPHA",
                        "Desmosome","ADGRA","EDN","CD45","BMP","PECAM1","MPZ",
                        "ESAM","ANGPT","TWEAK","Glutamate","CDH5","GAS",
                        "EPHB","IL1")

dif_ctrl_incoming=lipo.object.list.all$Ctrl@netP$pattern$incoming$pattern$signaling[lipo.object.list.all$Ctrl@netP$pattern$incoming$pattern$signaling$Signaling %in% signaling_pattern_list,]
dif_PA_incoming=lipo.object.list.all$PA@netP$pattern$incoming$pattern$signaling[lipo.object.list.all$PA@netP$pattern$incoming$pattern$signaling$Signaling %in% signaling_pattern_list,]

dif_ctrl_outgoing=lipo.object.list.all$Ctrl@netP$pattern$outgoing$pattern$signaling[lipo.object.list.all$Ctrl@netP$pattern$outgoing$pattern$signaling$Signaling %in% signaling_pattern_list,]
dif_PA_outgoing=lipo.object.list.all$PA@netP$pattern$outgoing$pattern$signaling[lipo.object.list.all$PA@netP$pattern$outgoing$pattern$signaling$Signaling %in% signaling_pattern_list,]

gg2 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = c("alpha","beta","delta","pp"), targets.use ="alpha", stacked = T, do.stat = TRUE)
gg3 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = c("alpha","beta","delta","pp"), targets.use ="beta", stacked = T, do.stat = TRUE)
gg4 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = c("alpha","beta","delta","pp"), targets.use ="delta", stacked = T, do.stat = TRUE)
gg5 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = c("alpha","beta","delta","pp"), targets.use ="pp", stacked = T, do.stat = TRUE)

grid.arrange(gg2, gg3, gg4, gg5, ncol = 4)

lipo.cellchat.all@meta$datasets = factor(lipo.cellchat.all@meta$datasets, levels = c("Ctrl", "PA"))
plotGeneExpression(lipo.cellchat.all, signaling = "NPY", split.by = "datasets", colors.ggplot = T, type = "violin")


gg12 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "count", sources.use = "alpha", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE,title="alpha")
gg13 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "count", sources.use = "beta", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE,title="beta")
gg14 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "count", sources.use = "delta", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE,title = "delta")
gg15 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "count", sources.use = "pp", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE, title="pp")
gg16 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "count", sources.use = "epsilon", targets.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE,title = "epsilon")

library(gridExtra)
grid.arrange(gg12, gg13, gg14,gg16,gg15, nrow=1)

gg121 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", targets.use = "alpha", sources.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE,title="alpha")
gg131 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", targets.use = "beta", sources.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE,title="beta")
gg141 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", targets.use = "delta", sources.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE,title = "delta")
gg151 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", targets.use = "pp", sources.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE, title="pp")
gg161 <- rankNet(lipo.cellchat.all, mode = "comparison", measure = "weight", targets.use = "epsilon", sources.use =c("alpha","beta","delta","epsilon","pp"), stacked = T, do.stat = TRUE,title = "epsilon")

grid.arrange(gg121, gg131, gg141,gg161,gg151, nrow=1)

load("/Users/liguo/Desktop/lipoglucotoxicity/all_samples/RNA_cellchat/abcde_only_object.list.RData")
load("/Users/liguo/Desktop/lipoglucotoxicity/all_samples/RNA_cellchat/abcde_only.RData")

gg2 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "alpha", targets.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="alpha")
gg3 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "beta", targets.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="beta")
gg4 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "delta", targets.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="delta")
gg5 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "pp", targets.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="pp")
gg6 <-rankNet(merged.cellchat, mode = "comparison", measure = "weight", sources.use = "epsilon", targets.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="epsilon")

grid.arrange(gg2, gg3, gg4, gg5,gg6, nrow=1)

gg22 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", targets.use = "alpha", sources.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="alpha")
gg32 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", targets.use = "beta", sources.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="beta")
gg42 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", targets.use = "delta", sources.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="delta")
gg52 <- rankNet(merged.cellchat, mode = "comparison", measure = "weight", targets.use = "pp", sources.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="pp")
gg62 <-rankNet(merged.cellchat, mode = "comparison", measure = "weight", targets.use = "epsilon", sources.use =c("alpha","beta","delta","pp","epsilon"), stacked = T, do.stat = TRUE,title="epsilon")

grid.arrange(gg22, gg32, gg42, gg52,gg62, nrow=1)



indeg_control <- lapply(lipo.object.list.all$Ctrl@netP$centr, function(signaling) {
  unlist(signaling$indeg)
})

indeg_control_df <- do.call(rbind,indeg_control)

indeg_control_df_endocrine=indeg_control_df[,c("alpha","beta","delta","epsilon","pp")]

indeg_PA <- lapply(lipo.object.list.all$PA@netP$centr, function(signaling) {
  unlist(signaling$indeg)
})

indeg_PA_df <- do.call(rbind,indeg_PA)

indeg_PA_df_endocrine=indeg_PA_df[,c("alpha","beta","delta","epsilon","pp")]


colnames(indeg_control_df_endocrine)=c("alpha_ctrl","beta_ctrl",
                             "delta_ctrl","epsilon_ctrl","pp_ctrl")
colnames(indeg_PA_df_endocrine)=c("alpha_GL","beta_GL",
                        "delta_GL","epsilon_GL","pp_GL")

indeg_combined_df <- merge(indeg_control_df_endocrine, indeg_PA_df_endocrine, by = "row.names", all = TRUE)

indeg_combined_df[is.na(indeg_combined_df)] <- 0

rownames(indeg_combined_df)=indeg_combined_df$Row.names

indeg_combined_df=indeg_combined_df[,-1]

indeg_combined_df_reorder=indeg_combined_df[,c("alpha_ctrl","alpha_GL",
                                               "beta_ctrl","beta_GL",
                                               "delta_ctrl","delta_GL",
                                               "epsilon_ctrl","epsilon_GL",
                                               "pp_ctrl","pp_GL")]

indeg_combined_df_reorder_filtered <- indeg_combined_df_reorder[rowSums(indeg_combined_df_reorder) != 0, ]


indeg_mat <- sweep(indeg_combined_df_reorder_filtered, 1L, apply(indeg_combined_df_reorder_filtered, 1, max), '/', check.margin = FALSE)

color.heatmap.use =colorRampPalette((brewer.pal(n = 9, name = "GnBu")))(100)

Heatmap(indeg_mat,cluster_columns = FALSE,col = color.heatmap.use,
        name = "Relative strength")


outdeg_control <- lapply(lipo.object.list.all$Ctrl@netP$centr, function(signaling) {
  unlist(signaling$outdeg)
})

outdeg_control_df <- do.call(rbind,outdeg_control)

outdeg_control_df_endocrine=outdeg_control_df[,c("alpha","beta","delta","epsilon","pp")]

outdeg_PA <- lapply(lipo.object.list.all$PA@netP$centr, function(signaling) {
  unlist(signaling$outdeg)
})

outdeg_PA_df <- do.call(rbind,outdeg_PA)

outdeg_PA_df_endocrine=outdeg_PA_df[,c("alpha","beta","delta","epsilon","pp")]


colnames(outdeg_control_df_endocrine)=c("alpha_ctrl","beta_ctrl",
                                       "delta_ctrl","epsilon_ctrl","pp_ctrl")
colnames(outdeg_PA_df_endocrine)=c("alpha_GL","beta_GL",
                                  "delta_GL","epsilon_GL","pp_GL")

outdeg_combined_df <- merge(outdeg_control_df_endocrine, outdeg_PA_df_endocrine, by = "row.names", all = TRUE)

outdeg_combined_df[is.na(outdeg_combined_df)] <- 0

rownames(outdeg_combined_df)=outdeg_combined_df$Row.names

outdeg_combined_df=outdeg_combined_df[,-1]

outdeg_combined_df_reorder=outdeg_combined_df[,c("alpha_ctrl","alpha_GL",
                                               "beta_ctrl","beta_GL",
                                               "delta_ctrl","delta_GL",
                                               "epsilon_ctrl","epsilon_GL",
                                               "pp_ctrl","pp_GL")]

outdeg_combined_df_reorder_filtered <- outdeg_combined_df_reorder[rowSums(outdeg_combined_df_reorder) != 0, ]


outdeg_mat <- sweep(outdeg_combined_df_reorder_filtered, 1L, apply(outdeg_combined_df_reorder_filtered, 1, max), '/', check.margin = FALSE)

color.heatmap.use =colorRampPalette((brewer.pal(n = 9, name = "BuGn")))(100)

Heatmap(outdeg_mat,cluster_columns = FALSE,col = color.heatmap.use,
        name = "Relative strength")

plotGeneExpression(lipo.cellchat.all, signaling = "INSULIN", split.by = "datasets", colors.ggplot = T, type = "violin")
plotGeneExpression(lipo.cellchat.all, signaling = "GCG", split.by = "datasets", colors.ggplot = T, type = "violin")
plotGeneExpression(lipo.cellchat.all, signaling = "SOMATOSTATIN", split.by = "datasets", colors.ggplot = T, type = "violin")
plotGeneExpression(lipo.cellchat.all, signaling = "NPY", split.by = "datasets", colors.ggplot = T, type = "violin")
plotGeneExpression(lipo.cellchat.all, signaling = "GHRELIN", split.by = "datasets", colors.ggplot = T, type = "violin")

lipo$cell.type.final<-factor(lipo$cell.type.final,levels=c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune","doublets","unknown"))
Idents(lipo) = "cell.type.final"

VlnPlot(lipo, features = c("GLP1R","INSR", "SSTR2","GHSR","NPY1R"),idents = c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune"), split.by = "condition",ncol=2)
VlnPlot(lipo, features = c("GLP1R","INSR", "SSTR2","GHSR","NPY1R"),idents = c("alpha","beta","delta","epsilon","pp"), split.by = "condition",ncol=2)

DefaultAssay(lipo) <- "RNA"
VlnPlot(lipo, features = c("GLP1R","INSR", "SSTR2","GHSR","NPY1R"),idents = c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune"), split.by = "condition",ncol=5)
VlnPlot(lipo, features = c("GLP1R","INSR", "SSTR2","GHSR","NPY1R"),idents = c("alpha","beta","delta","epsilon","pp"), split.by = "condition",ncol=1)

VlnPlot(lipo, features = c("GCGR","GLP1R", "PPYR1","NPY4R","NPY1R"),idents = c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune"), split.by = "condition",ncol=2)

VlnPlot(lipo, features = c("SSTR1","SSTR2", "SSTR3","SSTR4","SSTR5"),idents = c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune"), split.by = "condition",ncol=4)

VlnPlot(lipo, features = c("GCGR","GLP1R", "GLP2R"),idents = c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune"), split.by = "condition",ncol=3)

VlnPlot(lipo, features = c("NPY4-R","PPYR1", "PP1","Y4","NPY1R"),idents = c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune"), split.by = "condition",ncol=2)

lipo$cell.condition=paste(lipo$cell.type.final,lipo$condition,sep="_")
lipo$cell.condition=factor(lipo$cell.condition,levels=c(
  "alpha_Ctrl","alpha_PA","beta_Ctrl","beta_PA","delta_Ctrl","delta_PA",
  "epsilon_Ctrl","epsilon_PA","pp_Ctrl","pp_PA","ductal_Ctrl","ductal_PA",
  "acinar_Ctrl","acinar_PA","fibroblast_Ctrl","fibroblast_PA",
  "endothelial_Ctrl","endothelial_PA","immune_Ctrl","immune_PA","doublets_Ctrl","doublets_PA",
  "unknown_Ctrl","unknown_PA"
))

Idents(lipo) = "cell.condition"

DoHeatmap(lipo,
          features = c("GCGR","GLP1R","INSR","SSTR1","SSTR2","SSTR3","SSTR5",
                       "GHSR","NPY1R"),slot="data")

DoHeatmap(lipo,
          features = c("GCGR"))

receptor_RNA_data=lipo@assays$RNA$data[c("GCGR","GLP1R","INSR","SSTR1","SSTR2","SSTR3","SSTR5",
                                         "GHSR","NPY1R"),]

receptor_RNA_data_sum=rowSums(receptor_RNA_data)
