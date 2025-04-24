# Molecular-Mechanisms-of-Human-Pancreatic-Islet-Dysfunction-Under-Overnutrition-Metabolic-Stress

This repository hosts the codes to generate analysis and figures for Molecular Mechanisms of Human Pancreatic Islet Dysfunction Under Overnutrition Metabolic stress, currently under review in Diabetes

## Pre-processing

Before making the actual analysis, begin with ambient mRNA decontamination by SoupX, Modified from https://github.com/Gaulton-Lab/HPAP-scRNA-seq/blob/main/HPAP-SoupX.R

```r
library(Seurat)
library(dplyr)
library(SoupX)
library(Azimuth)
library(SeuratData)
library(ggplot2)
library(patchwork)

samples <- c('CITH107_Ctrl','CITH107_PA','HP22286_Ctrl','HP22286_PA', 'HP23098_Julia_Ctrl',
              'HP23098_Julia_PA','HP23098_Snow_Ctrl','HP23098_Snow_PA','HP23166_Ctrl','HP23166_PA',
              'SAMN32641505_Ctrl','SAMN32641505_PA','SAMN36705973_Ctrl','SAMN36705973_PA')

sc_before_soupx_list <- list() # save pre SoupX seurat objects
sc_after_soupx_list <- list() # save after SoupX seurat objects

for (sample in samples){
  # Load seurat objects, add sample predix to barcode 
  sc_object <- readRDS(paste0(sample,'.rds'))
  sc_object <- RenameCells(object=sc_object,add.cell.id=sample)
  sc_before_soupx_list[[sample]] <-sc_object

  DefaultAssay(sc_object) <- "RNA"

  # Get raw counts (raw_feature_bc_matrix.h5) and filtered counts (from seurat)
  toc <-GetAssayData(object=sc_object,assay="RNA",slot="counts")
  tod <-Seurat::Read10X_h5(file.path(sample,"raw_feature_bc_matrix.h5"))
  tod_1 <- tod[rownames(toc),] # common set

  # Get metadata (UMAP and clusters)
  metadata <- (cbind(as.data.frame(sc_object@reductions$umap@cell.embeddings),
                     as.data.frame(Idents(sc_object))
  ))
  colnames(metadata) <- c('RD1','RD2','Cluster')

  # Run SoupX
  sc <- SoupChannel(tod_1,toc)
  sc <- setDR(sc,metadata[colnames(sc$toc),c('RD1','RD2')])
  sc <- setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))
  sc <- autoEstCont(sc)

  # save ambient RNA contamination estimates
  contamination_data <- sc$fit$dd[,c('gene', 'soupExp')]
  colnames(contamination_data) <- c('gene', sample)
  write.table(contamination_data, paste0(sample, '_gene_level_soup_exp.tsv'), sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Stochastically adjusts counts while rounding to maintain overall contamination fraction while outputting integer counts 
  out <- adjustCounts(sc, roundToInt=TRUE)

  # Create new seurat object with corrected counts
  sc_object_new <- CreateSeuratObject(out)
  sc_object_new[['percent.mt']] <- PercentageFeatureSet(sc_object_new, pattern = '^MT-')
  sc_object_new <- SCTransform(sc_object_new, vst.flavor = "v2", verbose = FALSE) 
  sc_object_new <- RunPCA(sc_object_new, npcs=30) 
  sc_object_new <- RunUMAP(sc_object_new,reduction="pca",dims=1:30,verbose=FALSE) 
  sc_object_new <- FindNeighbors(sc_object_new,reduction = "pca", dims = 1:30, verbose = FALSE)
  sc_object_new <- FindClusters(sc_object_new,resolution=1)
  saveRDS(sc_object_new, file = paste0(sample, "_SoupX.rds"))

  sc_after_soupx_list[[sample]] <- sc_object_new

}

# Create a merged Seurat object lipo from the individual sample post-SoupX Seurat objects

lipo <- merge(
  sc_after_soupx_list[[1]],  # The first Seurat object (base for merging)
  y = sc_after_soupx_list[-1],  # All remaining Seurat objects
  add.cell.ids = names(sc_after_soupx_list), 
  project = "lipoglucotoxicity"
)

lipo$sample <- sub('_[^_]*$', '',rownames(lipo@meta.data))
lipo$donor <- sub('_[^_]*$', '',lipo$sample)
lipo$donor <- sub('_[^_]*$', '',lipo$donor)
lipo$condition <- sub('.*_', '',lipo$sample)
lipo$sex <- ifelse(lipo$donor%in%c("SAMN36705973","HP23166","HP23098"),"female","male")
lipo[['percent.mt']] <- PercentageFeatureSet(lipo, pattern = '^MT-')

```

## Filtering, integration and cell calling

Seurat tutorial is followed.

```R
library(Seurat)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)

lipo[["RNA"]] <- split(lipo[["RNA"]], f = lipo$sample)

lipo <- subset(lipo,subset=nFeature_RNA>200 & ncount_RNA <=10000 & percent.mt<=10)

lipo <- NormalizeData(lipo)
lipo <- FindVariableFeatures(lipo)
lipo <- ScaleData(lipo)
lipo <- RunPCA(lipo)
lipo <- FindNeighbors(lipo, dims = 1:30, reduction = "pca")
lipo <- FindClusters(lipo, resolution = 2, cluster.name = "unintegrated_clusters")
lipo <- RunUMAP(lipo, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(lipo, reduction = "umap.unintegrated", group.by = c("condition", "donor"))

# Comparison between different integration method

## cca
lipo <- IntegrateLayers(
  object = lipo, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
lipo <- FindNeighbors(lipo, reduction = "integrated.cca", dims = 1:30)
lipo <- FindClusters(lipo, resolution = 2, cluster.name = "cca_clusters")
lipo <- RunUMAP(lipo, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(
  lipo,
  reduction = "umap.cca",
  group.by = c("condition", "donor", "cca_clusters"),
  label.size = 2
)

## rpca

lipo <- IntegrateLayers(
  object = lipo, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
lipo <- FindNeighbors(lipo, reduction = "integrated.rpca", dims = 1:30)
lipo <- FindClusters(lipo, resolution = 2, cluster.name = "rpca_clusters")
lipo <- RunUMAP(lipo, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

DimPlot(
  lipo,
  reduction = "umap.rpca",
  group.by = c("condition", "donor", "rpca_clusters"),
  label.size = 2
)

## harmony
lipo <- IntegrateLayers(
  object = lipo, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

lipo <- FindNeighbors(lipo, reduction = "harmony", dims = 1:30)
lipo <- FindClusters(lipo, resolution = 2, cluster.name = "harmony_clusters")
lipo <- RunUMAP(lipo, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(
  lipo,
  reduction = "umap.harmony",
  group.by = c("condition", "donor", "harmony_clusters"),
  label.size = 2
)

## mnn

lipo <- IntegrateLayers(
  object = lipo, method = FastMNNIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

lipo <- FindNeighbors(lipo, reduction = "integrated.mnn", dims = 1:30)
lipo <- FindClusters(lipo, resolution = 2, cluster.name = "mnn_clusters")
lipo <- RunUMAP(lipo, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")

DimPlot(
  lipo,
  reduction = "umap.mnn",
  group.by = c("condition", "donor", "mnn_clusters"),
  label.size = 2
)

features=c("GCG", "INS", "SST", "GHRL","PPY","CFTR","CPA2","SPARC", "VWF","PTPRC","BRCA1")
reductions=c("umap.cca","umap.rpca","umap.harmony","umap.mnn")
DefaultAssay(lipo) = "SCT"

pdf("FeaturePlots.pdf")
for (reduction in reductions){
  plot <- FeaturePlot(lipo, features = features, cols = c("lightgrey", "red"), reduction=reduction)+
  plot_annotation(title=paste("Reduction",reduction))
  print(plot)
}
dev.off()

##cca method has the best integration results based on visualization from feature plots

# Next, cell type calling based on marker gene expressions 

VlnPlot(lipo,features=features,group_bym, group_by="cca_clusters",ncol=4)

new.cluster.ids <- c("ductal","ductal","alpha","ductal","acinar","fibroblast","alpha","beta",
                     "alpha","doublets","alpha","alpha","ductal","beta","beta","ductal",
                     "delta","alpha","pp","acinar","ductal","beta","alpha","alpha","acinar",
                     "acinar","unknown","pp","beta","unknown","fibroblast","beta","endothelial",
                     "alpha","delta","alpha","doublets","fibroblast","beta","ductal","immune",
                     "immune","endothelial","endothelial","immune"
)

names(new.cluster.ids) <- levels(lipo$cca_clusters)
Idents(lipo)=lipo$cca_clusters
lipo <- RenameIdents(lipo, new.cluster.ids)
lipo$cell.type <- Idents(lipo)

lipo$cca_clusters <- factor(lipo$cca_clusters,levels = c(0:44))

DimPlot(lipo, reduction = "umap.cca", label = TRUE, pt.size = 0.5) 
FeaturePlot(lipo, reduction="umap.cca",features = c("percent.mt","nCount_RNA", "nFeature_RNA"))
lipo <- JoinLayers(lipo)

## Epsilon population has too few cells to be acurately annotated. Use Azimuth reference mapping to enhance epsilon cell annotation. 
### Final cell type annotation is manual annotation updated with epsilon cell lable from azimuth, azimuth reference files could be downloaded to local from Zenodo from Azimuth Website: https://azimuth.hubmapconsortium.org/references/ or loaded from SeuratData

InstallData("pancreasref.SeuratData")
lipo <- RunAzimuth(lipo,reference="pancreasref")
DimPlot(lipo, reduction = "umap.cca", group.by="predicted.annotation.l1",label = TRUE, pt.size = 0.5)
lipo$cell.type <- as.character(lipo$cell.type)
lipo$cell.type.final <- ifelse(lipo$predicted.annotation.l1=="epsilon","epsilon",lipo$cell.type)
Idents(lipo) <- "cell.type.final"
table(Idents(lipo))
lipo <- subset(lipo,idents=c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune","doublets"))

## UMAP, color code for different cell types
colors <- c("#4682B4","#CD5C5C", "#5F9EA0", "firebrick","#87CEEB", "#FF8C00", "#48D1CC","#FFD700", "#7B68EE","#FF6347","darkgrey")
DimPlot(lipo, reduction = "umap.cca",cols = colors,label = TRUE)+labs(x="UMAP1",y="UMAP2")

save(lipo,file="lipo_integrated.rds")

## features dot plot
Idents(lipo) <- factor(x=Idents(lipo),levels=rev(c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune","doublets")))

features=rev(c("GCG","INS","SST","GHRL","PPY","CFTR","CPA2","SPARC","VWF","PTPRC"))
DotPlot(lipoglucotoxicity,features = features)

## cell percentile and plotting

counts_data=lipo@meta.data %>%
   group_by(donor,condition, cell.type.final) %>%
   summarize(count=n()) %>%
   ungroup() %>%
   filter(cell.type.final != "doublets")

ggplot(counts_data,aes(x=condition,y=count,fill=cell.type.final))+
  facet_grid(rows=vars(donor))+geom_bar(position = "fill",stat = "identity")+
  scale_fill_manual("cell.type.final",values = c("alpha"="#DE8C00",
  "beta"="#B79F00","delta"="#7CAE00","epsilon"="#00B4F0","pp"="#F564E3"))+
  theme_classic()+
  geom_text(aes(label=count),position=position_fill(vjust = 0.5),size=3)+coord_flip()

```

## DEG Analysis
DEG analysis for all cell types comparing the conditions `GL` vs `Ctrl`. The `FindMarkers` function from the Seurat package is used for differential expression analysis, refer to the [Seurat Differential Expression Vignette](https://satijalab.org/seurat/articles/de_vignette).

```r
library(Seurat)

lipo$celltype.condition=paste(lipo$cell.type.final,lipo$condition,sep="_")
Idents(lipo)="celltype.condition"
DefaultAssay(lipo) <- "RNA"

## Run analysis on all the cell types

cell_types=c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune")
DEG_lipo=list()
for (celltype in cell_types){
  ident_1 <- paste0(cell_type, "_GL")
  ident_2 <- paste0(cell_type, "_Ctrl")
  DEG_lipo[[celltype]] <- FindMarkers(lipo,assay="RNA",
                          ident.1=ident_1,ident.2=ident_2,
                          verbose=FALSE,logfc.threshold=0)
}

## Pseudobulking analysis on all the cell cell_types

pseudo_lipo=AggregateExpression(lipo,assays="RNA",return.seurat = T, group.by = c("cell.type.final", "condition", "donor"))
pseudo_lipo$celltype.condition <- paste(pseudo_lipo$cell.type.final, pseudo_lipo$condition, sep = "_")
Idents(pseudo_lipo)="celltype.condition"
DEG_lipo_bulk=list()
for (celltype in cell_types){
  ident_1 <- paste0(cell_type, "_GL")
  ident_2 <- paste0(cell_type, "_Ctrl")
  DEG_lipo[[celltype]] <- FindMarkers(pseudo_lipo,
                          ident.1=ident_1,ident.2=ident_2,
                          test.use="DESeq2")
```
Then,

* DEG Filtering |FC|>1.2 & p.adj < 0.05
* calculate DEG number of each cell type
* identify the unique DEG or common shared DEG across the cell types

```r
DEG_1.2=list()

for (j in seq_along(DEG_lipo)){
  DEG_1.2[[j]] <-
    subset(DEG_lipo[[j]], abs(avg_log2FC)>log2(1.2) &p_val_adj<0.05)
}

names(DEG_1.2)=names(DEG_lipo)

## get the number of DEGs for each cell types after filtering
DEG_1.2_count=lapply(DEG_1.2,nrow)

## find the unique or common (combination from 2 to 9) DEG across all the cell types

up_genes <- list()
down_genes <- list()

for (name in names(DEG_1.2)) {
  up_genes[[paste0(name, "_up")]] <- rownames(subset(DEG_1.2[[name]], avg_log2FC > 0))
  down_genes[[paste0(name, "_down")]] <- rownames(subset(DEG_1.2[[name]], avg_log2FC < 0))
}

## for down regulated genes
down_shared_gene_counts <- numeric(length(down_genes) - 1)

for (n in 2:length(down_genes)) {
  combs <- combn(names(down_genes), n, simplify = FALSE)
  for (comb in combs) {
    common_genes <- Reduce(intersect, down_genes[comb])
    if (n < length(down_genes)) {
      other_genes <- setdiff(names(down_genes), comb)
      for (other in other_genes) {
        common_genes <- setdiff(common_genes, down_genes[[other]])
      }
    }
    down_shared_gene_counts[n-1] <- down_shared_gene_counts[n-1] + length(common_genes)
  }
}

down_unique_N <- length(unique(unlist(down_genes)))-sum(down_shared_gene_counts)

## for up regulated genes

up_shared_gene_counts <- numeric(length(up_genes) - 1)

# Loop through combinations from 2 to 9
for (n in 2:length(up_genes)) {
  combs <- combn(names(up_genes), n, simplify = FALSE)

  # For each combination, find the common genes
  for (comb in combs) {
    common_genes <- Reduce(intersect, up_genes[comb])

    # Exclude genes that are shared by more than n data frames
    if (n < length(up_genes)) {
      other_genes <- setdiff(names(up_genes), comb)
      for (other in other_genes) {
        common_genes <- setdiff(common_genes, up_genes[[other]])
      }
    }

    # Count the unique genes shared by exactly n data frames
    up_shared_gene_counts[n-1] <- up_shared_gene_counts[n-1] + length(common_genes)
  }
}

up_unique_N <- length(unique(unlist(up_genes)))-sum(up_shared_gene_counts)

## generate the plot

df <- data.frame(cell_type_N=rep(1:9,2),
  Direction=c(rep("up",9),rep("down",9)),
  Number=c(up_unique_N,up_shared_gene_counts,down_unique_N,down_shared_gene_counts)

ggplot(df, aes(x = Celltype_N, y = Number, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge to place bars side by side
  labs(x = "Number of Cell Types", y = "Gene Count", fill = "Direction") +
  scale_fill_manual(values = c("up" = "firebrick", "down" = "steelblue"))+
  theme_classic()

```

## Augur Analysis
Augur is a R package to estimate the cell sensitivity to the env disturbation, for the details please refer to https://github.com/neurorestore/Augur

```r
library(tidyverse)
library(Seurat)
library(Augur)

lipo_subset=subset(lipo,idents=c("acinar","alpha","beta","delta","ductal",
                            "endothelial","fibroblast","immune","pp"))
DefaultAssay(lipo_subset)="RNA"

## In our case, the ranking is not sensitive to the parameter change. tree number only affect the AUC value.
lipo.augur=calculate_auc(lipoglucotoxicity,
  label_col = "condition",cell_type_col = "cell.type.final",
  rf_params=list(trees=500))

## plotting
plot_lollipop(lipo.augur)
plot_umap(lipo.augur, lipo,reduction = "umap.cca",cell_type_col = "cell.type.final")

```

## Pathway enrichment

```r
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)

enrich_GO_for_genes <- function(gene_list) {
  enrichGO(
    gene = gene_list,
    OrgDb = "org.Hs.eg.db",
    ont = "BP",
    keyType = "SYMBOL",
    qvalueCutoff = 0.05
  )
}

up_GO_results <- lapply(up_genes, function(gene_list) {
  lapply(gene_list, enrich_GO_for_genes)
})

down_GO_results <- lapply(down_genes, function(gene_list) {
  lapply(gene_list, enrich_GO_for_genes)
})

## Then each GO results were exported and further used revigo(http://revigo.irb.hr/) to pick the parental terms
```
## UPR module scores calculation

The `AddModuleScore` function is used to calculate the average expression level of UPR pathways in each cells, refer to https://satijalab.org/seurat/reference/addmodulescore

```r
library(Seurat)
library(dplyr)
library(tibble)
library(ComplexHeatmap)

lipo_subset <- subset(lipo,idents=c("alpha","beta","delta","pp","ductal","acinar"))

## resources from GO
## ATF https://amigo.geneontology.org/amigo/term/GO:0036500,
## IRE1 https://amigo.geneontology.org/amigo/term/GO:0036498
## PERK https://amigo.geneontology.org/amigo/term/GO:0036499

features=list(ATF_UPR=c("DDIT3","CREBZF","XBP1","ATF6","ATF6B","MBTPS2","MBTPS1","MANF"),
              IRE1_UPR=c("PTPN1",,"XBP1","ERN1","VAPB"),
              PERK_UPR=c("RPAP2","DDIT3","EIF2AK3","NFE2L2","QRICH1","ATF4"))

lipo_subset <- AddModuleScore(
  object = lipo_subset,
  features = features,
  name = c("ATF_UPR_Score", "IRE1_UPR_Score", "PERK_UPR_Score")
)

lipo_subset.module <- lipo_subset@meta.data[,c("cell.type.final","condition","ATF_UPR_Score1","IRE1_UPR_Score2","PERK_UPR_Score3")]

module_average <- lipo_subset.module %>%
  group_by(cell.type.final,condition) %>%
  summarise(across(c("ATF_UPR_Score1", "IRE1_UPR_Score2", "PERK_UPR_Score3"),
                   mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(cell_condition = paste(cell.type.final, condition, sep = "_")) %>%
  select(-cell.type.final,-condition) %>%
  column_to_rownames("cell_condition") %>%
  .[c("alpha_Ctrl", "alpha_GL",
      "beta_Ctrl", "beta_GL",
      "delta_Ctrl", "delta_GL",
      "pp_Ctrl", "pp_GL",
      "ductal_Ctrl", "ductal_GL",
      "acinar_Ctrl", "acinar_GL"),] %>%
  t()


##### PDF plot
colors <- colorRampPalette(c("grey","white", "orange"))(100)

pdf("UPR_score_20241105version.pdf", width = 10, height = 3)
pheatmap(as.matrix(module_average),cluster_cols=F,cluster_rows=F,
         color=colors)
dev.off()
```
## Cell-Cell communication

CellChat R package is used for this section, refer to https://github.com/jinworks/CellChat

First, generate the cellchat for downstream visualization

```r
library(CellChat)
library(patchwork)
library(Seurat)

lipo.subset=subset(lipo,idents=c("acinar","alpha","beta","delta","ductal","epsilon",
                            "endothelial","fibroblast","immune","pp"))
set.seed(123)
options(future.seed=TRUE)

## Create cellchat object, GL and Ctrl need to seperate

create_cellchat_object <- function(lipo_data, condition) {
  cell.use <- rownames(lipo_data@meta.data)[lipo_data$condition == condition]
  data_input <- lipo_data[["RNA"]]$data[, cell.use]
  metadata <- lipo_data@meta.data[cell.use, ]
  metadata_cellchat <- data.frame(labels = metadata$cell.type.final, row.names = rownames(metadata))
  createCellChat(object = data_input, meta = metadata_cellchat, group.by = "labels")
}

PA_cellchat <- create_cellchat_object(lipo.subset, "GL")
Ctrl_cellchat <- create_cellchat_object(lipo.subset, "Ctrl")

## reorder the cell order before further processing

cell.type.order=c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune")
Ctrl_cellchat=setIdent(Ctrl_cellchat,ident.use = "labels",levels =cell.type.order)
PA_cellchat=setIdent(PA_cellchat,ident.use = "labels",levels =cell.type.order)

CellchatDB=CellChatDB.human

CellchatDB.use <- CellchatDB
Ctrl_cellchat@DB <- CellchatDB.use

## Process CellChat objects

process_cellchat <- function(cellchat_obj) {
  cellchat_obj <- subsetData(cellchat_obj)
  future::plan("multisession", workers = 6)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj) %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = "triMean") %>%
    computeCommunProbPathway() %>%
    aggregateNet() %>%
    netAnalysis_computeCentrality(slot.name = "netP")
  return(cellchat_obj)
}

Ctrl_cellchat <- process_cellchat(Ctrl_cellchat)
PA_cellchat <- process_cellchat(PA_cellchat)

## Merge the results
object.list <- list(Ctrl=Ctrl_cellchat,PA=PA_cellchat)
cellchat <- mergeCellChat(object.list,add.names = names(object.list))

```
Downstream visualization follow the CellChat tutorial, except the heatmap of outgoing and incoming signaling pattern, which is modified based on the original code of `netAnalysis_signalingRole_heatmap`.

```r
## Process in-degree and out-degree for both Ctrl and GL
process_degrees <- function(object.list, degree_type, cell_types) {
  degrees <- lapply(object.list@netP$centr, function(signaling) {
    unlist(signaling[[degree_type]])
  })
  degree_df <- do.call(rbind, degrees)
  degree_df_endocrine <- degree_df[, cell_types]

  colnames(degree_df_endocrine) <- paste0(cell_types, "_", degree_type)

  return(degree_df_endocrine)
}

cell_types <- c("alpha", "beta", "delta", "epsilon", "pp")

indeg_control_df_endocrine <- process_degrees(lipo.object.list.all$Ctrl, "indeg", cell_types)
indeg_PA_df_endocrine <- process_degrees(lipo.object.list.all$PA, "indeg", cell_types)
outdeg_control_df_endocrine <- process_degrees(lipo.object.list.all$Ctrl, "outdeg", cell_types)
outdeg_PA_df_endocrine <- process_degrees(lipo.object.list.all$PA, "outdeg", cell_types)

merge_degrees <- function(control_df, PA_df) {
  combined_df <- merge(control_df, PA_df, by = "row.names", all = TRUE)
  combined_df[is.na(combined_df)] <- 0
  rownames(combined_df) <- combined_df$Row.names
  combined_df <- combined_df[, -1]

  combined_df_reorder <- combined_df[, c(
    "alpha_ctrl", "alpha_GL", "beta_ctrl", "beta_GL", "delta_ctrl", "delta_GL",
    "epsilon_ctrl", "epsilon_GL", "pp_ctrl", "pp_GL"
  )]
  combined_df_reorder_filtered <- combined_df_reorder[rowSums(combined_df_reorder) != 0, ]

  return(combined_df_reorder_filtered)
}

indeg_combined_df_reorder_filtered <- merge_degrees(indeg_control_df_endocrine, indeg_PA_df_endocrine)
outdeg_combined_df_reorder_filtered <- merge_degrees(outdeg_control_df_endocrine, outdeg_PA_df_endocrine)

generate_heatmap <- function(degree_df, color_palette) {
  degree_mat <- sweep(degree_df, 1L, apply(degree_df, 1, max), '/', check.margin = FALSE)
  color.heatmap.use <- colorRampPalette(brewer.pal(n = 9, name = color_palette))(100)
  Heatmap(degree_mat, cluster_columns = FALSE, col = color.heatmap.use, name = "Relative strength")
}

generate_heatmap(indeg_combined_df_reorder_filtered, "GnBu")
generate_heatmap(outdeg_combined_df_reorder_filtered, "BuGn")

```
## Regulon Analysis

This section used pySCENIC python package for analysis, and later AUC score is loaded into seurat object to find DE regulon. Refer to https://github.com/aertslab/pySCENIC, https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic-differential-regulons.html

Due to the inherent stochasticity of grnboost2, a consensus GRN was generated by peforming 11 independent pySCENIC runs. Only regulons that appeared in at least 5 runs were retained for the downstream analysis, and their average AUC values were used.

Firstly, generate loom from seurat for the input to pySCENIC

```r
library(Seurat)
library(SeuratDisk)
library(loomR)

lipo_subset <- subset(lipo,idents=c("alpha","beta","delta","epsilon","pp"))
lipo.endocrine.loom <- as.loom(lipo_subset,filename='lipo.endocrine.loom',verbose=FALSE)
lipo.endocrine.loom$close_all()

```
Secondly, Identify regulons and calculate the cellular regulon enrichment score (AUC).

```bash

## first, get the co-expression modules by 11 runs 
 
NUM_RUNS=11

for i in $(seq 1 $NUM_RUNS); do
    OUTPUT_FILE="lipo_endocrine_new_run_10/lipo_integrated_endocrine_new.adj.run${i}.tsv"
    echo "Running iteration $i, saving output to $OUTPUT_FILE"

    pyscenic grn \
        --num_workers 63 \
        -o $OUTPUT_FILE \
        lipo.endocrine.loom \
        allTFs_hg38.txt

## Second, filter modules for targets with cis regulatory motifs. (TF binding motifs)
NUM_RUNS=11

for i in $(seq 1 $NUM_RUNS); do
    ADJ_FILE="lipo_endocrine_new_run_10/lipo_integrated_endocrine_new.adj.run${i}.tsv"
    OUTPUT_FILE="../lipo_integrated_endocrine_regulons_new.run${i}.csv"

    echo "Running iteration $i, using input $ADJ_FILE and saving output to $OUTPUT_FILE"

    pyscenic ctx \
        "$ADJ_FILE" \
        hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname lipo.endocrine.loom \
        --mode "dask_multiprocessing" \
        --output "$OUTPUT_FILE" \
        --num_workers 51

## Calculate the cellular regulon enrichment matrix
NUM_RUNS=11
 
for i in $(seq 1 $NUM_RUNS); do
    INPUT_REGULONS="../lipo_integrated_endocrine_regulons_new.run${i}.csv"
    OUTPUT_LOOM="lipo_endocrine_new_run_10/lipo_integrated_endocrine_SCENIC.run${i}.loom"

    echo "Running iteration $i for aucell..."
    
    pyscenic aucell \
        --num_workers 51 \
        -o "$OUTPUT_LOOM" \
        lipo.endocrine.loom \
        "$INPUT_REGULONS"
done

echo "All aucell runs completed."

## Create a Scope-compatible loom file, add_visualization.py and export_to_loom.py are neeed, which can be found in the Regulon_Analysis folder
python add_visualization.py \
    --loom_input lipo_endocrine_new_run_10/lipo_integrated_endocrine_SCENIC.run1.loom \
    --loom_output lipo_visualization.loom \
    --num_workers 25

```
After getting the AUC values from 11 runs, the data were imported back to seurat to perform the downstream analysis, including

* select the high frequently occured regulons (>=5 runs)
* identify the differnetially expressed regulons in each endocrine cell types

```r
library(SCopeLoomR)
library(Seurat)
library(SeuratDisk)
library(SCENIC)
library(ComplexHeatmap)
library(ggvenn)
library(scales)

lipo_scenic <- subset(lipo,idents=c("alpha","beta","delta","epsilon","pp"))

## load the AUC scores from 11 runs into Seurat and create assay for each
loom_files=list.files(pattern = "\\.loom$")

for (loom_file in loom_files) {
   print(paste("Processing:", loom_file))
   
  # Extract run number from filename
   run_id <- gsub(".*run([0-9]+).loom", "\\1", loom_file)
   
   # Load loom file
   loom <- open_loom(loom_file)
   
   # Extract regulons
   regulons_incidMat <- get_regulons(loom, column.attr.name = 'Regulons')
   regulons <- regulonsToGeneLists(regulons_incidMat)
   
   # Extract AUC matrix
   regulonsAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
   AUCmat <- AUCell::getAUC(regulonsAUC)
   rownames(AUCmat) <- gsub("[_(+)]", "", rownames(AUCmat))
   
   # Create assay in Seurat object
   assay_name <- paste0("AUC_run", run_id)
   lipo_scenic[[assay_name]] <- CreateAssayObject(data = AUCmat)
   DefaultAssay(lipo_scenic) <- assay_name
}

## Get the AUC score for each cell type condition, and select regulons occur >=5 runs
lipo_scenic$cell.type.condition <- paste(lipo_scenic$cell.type.final, lipo_scenic$condition, sep = "_")
Idents(lipo_scenic) <- "cell.type.condition"
auc_dfs <- AverageExpression(lipo_scenic)
auc_dfs_auc <- auc_dfs[c("AUC",paste0("AUC_run", 1:11))]
regulons_list <- lapply(auc_dfs_auc, rownames)
regulon_counts <- table(unlist(regulons_list))
selected_regulon <- names(regulon_counts[regulon_counts >= 5])

## Standardize all AUC matrices to have the same regulons, filling missing values with NA
auc_dfs_auc_standardized <- lapply(auc_dfs_auc, function(df) {
  df=as.matrix(df)
  # Create a new matrix/dataframe with selected regulons
  standardized_df <- matrix(NA, nrow = length(selected_regulon), ncol = ncol(df))
  rownames(standardized_df) <- selected_regulon
  colnames(standardized_df) <- colnames(df)
  
  # Fill in values where the regulon exists in the original dataframe
  common_regulons <- intersect(selected_regulon, rownames(df))
  standardized_df[common_regulons, ] <- df[common_regulons, , drop = FALSE]
  
  return(as.data.frame(standardized_df))  # Convert to data frame if needed
})

## Average AUC score across runs
calculate_mean_auc <- function(auc_list) {
  # Convert each data frame in the list to a numeric matrix
  auc_matrices <- lapply(auc_list, function(df) as.matrix(df))
  
  # Convert the list of matrices into an array
  auc_array <- simplify2array(auc_matrices)
  
  # Compute mean while ignoring NA values
  mean_auc <- apply(auc_array, c(1,2), function(x) mean(as.numeric(x), na.rm = TRUE))
  
  return(as.data.frame(mean_auc))  # Convert to data frame
}

mean_auc_per_regulon_df <- calculate_mean_auc(auc_dfs_auc_standardized)

## AUC heatmap plot
mean_auc_per_regulon_df_subset=mean_auc_per_regulon_df[,c("alpha-Ctrl","alpha-PA","beta-Ctrl","beta-PA","delta-Ctrl","delta-PA","pp-Ctrl","pp-PA")]

pheatmap(as.matrix(mean_auc_per_regulon_df_subset),cluster_rows = TRUE,cluster_cols = FALSE,
        clustering_method = "ward.D2",clustering_distance_rows="correlation",scale = "row")

## Identify DE regulon 

auc_long_df <- lapply(names(auc_dfs_auc_standardized), function(run) {
  df <- auc_dfs_auc_standardized[[run]]
  df <- as.data.frame(df)  # Ensure it's a dataframe
  df$regulon <- rownames(df)  # Add regulon names
  df$run <- run  # Add run ID
  return(df)
}) %>%
  bind_rows() %>%  # Bind all dataframes together
  pivot_longer(cols = -c(regulon, run), names_to = "cell_type_condition", values_to = "AUC") %>%
  separate(cell_type_condition, into = c("cell_type", "condition"), sep = "-") %>%
  drop_na(AUC)  # Remove rows with NA values

# Ensure condition is a factor with correct ordering
auc_long_df <- auc_long_df %>%
  mutate(condition = factor(condition, levels = c("Ctrl", "PA"))) 

# Pivot wider to get separate columns for PA and Ctrl
auc_wide_df <- auc_long_df %>%
  pivot_wider(names_from = condition, values_from = AUC) %>%
  drop_na()  # Ensure paired values exist

# Perform the paired Wilcoxon signed-rank test for each regulon within each cell type
wilcoxon_results <- auc_wide_df %>%
  group_by(regulon, cell_type) %>%
  summarise(
    p_value = wilcox.test(PA, Ctrl, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "fdr"))  # FDR correction

# Compute log2 fold change (avg_log2FC) per regulon-cell type
log2FC_results <- auc_wide_df %>%
  group_by(regulon, cell_type) %>%
  summarise(avg_log2FC = log2(mean(PA, na.rm = TRUE) / mean(Ctrl, na.rm = TRUE)), .groups = "drop")

# Merge Wilcoxon results with log2FC
final_results <- wilcoxon_results %>%
  left_join(log2FC_results, by = c("regulon", "cell_type"))

# Define significance threshold
significance_threshold <- 0.05  # Adjusted p-value cutoff

# Filter for significantly upregulated regulons (FDR-adjusted p-value < 0.05 and log2FC > 0)
upregulated_list_0 <- final_results %>%
  filter(adj_p_value < significance_threshold, avg_log2FC > log2(1.1)) %>%
  split(.$cell_type)

# Filter for significantly downregulated regulons (FDR-adjusted p-value < 0.05 and log2FC < 0)
downregulated_list_0 <- final_results %>%
  filter(adj_p_value < significance_threshold, avg_log2FC < -log2(1.1)) %>%
  split(.$cell_type)

upregulated_list_0 <- upregulated_list_0[setdiff(names(downregulated_list), "epsilon")]
downregulated_list_0 <- downregulated_list_0[setdiff(names(downregulated_list), "epsilon")]

# Venn diagram for upregulated regulons
ggvenn(
  data=lapply(upregulated_list_0,function(df) {df$regulon}), 
  fill_color = c("#4783b5", "#ce5d5d", "#609e9f", "#88cdea"),
  stroke_size = 0.5, set_name_size = 4
  )

ggvenn(
  data=lapply(downregulated_list_0,function(df) {df$regulon}), 
  fill_color = c("#4783b5","#ce5d5d", "#609e9f", "#88cdea"),
  stroke_size = 0.5, set_name_size = 4
  )

# save AUC_mean into seurat
# Get all AUC-related assay names
auc_assays <- grep("^AUC", names(lipo_scenic@assays), value = TRUE)

# Initialize a matrix to store mean AUC values for selected regulons
auc_mean_matrix <- matrix(NA, nrow = length(selected_regulon), ncol = ncol(lipo_scenic),dimnames = list(selected_regulon, colnames(lipo_scenic)))

# Loop through each selected regulon and calculate mean AUC across assays
for (regulon in selected_regulon) {
  # Extract values across assays for the given regulon
  regulon_values <- lapply(auc_assays, function(assay_name) {
    if (regulon %in% rownames(lipo_scenic[[assay_name]]@data)) {
      return(lipo_scenic[[assay_name]]@data[regulon, , drop = FALSE])  # Keep as a matrix
    } else {
      return(matrix(NA, nrow = 1, ncol = ncol(lipo_scenic)))  # Fill with NA if missing
    }
  })
  
  # Combine the extracted values into a single matrix
  regulon_matrix <- do.call(rbind, regulon_values)

  # Compute mean across assays while ignoring NA values
  auc_mean_matrix[regulon, ] <- colMeans(regulon_matrix, na.rm = TRUE)
}

# Create a new Seurat assay for AUC_mean and add it to the object
lipo_scenic[["AUC_mean"]] <- CreateAssayObject(data = auc_mean_matrix)

saveRDS(lipo_scenic,file = "lipo_scenic_11runs.rds")

## get heatmap mean_AUC for each cell
DefaultAssay(lipo_scenic) <- "AUC_mean"
lipo_scenic <-ScaleData(lipo_scenic,assay = "AUC_mean", features = rownames(lipo_scenic@assays$AUC_mean))

Idents(lipo_scenic) <- "cell.type.condition"
Idents(lipo_scenic) <- factor(Idents(lipo_scenic), levels = c("alpha_Ctrl","alpha_PA","beta_Ctrl","beta_PA","delta_Ctrl","delta_PA","epsilon_Ctrl", "epsilon_PA","pp_Ctrl","pp_PA"))

regulons_reordered <- c("GLIS3", "NR3C1", "XBP1", "SMC3", "CREB3", "ILF2", "YY1", "ATF4", "MAFG", "NFE2L1", "CEBPB", "CEBPG", "DDIT3", "CREM", "ETS2", "NFIC","PDX1", "RXRG", "NKX6-1", "MAFA", "MNX1", "ASCL2", "SREBF1", "GTF3A", "FOXA2", "HLTF", "ESRRA", "CHURC1", "JUND", "THAP11", "ELK1", "TAF6", "NFYB", "DRAP1", "HTATIP2", "ZNF639", "MLX", "RFXANK", "MAZ", "MXD4", "MZF1", "ATF1", "BHLHE41","HSF1", "E2F4", "PURA", "USF2", "ZNF580", "CREB3L2", "ISL1", "PAX6", "MAFB", "NEUROD1", "JUNB", "FOS", "ARX", "FEV", "GTF2I", "GATA6", "PGR", "STAT3","BCL6", "STAT1", "CEBPD", "IRF1", "ELF3", "KLF6", "HNF4A", "IRF6", "ARNT2","ETV6", "FOXO1", "PBX1", "ELF2", "BACH2", "KLF7", "RFX3", "TEAD1", "ZNF148","SOX6", "TCF12", "NFIL3", "NFE2L2", "JUN", "MAFF", "NFKB2", "RELB", "BPTF","CHD2", "REL", "CREB5", "FOSL2", "ZEB1", "ELF1", "FOXO3", "NFKB1", "MXD1","RELA", "ATF3", "ZBTB21", "HMGA1", "MAX", "ZBTB11", "ZNF281", "SOX4", "EBF1","ETV1", "BCLAF1", "SMAD3", "BACH1", "NR1H4", "ATF6", "ZBTB43", "GABPA", "HLF","NR2F2", "ARID3A", "KMT2B", "MXI1", "FOXA3", "IRF2", "MEIS1", "EP300", "EGR1","SMAD1", "RXRA", "SMARCC2", "FOXP4", "IRF3", "ZNF76", "GATAD1", "NKX2-3")

lipo_scenic_1=subset(lipo_scenic,idents = c("alpha_Ctrl","alpha_PA","beta_Ctrl","beta_PA","delta_Ctrl","delta_PA","pp_Ctrl","pp_PA"))

DoHeatmap(lipo_scenic_1,features = regulons_reordered,raster = F)

## beta ranking based on p value vs magma ranking from
beta_DE_regulon <- final_results %>%
  filter(cell_type=="beta") %>%
  arrange(p_value) %>%
  mutate(ranking_beta = c(1:length(beta_DE_regulon$regulon))) %>%
  select(regulon,ranking_beta)

magma=read.csv("/Users/liguo/Desktop/lipoglucotoxicity/all_samples/SCENIC/DATA/magma.csv") %>%
  select(SYMBOL,MAGMA_RNK)

colnames(beta_DE_regulon)=c("SYMBOL","Beta_RNK")
common_data <- inner_join(beta_DE_regulon, magma, by = "SYMBOL")

common_data <- common_data %>%
  mutate(MAGMA_RNK_RANK = rank(MAGMA_RNK, ties.method = "first"))

highlight_data <- subset(common_data, MAGMA_RNK_RANK <= 50 & Beta_RNK <=50)

ggplot(common_data, aes(x = MAGMA_RNK_RANK, y = Beta_RNK, label = SYMBOL)) +
  geom_point(color = "black",size=1)+
  geom_point(data = highlight_data, color = "firebrick",size=1) +
  geom_text(data = highlight_data,aes(x = MAGMA_RNK_RANK, y = Beta_RNK, label = SYMBOL),
    vjust = -1, hjust = 1,color="firebrick") +
  labs(x = "MAGMA_RNK",
       y = "Beta_p_value") +
  scale_x_continuous(trans = log2_trans(),
                     breaks = c(1,50,100,120))+
  scale_y_continuous(trans = log2_trans(),
                     breaks = c(1,50,100,120))+
  geom_vline(xintercept = 50, linetype = "dashed") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  coord_fixed(ratio = 1)+
  theme_classic()

```
THe conventional post-pySCENIC analysis was also done in python following the tutorial https://pyscenic.readthedocs.io/en/latest/tutorial.html,
the code is included in the *Regulon_Analysis* folder in file *post_pyscenic.ipynb*.

## Pseudotime analysis

`slingshot` was used for pseudotime analysis, modified from tutorial https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html 
and tutorial https://hectorrdb.github.io/condimentsPaper/articles/TGFB.html

``` r
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(scales)
library(viridis)
library(UpSetR)
library(pheatmap)
library(msigdbr)
library(fgsea)
library(knitr)
library(ggplot2)
library(gridExtra)
library(condiments)
library(dplyr)
library(tradeSeq)
library(cowplot)

# Selecting beta cells
Idents(lipo)=lipo$cell.type.final
beta = subset(lipo,idents="beta")

# Converting to singleCellExperiment
beta_sce <- as.SingleCellExperiment(beta, assay = "RNA")
beta_sce = slingshot(beta_sce,reducedDim='UMAP.CCA', clusterLabels=colData(beta_sce)$condition,start.clus="Ctrl", approx_points=150)

# Differnetial topology test
set.seed(821)
library(tradeSeq)
BPPARAM <- BiocParallel::MulticoreParam(workers = 20)

## fit negative binomial GAM
## filter genes to have minimum shared counts of 20, and minimum expressed in 3 cells, in line with scVelo default, and to improve fitGAM execusion speed
gene_cells <- rowSums(counts(beta_sce)!=0)
beta_sce_filtered <- beta_sce[gene_cells > 3,]  
gene_sums <- rowSums(counts(beta_sce_filtered))
beta_sce_filtered <-  beta_sce_filtered[gene_sums > 20,]

icMat <- evaluateK(counts=beta_sce_filtered,
                   nGenes = 500,
                   k = 3:10, 
                   conditions=factor(beta_sce$condition),
                   parallel=TRUE)

plot_evalutateK_results(icMat,k=3:10)

beta_sce_filtered <- fitGAM(beta_sce_filtered, nknots=7, parallel=TRUE) 


## Differential expression along pseudotime
ATres <- associationTest(beta_sce_filtered) # testing whether the average gene expression is significantly changed along pseudotime
ATres$p.adjust <- p.adjust(ATres$pvalue,"fdr")

## Plot marker genes in the UPR term
genelist= c("TMED2","ATF6","HERPUD1","PTPN1","NFE2L2","XBP1","DNAJC10","ATF4","DDIT3","HSPA5",
            "ERN1","CREBZF","EIF2S1","RPAP2","QRICH1","EIF2AK3","MBTPS2","PARP16","ATF6B","MBTPS1","VAPB")
scales <- brewer.pal(7,"Accent")[1:2]
for (i in genelist)
{
  gene = i
  print(plotSmoothers(beta_sce_filtered, assays(beta_sce_filtered)$counts,
                      gene = gene,
                      alpha = 1, border = TRUE, curvesCols=scales)+
          scale_color_manual(values=scales)+
          ggtitle(gene))
}

## Heatmaps of genes whose expression vary over pseudotime
pseudotime_genes <- rownames(ATres)[
  which((ATres$p.adjust < 0.01)  )
]

## based on mean smoother
yhatSmooth <- 
  predictSmooth(beta_sce_filtered, gene = pseudotime_genes, nPoints = 100, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth,1, scales::rescale))

heatSmooth <- pheatmap(yhatSmoothScaled,
                       cluster_cols = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE, 
                       clustering_distance_rows="correlation",
                       clustering_method="ward.D2",
                       legend = TRUE
)


my_gene_col <- cutree(heatSmooth$tree_row, k=3)
my_gene_col <- data.frame(my_gene_col)
my_colour = list(
  my_gene_col = brewer.pal(7,"Accent")[1:3]
)


pheatmap(yhatSmoothScaled,
         cluster_cols = FALSE,
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         clustering_distance_rows="correlation",
         clustering_method="ward.D2",
         annotation_row=my_gene_col,
         annotation_colors=my_colour,
         annotation_names_row = FALSE,
         legend = TRUE
)
```


