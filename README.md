# Molecular-Mechanisms-of-Human-Pancreatic-Islet-Dysfunction-Under-Overnutrition-Metabolic-Stress

This repository hosts the codes to generate analysis and figures for Molecular Mechanisms of Human Pancreatic Islet Dysfunction Under Overnutrition Metabolic stress, currently under review in Diabetes

## Pre-processing

Before making the actual analysis, begin with ambient mRNA decontamination

## Filtering

## Integration and cell calling

```r
library(dplyr)
library(Seurat)
library(ggplot2)

## features dot plot
Idents(lipo) <- factor(x=Idents(lipoglucotoxicity),levels=rev(c("alpha","beta","delta","epsilon","pp","ductal","acinar","fibroblast","endothelial","immune","doublets","unknown")))

features=rev(c("GCG","INS","SST","GHRL","PPY","CFTR","CPA2","SPARC","VWF","PTPRC"))
DotPlot(lipoglucotoxicity,features = features)

## color code for different cell types
colors <- c("#4682B4","#CD5C5C", "#5F9EA0", "firebrick","#87CEEB", "#FF8C00", "#48D1CC","#FFD700", "#7B68EE","#FF6347","darkgrey")

DimPlot(lipo, reduction = "umap.cca",cols = colors,label = TRUE)+labs(x="UMAP1",y="UMAP2")

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

## resources from GO, updated in 2024-11-05
## ATF https://amigo.geneontology.org/amigo/term/GO:0036500,
## IRE1 https://amigo.geneontology.org/amigo/term/GO:0036498
## PERK https://amigo.geneontology.org/amigo/term/GO:0036499

features=list(ATF_UPR=c("DDIT3","CREBZF","XBP1","ATF6","ATF6B","MBTPS2","MBTPS1"),
              IRE1_UPR=c("PARP16","PTPN1","DNAJC10","XBP1","ERN1","ERN2","VAPB"),
              PERK_UPR=c("TMED2","RPAP2","DDIT3","EIF2AK3","EIF2S1","NFE2L2","QRICH1","ATF4"))

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

Firstly, generate loom from seurat for the input to pySCENIC

```r
library(Seurat)
library(SeuratDisk)
library(loomR)

lipo_subset <- subset(lipo,idents=c("alpha","beta","delta","pp"))
lipo.endocrine.loom <- as.loom(lipo_subset,filename='lipo.endocrine.loom',verbose=FALSE)
lipo.endocrine.loom$close_all()

```
Secondly, Define and regulons and calculate the cellular regulon enrichment score (AUC).

```bash

## first, get the co-expression modules
pyscenic grn \
    --num_workers 63 \
    -o /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo.endocrine.adj.tsv \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo.endocrine.loom\
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/allTFs_hg38.txt

## Second, filter modules for targets with cis regulatory motifs. (TF binding motifs)
pyscenic ctx \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo.endocrine.adj.tsv \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo.endocrine.loom \
    --mode "dask_multiprocessing" \
    --output /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine_regulons.csv \
    --num_workers 51

## Calculate the cellular regulon enrichment matrix
pyscenic aucell \
    --num_workers 51 \
    -o /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine_SCENIC.loom \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo.endocrine.loom \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine_regulons.csv

## Create a Scope-compatible loom file, add_visualization.py and export_to_loom.py are neeed, which can be found in the Regulon_Analysis folder
python add_visualization.py \
    --loom_input lipo_integrated_endocrine_SCENIC.loom \
    --loom_output lipo_visualization.loom \
    --num_workers 25

```
After getting the AUC for each cells, let us import the data back to Seurat and using `FindMarkers` to find the differentially expressed regulons.

```r
library(SCENIC)
library(SCopeLoomR)

lipo_subset <- subset(lipo,idents=c("alpha","beta","delta","pp"))

## load the AUC score
loom <- open_loom("lipo_visualization.loom") # the scope-compatile loom is not necessary for the downstream analysis, output from aucell is enough.
regulons_incidMat <- get_regulons(loom,column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("[_(+)]", "", rownames(AUCmat))
lipo_subset[['AUC']] <- CreateAssayObject(data = AUCmat)

## identify the DE regulons
DefaultAssay(lipo_subset) <- "AUC"
lipo_subset$cell.type.condition=paste(lipo_subset$cell.type.final,lipo_subset$condition,sep="_")

cell_types=c("alpha","beta","delta","pp")
DE_regulon=list()
for (celltype in cell_types){
  ident_1 <- paste0(cell_type, "_GL")
  ident_2 <- paste0(cell_type, "_Ctrl")
  DE_regulon[[celltype]] <- FindMarkers(lipo_subset,
                          ident.1=ident_1,ident.2=ident_2,
                          verbose=FALSE,logfc.threshold=0)
}

```
THe conventional post-pySCENIC analysis was also done in python following the tutorial https://pyscenic.readthedocs.io/en/latest/tutorial.html,
the code is included in the *Regulon_Analysis* folder in file *post_pyscenic.ipynb*.

