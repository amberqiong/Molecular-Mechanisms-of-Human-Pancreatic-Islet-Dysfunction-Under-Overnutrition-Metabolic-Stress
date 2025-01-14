# Molecular-Mechanisms-of-Human-Pancreatic-Islet-Dysfunction-Under-Overnutrition-Metabolic-Stress

This repository hosts the codes to generate analysis and figures for Molecular Mechanisms of Human Pancreatic Islet Dysfunction Under Overnutrition Metabolic stress, currently under review in Diabetes

## Pre-processing

Before making the actual analysis, begin with ambient mRNA decontamination

## Filtering

## Integration and cell calling

## DEG Analysis
This section outlines the steps to perform DEG analysis for all cell types comparing the conditions `GL` vs `Ctrl`. The `FindMarkers` function from the Seurat package is used for differential expression analysis.

For a detailed explanation of `FindMarkers`, refer to the [Seurat Differential Expression Vignette](https://satijalab.org/seurat/articles/de_vignette).

```
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

###DEG Filtering |FC|>1.2 & p.adj < 0.05
```
DEG_1.2=list()

for (j in seq_along(DEG_lipo)){
  DEG_1.2[[j]] <-
    subset(DEG_lipo[[j]], abs(avg_log2FC)>log2(1.2) &p_val_adj<0.05)
}

names(DEG_1.2)=names(DEG_lipo)

## get the number of DEGs for each cell types after filtering
DEG_1.2_count=lapply(DEG_1.2,nrow)
```

## Augur Analysis
Augur is a R package to estimate the cell sensitivity to the env disturbation, for the details please refer to https://github.com/neurorestore/Augur

```
lipo_subset=subset(lipo,idents=c("acinar","alpha","beta","delta","ductal",
                            "endothelial","fibroblast","immune","pp"))
DefaultAssay(lipo_subset)="RNA"

## In our case, the ranking is not sensitive to the parameter change. tree number only affect the AUC value.
lipo.augur=calculate_auc(lipoglucotoxicity,
  label_col = "condition",cell_type_col = "cell.type.final",
  rf_params=list(trees=500))

```



Summary:

* GSIS (Fig.1)

* Preprocessing and cell calling (Fig.1)

* DEG and GO enrichment, cell sensitivity to env (Fig.2)

* Cell-Cell communication(Fig.3)

* Regulon Analysis (Fig.4)

* Pseudotime Analysis (Fig.5)
