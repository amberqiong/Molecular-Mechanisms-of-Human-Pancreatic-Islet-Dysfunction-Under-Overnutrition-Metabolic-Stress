###Unfolded Protein Response add module score

# load the data and subset endocrine and exocrine cells
library(Seurat)
library(dplyr)
library(tibble)
library(ComplexHeatmap)


lipo <- readRDS("lipo.rds")
Idents(lipo) <- "cell.type.final"
lipo_subset <- subset(lipo,idents=c("alpha","beta","delta","pp","ductal","acinar"))

#resources from GO, updated in 2024-11-05
#ATF https://amigo.geneontology.org/amigo/term/GO:0036500, 
#IRE1 https://amigo.geneontology.org/amigo/term/GO:0036498
#PERK https://amigo.geneontology.org/amigo/term/GO:0036499

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
  .[c("alpha_Ctrl", "alpha_PA", 
      "beta_Ctrl", "beta_PA", 
      "delta_Ctrl", "delta_PA", 
      "pp_Ctrl", "pp_PA", 
      "ductal_Ctrl", "ductal_PA", 
      "acinar_Ctrl", "acinar_PA"),] %>%
  t()
  

##### PDF plot
colors <- colorRampPalette(c("grey","white", "orange"))(100)

pdf("UPR_score_20241105version.pdf", width = 10, height = 3)

pheatmap(as.matrix(module_average),cluster_cols=F,cluster_rows=F,
         color=colors)
dev.off()