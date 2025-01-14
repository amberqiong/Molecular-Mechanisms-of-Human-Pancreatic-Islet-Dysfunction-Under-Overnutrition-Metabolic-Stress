##### script for loom in cpgm

library(Seurat)


ER.stress=readRDS("ERstress.integrated.042424.rds")
Idents(ER.stress)="cell.type.final"

ER.stress.endocrine=subset(ER.stress,idents=c("alpha","beta","delta","PP","epsilon"))

saveRDS(ER.stress.endocrine,file="ER.endocrine.0829.rds")

###loom this loom does not contain all the attribute for the pyscenic analysis, need to be modified to avoid the further errors.
library(loomR)
library(SeuratDisk)

ER.stress.endocrine=readRDS("/Users/liguo/Desktop/ERstress/ER.endocrine.0829.rds")

ER.stress.endocrine.loom <- as.loom(ER.stress.endocrine, filename = "ER.stress.endocrine.loom", verbose = TRUE)
ER.stress.endocrine.loom$close_all()


