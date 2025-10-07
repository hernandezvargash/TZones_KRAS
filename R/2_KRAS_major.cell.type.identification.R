# KRAS model - Cell Type Identification

# libraries ---------------------------------------------------------------

rm(list = ls())

setwd()

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(SingleCellExperiment)
  library(dplyr)
  library(xlsx)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(gridExtra)
  library(RColorBrewer)
  library(patchwork)
  library(hrbrthemes)
  library(pals)
  library(SingleR)
  library(scRNAseq)
  library(scibet)
  library(slingshot, quietly = T)
  library(mclust, quietly = T)
  library(celldex)
  library(KernSmooth)
  library(dendextend)
  library(circlize)
  library(enrichR)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(pathview)
  library(ReactomePA)
  library(ggpubr)
  library(eulerr)
  library(alluvial)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(biomaRt)
  library(Scillus)
  library(snow)
})

set.seed(1)


# initial markers -------------------------------------------------

# rm(list=ls())

int.sct <- readRDS("KRAS_model/objects/integrated.rds")

DefaultAssay(int.sct) <- "SCT"
# DimPlot(int.sct, label = T)
int.sct <- PrepSCTFindMarkers(int.sct)

all.markers <- FindAllMarkers(int.sct, only.pos = T, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]
head(all.markers)
table(all.markers$cluster)
tail(all.markers)

all.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "KRAS_model/results/cell_types/positive.markers.csv")
# write.xlsx(all.markers, file="KRAS_model/results/cell_types/positive.markers.xlsx")



# cell cycle adjustment ---------------------------------------------------

DefaultAssay(int.sct) <- "integrated"
DimPlot(int.sct, label = T)

RunPCA(int.sct, features = VariableFeatures(int.sct), ndims.print = 1:10, nfeatures.print = 10)
# PC 10 contains cell cycle genes
VizDimLoadings(int.sct, dims = 1:12, reduction = "pca", ncol = 6)
# DimHeatmap(int.sct, dims = c(1:12))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

int.sct <- CellCycleScoring(int.sct, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(int.sct[[]])
RidgePlot(int.sct, features = c("Mki67", "Top2a", "Birc5", "Mcm6"), ncol = 2)

# regress out cell cycle scores

# cc.int.sct <- ScaleData(int.sct, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(int.sct))
# RunPCA(cc.int.sct, features = VariableFeatures(cc.int.sct), ndims.print = 1:10, nfeatures.print = 10)

# Alternate Workflow (from Seurat's website):
# The procedure above removes all signal associated with cell cycle.
# In some cases, we’ve found that this can negatively impact downstream analysis,
# particularly in differentiating processes (like murine hematopoiesis),
# where stem cells are quiescent and differentiated cells are proliferating (or vice versa).
# In this case, regressing out all cell cycle effects can blur the distinction between stem and progenitor cells as well.
# As an alternative, we suggest regressing out the difference between the G2M and S phase scores.
# This means that signals separating non-cycling cells and cycling cells will be maintained,
# but differences in cell cycle phase among proliferating cells (which are often uninteresting), will be regressed out of the data

# check also: https://github.com/satijalab/seurat/issues/1679

int.sct$CC.Difference <- int.sct$S.Score - int.sct$G2M.Score

cc.int.sct <- SCTransform(int.sct, vars.to.regress = "CC.Difference")
cc.int.sct <- RunPCA(cc.int.sct, features = VariableFeatures(cc.int.sct), ndims.print = 1:10, nfeatures.print = 10) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5)

VizDimLoadings(cc.int.sct, dims = 1:12, reduction = "pca", ncol = 6)
# DimHeatmap(cc.int.sct, dims = c(1:12))

DimPlot(cc.int.sct, label = T)

p1 <- DimPlot(cc.int.sct)
p2 <- DimPlot(cc.int.sct, group.by = "sample")
p1 + p2

DimPlot(cc.int.sct, split.by = "sample")


saveRDS(cc.int.sct, file = "KRAS_model/objects/integrated_cellcycle.rds")


## new markers ----

table(cc.int.sct$sample)

DefaultAssay(cc.int.sct) <- "SCT"

all.markers <- FindAllMarkers(cc.int.sct, only.pos = T, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]
head(all.markers)
table(all.markers$cluster)


write.csv(all.markers, file = "KRAS_model/results/positive.markers.cellcycle.adjusted.csv")



## re-label ----

p1 <- DimPlot(cc.int.sct, label = T, label.size = 4, label.box = T) + NoLegend()

# after checking with Chloe/Geraldine:
cc.int.sct <- RenameIdents(cc.int.sct,
  "0" = "Myeloid",
  "1" = "Stratified",
  "2" = "Rectum",
  "3" = "T cells",
  "4" = "Rectum",
  "5" = "Stratified",
  "6" = "Stratified", # Anal TZ
  "7" = "T cells",
  "8" = "Anal gland",
  "9" = "Stratified",
  "10" = "Hair follicle",
  "11" = "Rectum",
  "12" = "Myeloid",
  "13" = "Rectum",
  "14" = "Myeloid",
  "15" = "Rectum",
  "16" = "Stratified",
  "17" = "B cells",
  "18" = "Rectum",
  "19" = "Stratified",
  "20" = "Stromal",
  "21" = "Rectum"
)

p2 <- DimPlot(cc.int.sct, label = T, label.box = T) + NoLegend()

p1 + p2


cc.int.sct$major.labels <- Idents(cc.int.sct)
head(cc.int.sct[[]])



saveRDS(int.sct, file = "KRAS_model/objects/integrated_annotated.rds")






# end ---------------------------------------------------------------------
