# KRAS model - Preprocessing after ZsGreen alignment

# summary -----------------------------------------------------------------

# Integration of 5 different experiments:

# Normal (RUN 380 – Experiment 07/07/2021)
# Hyperplasia (RUN 381 – Experiment 20/07/2021)
# Hyperplasia (RUN 447 – Experiment 21/02/2023)
# Dysplasia (RUN 400 - Experiment 14/12/2021)
# Carcinoma (RUN 354 – Experiment 02/03/2021)

# Single Cell Integration + Velocity
# Fast integration using reciprocal PCA (RPCA)
# https://satijalab.org/seurat/articles/integration_rpca.html


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
  library(stringr) # library(xlsx);
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
  library(bumphunter)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(biomaRt)
  library(Scillus)
  library(Nebulosa)
  library(SingleCellSignalR)
  library(snow)
  library(DoubletFinder)
})

set.seed(1)


# loading from barcode matrices -----------------------------------------------------------------

dirs <- list.dirs(path = "/KRAS/matrices/ZsGreen", recursive = F, full.names = T)
names <- str_sub(dirs, 66, -1)

list.samples <- list()

for (i in 1:length(dirs)) {
  list.samples[[i]] <- Read10X(data.dir = dirs[i])
  list.samples[[i]] <- CreateSeuratObject(list.samples[[i]])
  list.samples[[i]]$sample <- names[[i]]
}

names(list.samples) <- names
lapply(list.samples, dim)

list.samples <- lapply(X = list.samples, FUN = PercentageFeatureSet, pattern = "^mt-", col.name = "percent.mt")
list.samples <- lapply(X = list.samples, FUN = PercentageFeatureSet, pattern = "^Rp[sl][[:digit:]]", col.name = "percent.ribo")
# list.samples <- lapply(X = list.samples, FUN = subset, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)

lapply(list.samples, dim)
head(list.samples[[5]])

object.list <- list.samples

object.list <- lapply(X = object.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = features)
object.list <- lapply(X = object.list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(
  object.list = object.list, normalization.method = "SCT",
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5
) # initially used k.anchor = 20, which may be way too high
int.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
int.sct$sample <- factor(int.sct$sample, levels = c("normal", "hyperplasia.1", "hyperplasia.2", "dysplasia", "carcinoma"))
table(int.sct$sample)
# normal hyperplasia.1 hyperplasia.2     dysplasia     carcinoma
# 2333           616          1190          5682          6352

int.sct <- RunPCA(int.sct, verbose = FALSE)
int.sct <- FindNeighbors(int.sct, dims = 1:10)
int.sct <- RunUMAP(int.sct, reduction = "pca", dims = 1:30)
int.sct <- FindClusters(int.sct, resolution = 0.5)


VlnPlot(int.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)
VlnPlot(int.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), group.by = "sample", ncol = 4)

DimPlot(int.sct, label = T)

p1 <- DimPlot(int.sct)
p2 <- DimPlot(int.sct, group.by = "sample")
p1 + p2

DimPlot(int.sct, split.by = "sample")

head(int.sct[[]])

# save(int.sct, file = "ZsGreen/ZsGreen.sc.sct.integrated.filtered.RData")
# save(int.sct, file = "ZsGreen/ZsGreen.sc.sct.integrated.RData")



# inspection --------------------------------------------------------------


DefaultAssay(int.sct) <- "SCT"

FeaturePlot(int.sct, c("Krt17", "ZsGreen"), order = T)
plot_density(int.sct, c("Krt17", "ZsGreen"))
FeatureScatter(int.sct, "Krt17", "ZsGreen")
VlnPlot(int.sct, c("Krt17", "ZsGreen"))

DefaultAssay(int.sct) <- "RNA"

length(WhichCells(int.sct, expression = Krt17 > 0)) # 6021 / 8614 (unfiltered)
length(WhichCells(int.sct, expression = ZsGreen > 0)) # 3519 / 4120 (unfiltered)

sel.cells <- WhichCells(int.sct, expression = ZsGreen > 0)
p1 <- plot_density(int.sct, "Krt17")
p2 <- plot_density(int.sct, "ZsGreen")
p3 <- DimPlot(int.sct,
  label = F, label.box = T, repel = T, # cols = "darkslateblue",
  cells.highlight = sel.cells, cols.highlight = "springgreen4", sizes.highlight = 0.2
) +
  ggtitle("ZsGreen+ cells") +
  NoLegend()

p1 + p2 + p3

# save(sel.cells, file = "ZsGreen/ZsGreen.cells.RData")


sel.cells.2 <- WhichCells(int.sct, expression = ZsGreen > 0) # unfiltered
test <- subset(int.sct, cells = sel.cells.2)

test.sub <- subset(test, sample == "normal")
sel.names.normal <- paste0("normal_130721:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.normal)
test.sub <- subset(test, sample == "hyperplasia.1")
sel.names.hyperplasia.1 <- paste0("KRAS_260721:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.hyperplasia.1)
test.sub <- subset(test, sample == "hyperplasia.2")
sel.names.hyperplasia.2 <- paste0("hyperplasia_220223:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.hyperplasia.2)
test.sub <- subset(test, sample == "dysplasia")
sel.names.dysplasia <- paste0("dysplasia_281221:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.dysplasia)
test.sub <- subset(test, sample == "carcinoma")
sel.names.carcinoma <- paste0("wound_2:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.carcinoma)

sel.cells.names <- c(sel.names.normal, sel.names.hyperplasia.1, sel.names.hyperplasia.2, sel.names.dysplasia, sel.names.carcinoma)

# save(sel.cells.names, file = "ZsGreen/ZsGreen.cells.prefilter.RData")



# add labels from loom integrated object ----------------------------------

load("ZsGreen/ZsGreen.sc.sct.integrated.filtered.RData")
head(int.sct)

# replace cell names in metadata
meta <- read.csv("metadata.sc.all.celltypes.csv", row.names = 1)
table(meta$sample)
table(as.numeric(as.factor(meta$sample)))
# meta$cellnames <- paste0(str_sub(rownames(meta), -17, -2), "-1_", as.numeric(as.factor(meta$sample))) # same result
meta$cellnames <- paste0(str_sub(str_split(rownames(meta), ":", simplify = T)[, 2], 1, -2), "-1_", as.numeric(as.factor(meta$sample)))
head(meta)
tail(meta)
intersect(colnames(int.sct), meta$cellnames)
setdiff(colnames(int.sct), meta$cellnames)
setdiff(meta$cellnames, colnames(int.sct))


DefaultAssay(int.sct) <- "SCT"
Idents(int.sct) <- int.sct$seurat_clusters
DimPlot(int.sct)
int.sct$major.labels <- "NA"
Idents(int.sct) <- int.sct$major.labels

meta$major.labels <- factor(meta$major.labels)
celltypes <- levels(meta$major.labels)
for (l in 1:length(celltypes)) {
  meta.temp <- meta[meta$major.labels == celltypes[l], ]
  cells.temp <- meta.temp$cellnames
  cells.temp <- intersect(colnames(int.sct), cells.temp)
  int.sct <- SetIdent(int.sct, cells = cells.temp, value = celltypes[l])
}

table(Idents(int.sct))
int.sct <- subset(int.sct, idents = "NA", invert = T)

p1 <- DimPlot(int.sct, cols = glasbey())
p2 <- plot_density(int.sct, "Krt17")
p3 <- plot_density(int.sct, "ZsGreen")
sel.cells <- WhichCells(int.sct, expression = ZsGreen > 0)
p4 <- DimPlot(int.sct,
  label = F, label.box = T, repel = T, # cols = "darkslateblue",
  cells.highlight = sel.cells, cols.highlight = "springgreen4", sizes.highlight = 0.2
) +
  ggtitle("ZsGreen+ cells") +
  NoLegend()

p1 + p2 + p3 + p4

DefaultAssay(int.sct) <- "RNA"
VlnPlot(int.sct, "ZsGreen", cols = glasbey()) +
  scale_y_continuous(breaks = 1:20)


# based on this, select ZsGreen > 1 as GFP+

# load("ZsGreen/ZsGreen.sc.sct.integrated.filtered.RData")

load(file = "ZsGreen/ZsGreen.sc.sct.integrated.RData")

DefaultAssay(int.sct) <- "RNA"

sel.cells.2 <- WhichCells(int.sct, expression = ZsGreen > 1) # unfiltered
test <- subset(int.sct, cells = sel.cells.2)

test.sub <- subset(test, sample == "normal")
sel.names.normal <- paste0("normal_130721:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.normal)
test.sub <- subset(test, sample == "hyperplasia.1")
sel.names.hyperplasia.1 <- paste0("KRAS_260721:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.hyperplasia.1)
test.sub <- subset(test, sample == "hyperplasia.2")
sel.names.hyperplasia.2 <- paste0("hyperplasia_220223:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.hyperplasia.2)
test.sub <- subset(test, sample == "dysplasia")
sel.names.dysplasia <- paste0("dysplasia_281221:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.dysplasia)
test.sub <- subset(test, sample == "carcinoma")
sel.names.carcinoma <- paste0("wound_2:", str_sub(colnames(test.sub), 1, -5), "x")
head(sel.names.carcinoma)

sel.cells.names <- c(sel.names.normal, sel.names.hyperplasia.1, sel.names.hyperplasia.2, sel.names.dysplasia, sel.names.carcinoma)

save(sel.cells.names, file = "ZsGreen/ZsGreen.cells.prefilter.1.RData")


# add GFP info into loom integrated object --------------------------------


load("sc.all.celltypes.RData")

table(cc.int.sct$sample)
head(cc.int.sct)
int1 <- intersect(sel.cells.names, colnames(cc.int.sct)) # done with ZsGreen > 0 ad ZsGreen > 1

# DimPlot(cc.int.sct, label = T, label.box = T) + NoLegend()
p1 <- DimPlot(cc.int.sct,
  label = T, label.box = T, repel = T, # cols = "darkslateblue",
  cells.highlight = sel.cells.names, cols.highlight = "springgreen4", sizes.highlight = 0.5
) + NoLegend()
# FeaturePlot(cc.int.sct, "Krt17")
p2 <- plot_density(cc.int.sct, "Krt17")

p1 + p2

# cc.int.sct$GFP <- "neg" for ZsGreen > 0
cc.int.sct$GFP.1 <- "neg" # for ZsGreen > 1
Idents(cc.int.sct) <- cc.int.sct$GFP.1
cc.int.sct <- SetIdent(cc.int.sct, cells = int1, value = "pos")
DimPlot(cc.int.sct)
cc.int.sct$GFP.1 <- Idents(cc.int.sct)
table(cc.int.sct$major.labels, cc.int.sct$GFP.1)

Idents(cc.int.sct) <- cc.int.sct$major.labels

save(cc.int.sct, file = "sc.all.celltypes.RData")

meta <- cc.int.sct@meta.data
head(meta)
write.csv(meta, file = "metadata.sc.all.celltypes.csv")



# end ---------------------------------------------------------------------

sessionInfo()
