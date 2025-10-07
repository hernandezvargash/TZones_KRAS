# KRAS model - subclustering within immune cells


# libraries ---------------------------------------------------------------

rm(list = ls())

setwd()

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(SeuratDisk) # library(SeuratWrappers); library(SingleCellExperiment)
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
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(biomaRt)
  library(Scillus)
  library(Nebulosa)
  library(snow)
})


set.seed(1)


# subsetting after re-integration ---------------------------------------------

cc.int.sct <- readRDS(file = "KRAS_model/objects/integrated_annotated.rds")
head(cc.int.sct[[]])
DimPlot(cc.int.sct)
plot_density(cc.int.sct, "Ptprc")

sc.sub <- subset(cc.int.sct, idents = c("Myeloid", "T cells", "B cells"))
table(sc.sub$sample)
# Normal Hyperplasia.1 Hyperplasia.2     Dysplasia     Carcinoma
# 1319           216           288          1751          2162
rm(cc.int.sct)

sc.sub <- SplitObject(sc.sub, split.by = "sample")
sc.sub <- lapply(sc.sub, SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(sc.sub, nfeatures = 3000)
sc.sub <- PrepSCTIntegration(sc.sub, anchor.features = features)
sc.sub <- lapply(X = sc.sub, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(
  object.list = sc.sub, normalization.method = "SCT",
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5
) # initially used k.anchor = 20, which may be way too high
sc.sub <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 50, dims = 1:30) # had to reduce k.weight probably because of low number of cells
sc.sub <- RunPCA(sc.sub)
sc.sub <- RunUMAP(sc.sub, reduction = "pca", dims = 1:10)
sc.sub <- FindNeighbors(sc.sub, dims = 1:30)
sc.sub <- FindClusters(sc.sub, resolution = 0.5)

sc.sub$sample <- gsub("Hyperplasia.1", "Hyperplasia", sc.sub$sample)
sc.sub$sample <- gsub("Hyperplasia.2", "Hyperplasia", sc.sub$sample)

sc.sub$sample <- factor(sc.sub$sample,
  levels = c(
    "Normal",
    "Hyperplasia",
    "Dysplasia",
    "Carcinoma"
  )
)

p1 <- DimPlot(sc.sub, label = T, label.box = T, label.size = 5)
p2 <- DimPlot(sc.sub, group.by = "sample")

p1 + p2

DimPlot(sc.sub, split.by = "sample") # , pt.size = 1)

table(sc.sub$sample)
# Normal Hyperplasia   Dysplasia   Carcinoma
# 1319         504        1751        2162
table(Idents(sc.sub))

saveRDS(sc.sub, file = "KRAS_model/objects/immune.rds")




# markers -----------------------------------------------------------------

options(Seurat.object.assay.version = "v3") # for backward compatibility with DoubletFinder

sc.sub <- readRDS("KRAS_model/objects/immune.rds")

all.markers <- FindAllMarkers(sc.sub, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]
head(all.markers)
table(all.markers$cluster)
tail(all.markers)

all.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "KRAS_model/results/immune.markers.reclustering.csv")

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

jpeg("KRAS_model/results/immune.heatmap.jpeg", quality = 100, width = 1200, height = 1500)
plot_heatmap(
  dataset = sc.sub,
  markers = unique(sel.markers$gene),
  anno_var = c("seurat_clusters", "sample"),
  anno_colors = list(
    alphabet2(16),
    rainbow(4)
  ),
  #               alphabet(8)),
  row_font_size = 10
)
dev.off()

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 2, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

jpeg("KRAS_model/results/immune.featureplots.jpeg", quality = 100, width = 1500, height = 1500)
FeaturePlot(sc.sub, features = sel.markers, ncol = 5)
dev.off()



# automated identification ------------------------------------------------

## ScType ----

# https://www.nature.com/articles/s41467-022-28803-w
# https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md


# load libraries and functions
lapply(c("dplyr", "Seurat", "HGNChelper", "ggraph", "igraph", "tidyverse", "data.tree"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# get cell-type-specific gene sets from our in-built database (DB)
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list <- gene_sets_prepare(db_, tissue)
lapply(gs_list, length)
head(gs_list[[1]])
head(gs_list[[2]])


# assign cell types
es.max <- sctype_score(
  scRNAseqData = sc.sub[["integrated"]]@scale.data,
  scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)
es.max[1:10, 1:5]

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(sc.sub@meta.data$seurat_clusters), function(cl) {
  es.max.cl <- sort(rowSums(es.max[, rownames(sc.sub@meta.data[sc.sub@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc.sub@meta.data$seurat_clusters == cl)), 10)
}))
sctype_scores <- cL_resutls %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"
print(sctype_scores[, 1:3])


# plot

sc.sub@meta.data$customclassif <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  sc.sub@meta.data$customclassif[sc.sub@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

DimPlot(sc.sub, reduction = "umap", label = TRUE, label.size = 4, label.box = T, repel = TRUE, group.by = "customclassif")
DimPlot(sc.sub, reduction = "umap", label = F, repel = TRUE, group.by = "customclassif", split.by = "sample")

head(sc.sub[])
saveRDS(sc.sub, file = "KRAS_model/objects/immune.rds")



# add GFP info ----

rm(list = ls())

sc <- readRDS(file = "KRAS_model/objects/integrated_annotated.rds")
head(sc) # 17721 cells
table(sc$sample)
DimPlot(sc, label = T, repel = T) + NoLegend()

# ZsGreen > 1 as GFP+
load(file = "~/Dropbox/BioInfo/Lab/TZones/KRAS_model/ZsGreen/ZsGreen.cells.prefilter.1.RData")
sel.cells.names.1 <- sel.cells.names
load(file = "~/Dropbox/BioInfo/Lab/TZones/KRAS_model/ZsGreen/ZsGreen.cells.prefilter.RData")
head(sel.cells.names)

int0 <- intersect(sel.cells.names, colnames(sc)) # 3542
int1 <- intersect(sel.cells.names.1, colnames(sc)) # 1586

p1 <- DimPlot(sc,
  cells.highlight = sel.cells.names,
  label = T, repel = T
) + NoLegend()
p2 <- DimPlot(sc,
  cells.highlight = sel.cells.names.1,
  label = T, repel = T
) + NoLegend()
p1 + p2

sc$GFP.0 <- "neg" # for ZsGreen > 0
Idents(sc) <- sc$GFP.0
sc <- SetIdent(sc, cells = int0, value = "pos")
DimPlot(sc)
sc$GFP.0 <- Idents(sc)

sc$GFP.1 <- "neg" # for ZsGreen > 1
Idents(sc) <- sc$GFP.1
sc <- SetIdent(sc, cells = int1, value = "pos")
DimPlot(sc)
sc$GFP.1 <- Idents(sc)

Idents(sc) <- sc$major.labels
table(sc$GFP.0, sc$GFP.1)


saveRDS(sc, file = "KRAS_model/objects/integrated_annotated.rds")



# end ---------------------------------------------------------------------
sessionInfo()
