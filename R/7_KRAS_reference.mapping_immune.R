# Reference mapping

# https://satijalab.org/seurat/articles/integration_mapping
# https://satijalab.org/seurat/articles/covid_sctmapping
# https://satijalab.org/seurat/articles/multimodal_reference_mapping

# labels are tansfered from Atlas
# next, all samples are integrated, and
# using Atlas as "Normal" for differential expression


# settings ----------------------------------------------------------------

rm(list = ls())

setwd()

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(BPCells)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(Nebulosa)
  library(hrbrthemes)
})

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

set.seed(100523)


# prepare data ----

ref <- readRDS("TZ_Atlas/data/normal_immune_reference.rds")
head(ref)
ref <- RunUMAP(ref, dims = 1:10, return.model = TRUE, reduction.name = "ref.umap")
DimPlot(ref)
DimPlot(ref, reduction = "ref.umap")


query <- readRDS(file = "KRAS_model/objects/immune.rds")
head(query)
DimPlot(query)
# DefaultAssay(query) <- "SCT" # this doesn't seem to affect the result


# label transfer from Atlas ----

anchors <- FindTransferAnchors(
  reference = ref, query = query,
  dims = 1:30, reference.reduction = "pca"
)
predictions1 <- TransferData(
  anchorset = anchors, refdata = ref$major.celltypes,
  dims = 1:30
)
predictions2 <- TransferData(
  anchorset = anchors, refdata = ref$minor.celltypes,
  dims = 1:30
)
query <- AddMetaData(query, metadata = predictions1$predicted.id, col.name = "Atlas.major.celltype")
query <- AddMetaData(query, metadata = predictions2$predicted.id, col.name = "Atlas.minor.celltype")
head(query)

d1 <- DimPlot(query,
  group.by = "Atlas.major.celltype",
  repel = T,
  alpha = 0.5, label = TRUE, ncol = 1, label.size = 4
) #+ NoLegend()
d2 <- DimPlot(query,
  group.by = "Atlas.minor.celltype",
  repel = T,
  alpha = 0.5, label = TRUE, ncol = 1, label.size = 4
) + NoLegend()
d3 <- DimPlot(query,
  group.by = "customclassif",
  repel = T,
  alpha = 0.5, label = TRUE, ncol = 1, label.size = 4
) #+ NoLegend()
d4 <- DimPlot(query,
  group.by = "major.labels",
  repel = T,
  alpha = 0.5, label = TRUE, ncol = 1, label.size = 4
) #+ NoLegend()
d1 + d2
d3 + d4
grid.arrange(d1, d2, d3, d4, ncol = 2)

table(query$predicted.id)

saveRDS(query, file = "KRAS_model/objects/immune_Atlas.anno.rds")

# inconsistency in the B cell cluster


# inspection ----

rm(list = ls())

sc <- readRDS(file = "KRAS_model/objects/immune_Atlas.anno.rds")
head(sc)

table(sc$Atlas.minor.celltype)
table(sc$Atlas.major.celltype, sc$sample)

DefaultAssay(sc) <- "SCT"
FeaturePlot(sc, features = c(
  "Cd3e", "Cd19", "Cd4", "Cd8a",
  "Foxp3", "Il2ra"
))

FeaturePlot(sc, c(
  "Trgv2", "Trgj2", "Trdc", "Trdv4", "Trdj1", "Trdmt1"
))



# cell type proportions ---------------------------------------------------

sc <- readRDS("KRAS_model/objects/immune_Atlas.anno.rds")

# label most abundant cell types for future analyses

Idents(sc) <- sc$Atlas.minor.celltype
Neutrophils <- WhichCells(sc, idents = "Neutrophils")
gd1 <- WhichCells(sc, idents = "gd cytotoxic T cells")
gd2 <- WhichCells(sc, idents = "gd17 T cells")

Idents(sc) <- sc$Atlas.major.celltype
sc <- SetIdent(sc, cells = Neutrophils, value = "Neutrophils")
sc <- SetIdent(sc, cells = gd1, value = "gd T cells")
sc <- SetIdent(sc, cells = gd2, value = "gd T cells")

table(Idents(sc))
table(sc$sample)
table(Idents(sc), sc$sample)

sc$new.cell.labels <- Idents(sc)

p1 <- DimPlot(sc, reduction = "umap", group.by = "sample")
p2 <- DimPlot(sc, reduction = "umap", group.by = "new.cell.labels", label = TRUE, repel = TRUE) +
  NoLegend()
p1 + p2

DimPlot(sc, split.by = "sample", group.by = "new.cell.labels")

Idents(sc) <- sc$sample
cells.n <- sample(WhichCells(sc, idents = "Normal"), 500)
cells.h <- sample(WhichCells(sc, idents = "Hyperplasia"), 500)
cells.d <- sample(WhichCells(sc, idents = "Dysplasia"), 500)
cells.c <- sample(WhichCells(sc, idents = "Carcinoma"), 500)

sc.subset <- subset(sc, cells = c(cells.n, cells.h, cells.d, cells.c))
table(sc.subset$sample)
DimPlot(sc.subset, split.by = "sample", group.by = "new.cell.labels", alpha = 0.5, pt.size = 1.5)



colors <- as.vector(glasbey(n = 24))
cluster.colors <- colors[1:20]
condition.colors <- c("brown", "orange", "yellow", "pink")

tab2 <- table(sc$sample, sc$new.cell.labels)
write.csv(tab2, file = "KRAS_model/results/immune/immune.major.celltype.proportions.csv")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("Cell_Types", "Sample", "Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = Cell_Types)) +
  geom_bar(position = "fill", stat = "identity") +
  #  scale_fill_manual(values = cluster.colors) +
  ggtitle("Major Cell Type Proportions") +
  theme_ipsum() +
  xlab("")

p2


saveRDS(sc, file = "KRAS_model/objects/immune_Atlas.anno.rds")





# Differential expression analysis ----

rm(list = ls())

output <- "KRAS_model/results/immune/DEGs/"

sc <- readRDS(file = "KRAS_model/objects/immune_Atlas.anno.rds")
head(sc)
table(sc$sample)
Idents(sc) <- sc$sample
DefaultAssay(sc) <- "integrated"
# DefaultAssay(sc) <- "SCT"
# sc <- PrepSCTFindMarkers(sc)


## global changes ----

comparisons <- list(
  Hyperplasia.vs.Normal = c("Hyperplasia", "Normal"),
  Dysplasia.vs.Normal = c("Dysplasia", "Normal"),
  Carcinoma.vs.Normal = c("Carcinoma", "Normal"),
  Dysplasia.vs.Hyperplasia = c("Dysplasia", "Hyperplasia"),
  Carcinoma.vs.Dysplasia = c("Carcinoma", "Dysplasia")
)

for (i in 1:length(comparisons)) {
  DEGs <- FindMarkers(sc, ident.1 = comparisons[[i]][1], ident.2 = comparisons[[i]][2])
  DEGs <- DEGs[DEGs$p_val_adj < 0.05, ]
  write.csv(DEGs, file = paste0(output, names(comparisons)[i], ".csv"))
}


## by cell type ----

sc$celltype.condition <- paste(sc$new.cell.labels, sc$sample, sep = "_")
Idents(sc) <- sc$celltype.condition
table(Idents(sc))

for (i in 1:length(comparisons)) {
  for (level in levels(factor(sc$new.cell.labels))) {
    try({
      ident1 <- paste0(level, "_", comparisons[[i]][1])
      ident2 <- paste0(level, "_", comparisons[[i]][2])
      degs <- FindMarkers(sc,
        ident.1 = ident1,
        ident.2 = ident2
      )
      write.csv(degs, file = paste0(output, "DEGs_", names(comparisons)[i], "_", level, ".csv"))
    })
    rm(ident1, ident2, degs)
  }
}


saveRDS(sc, file = "KRAS_model/objects/immune_Atlas.anno.rds")



# ssGSEA ------------------------------------------------------------------

rm(list = ls())

sc <- readRDS(file = "KRAS_model/objects/immune_Atlas.anno.rds")
head(sc)

library(escape)
library(dittoSeq)
library(UCell)
library(AUCell)

## global changes ----

# gene.sets <- getGeneSets(species = "Mus musculus", library = "H")
gene.sets <- getGeneSets(species = "Mus musculus", library = "C2") # , subcategory = "CP")
names(gene.sets)
length(names(gene.sets))
gene.sets <- gene.sets[grep("BIOCARTA_", names(gene.sets))]
# gene.sets <- gene.sets[grep("REACTOME_", names(gene.sets))]

ES <- enrichIt(
  obj = sc,
  method = "ssGSEA", # "Ucell",
  gene.sets = gene.sets,
  groups = 1000, cores = 12,
  min.size = 20
) # 20 for BioCarta and 50 for Reactome

head(ES)
sc2 <- AddMetaData(sc, ES)
head(sc2[[]])

ES2 <- data.frame(sc2[[]], Idents(sc2))
head(ES2)
colnames(ES2)[ncol(ES2)] <- "cluster"

PCA1 <- performPCA(enriched = ES2, gene.sets = names(gene.sets), groups = c("sample", "new.cell.labels"))
# pcaEnrichment(PCA1, PCx = "PC1", PCy = "PC2", contours = TRUE)
masterPCAPlot(ES2,
  gene.sets = names(gene.sets),
  PCx = "PC1", PCy = "PC2", top.contribution = 5
) +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(size = 16)
  ) +
  ggtitle(paste0("PC plot of BioCarta gene sets"))
ggsave(filename = paste0(output, "PCA_Biocarta.jpeg"))
#  ggtitle(paste0("PC plot of Reactome gene sets"))
# ggsave(filename = paste0(output, "PCA_Reactome.jpeg"))

enrich.scores <- getSignificance(ES2,
  group = "sample",
  gene.sets = names(gene.sets),
  fit = "ANOVA"
)
head(enrich.scores)
write.csv(enrich.scores, file = "KRAS_model/results/immune/pathways/BioCarta.csv")
# write.csv(enrich.scores, file = "KRAS_model/results/immune/pathways/Reactome.csv")



## new markers ----

# sc <- readRDS(file = "KRAS_model/objects/immune_Atlas.anno.rds")

DefaultAssay(sc) <- "integrated"
Idents(sc) <- sc$new.cell.labels
# sc <- PrepSCTFindMarkers(sc)

all.markers <- FindAllMarkers(sc, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]
head(all.markers)
table(all.markers$cluster)
tail(all.markers)

all.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "KRAS_model/results/immune/celltypes/celltype.markers.csv")

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

jpeg("KRAS_model/results/immune/celltypes/celltypes.heatmap.jpeg", quality = 100, width = 1200, height = 1500)
plot_heatmap(
  dataset = sc,
  markers = unique(sel.markers$gene),
  anno_var = c("new.cell.labels", "sample"),
  anno_colors = list(
    alphabet2(11),
    rainbow(5)
  ),
  #               alphabet(8)),
  row_font_size = 10
)
dev.off()

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 2, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

jpeg("KRAS_model/results/immune/celltypes/celltype.featureplots.jpeg", quality = 100, width = 1500, height = 1500)
FeaturePlot(sc, features = sel.markers, ncol = 5)
dev.off()


# add GFP info ----

sc <- readRDS(file = "KRAS_model/objects/immune_Atlas.anno.rds")
head(sc) # 5736 cells
table(sc$sample)
DimPlot(sc, group.by = "Atlas.major.celltype", label = T, repel = T) + NoLegend()

# ZsGreen > 1 as GFP+
load(file = "~/Dropbox/BioInfo/Lab/TZones/KRAS_model/ZsGreen/ZsGreen.cells.prefilter.1.RData")
sel.cells.names.1 <- sel.cells.names
load(file = "~/Dropbox/BioInfo/Lab/TZones/KRAS_model/ZsGreen/ZsGreen.cells.prefilter.RData")
head(sel.cells.names)

int0 <- intersect(sel.cells.names, colnames(sc)) # 86
int1 <- intersect(sel.cells.names.1, colnames(sc)) # 14

p1 <- DimPlot(sc,
  group.by = "Atlas.major.celltype", cells.highlight = sel.cells.names,
  label = T, repel = T
) + NoLegend()
p2 <- DimPlot(sc,
  group.by = "Atlas.major.celltype", cells.highlight = sel.cells.names.1,
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

Idents(sc) <- sc$Atlas.major.celltype
table(sc$GFP.0, sc$GFP.1)


saveRDS(sc, file = "KRAS_model/objects/immune_Atlas.anno.rds")



# end ---------------------------------------------------------------------

sessionInfo()
