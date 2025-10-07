
# KRAS model - subclustering within epithelial cells

# libraries ---------------------------------------------------------------

rm(list=ls())

setwd()

suppressPackageStartupMessages({
  library(Seurat);  library(sctransform);  library(SeuratDisk); #library(SeuratWrappers); library(SingleCellExperiment)
  library(dplyr);  library(stringr) #library(xlsx); 
  library(ggplot2); library(cowplot); library(viridis); library(gridExtra); library(RColorBrewer); library(patchwork); library(hrbrthemes); library(pals)
  library(SingleR); library(scRNAseq); library(scibet)
  library(slingshot, quietly = T)
  library(mclust, quietly = T); library(celldex); library(KernSmooth); library(dendextend); library(circlize)
  library(enrichR); library(clusterProfiler); library(enrichplot); library(DOSE); library(pathview); library(ReactomePA)
  library(ggpubr); library(eulerr); library(alluvial)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db); library(biomaRt)
  library(Scillus)
  library(snow)
  
})


set.seed(1)


# subsetting after re-integration ---------------------------------------------

cc.int.sct <- readRDS(file = "KRAS_model/objects/integrated_annotated.rds")
head(cc.int.sct[[]])

sc.sub <- subset(cc.int.sct, idents = "Stratified")
table(sc.sub$sample)
# Normal Hyperplasia.1 Hyperplasia.2     Dysplasia     Carcinoma 
# 151           121           183          2162          2661 
rm(cc.int.sct)

sc.sub <- SplitObject(sc.sub, split.by = "sample")
sc.sub <- lapply(sc.sub, SCTransform, method = "glmGamPoi")
features  <- SelectIntegrationFeatures(sc.sub, nfeatures = 3000)
sc.sub <- PrepSCTIntegration(sc.sub, anchor.features = features)
sc.sub <- lapply(X = sc.sub, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = sc.sub, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
sc.sub <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 50, dims = 1:30) # had to reduce k.weight probably because of low number of cells
sc.sub <- RunPCA(sc.sub)
sc.sub <- RunUMAP(sc.sub, reduction = "pca", dims = 1:10)
sc.sub <- FindNeighbors(sc.sub, dims = 1:30)
sc.sub <- FindClusters(sc.sub, resolution = 0.5)

sc.sub$sample <- gsub("Hyperplasia.1", "Hyperplasia", sc.sub$sample)
sc.sub$sample <- gsub("Hyperplasia.2", "Hyperplasia", sc.sub$sample)

sc.sub$sample <- factor(sc.sub$sample, 
                        levels = c("Normal",
                                   "Hyperplasia",
                                   "Dysplasia",
                                   "Carcinoma"))

p1 <- DimPlot(sc.sub, label = T, label.box = T, label.size = 5)
p2 <- DimPlot(sc.sub, group.by = "sample")

p1+p2

DimPlot(sc.sub, split.by = "sample")#, pt.size = 1)

table(sc.sub$sample)
# Normal Hyperplasia   Dysplasia   Carcinoma 
# 151         304        2162        2661 
table(Idents(sc.sub))


saveRDS(sc.sub, file = "KRAS_model/objects/stratified.rds")




## add GFP labels ----
## 

sc.sub <- readRDS("KRAS_model/objects/stratified.rds")

p0 <- DimPlot(sc.sub, label = T, label.box = T, repel = T, label.size = 6) +
  ggtitle("Seurat Clusters") + NoLegend()

load(file = "KRAS_model/ZsGreen/ZsGreen.cells.prefilter.RData")
head(sel.cells.names)
int1 <- intersect(sel.cells.names, colnames(sc.sub))

#DimPlot(cc.int.sct, label = T, label.box = T) + NoLegend()
p1 <- DimPlot(sc.sub, 
              label = T, label.box = T, repel = T, label.size = 6,
              cells.highlight = sel.cells.names, cols.highlight = "springgreen4", sizes.highlight = 0.5) + 
  ggtitle("ZsGreen+ cells") + NoLegend()
p2 <- plot_density(sc.sub, "Krt17")

p0+p1+p2

sc.sub$GFP <- "neg" # for ZsGreen > 0
Idents(sc.sub) <- sc.sub$GFP
sc.sub <- SetIdent(sc.sub, cells = int1, value = "pos")
p3 <- DimPlot(sc.sub) + ggtitle("ZsGreen+ [GFP > 0]")
sc.sub$GFP <- Idents(sc.sub)
table(sc.sub$GFP)


load(file = "KRAS_model/ZsGreen/ZsGreen.cells.prefilter.1.RData")
head(sel.cells.names)
int1 <- intersect(sel.cells.names, colnames(sc.sub))

sc.sub$GFP.1 <- "neg" # for ZsGreen > 1
Idents(sc.sub) <- sc.sub$GFP.1
sc.sub <- SetIdent(sc.sub, cells = int1, value = "pos")
sc.sub$GFP.1 <- Idents(sc.sub)
p4 <- DimPlot(sc.sub) + ggtitle("ZsGreen+ [GFP > 1]")
table(sc.sub$GFP.1)

Idents(sc.sub) <- sc.sub$seurat_clusters

p0+p2+p3+p4

head(sc.sub)


saveRDS(sc.sub, file = "KRAS_model/objects/stratified.rds")



# markers -----------------------------------------------------------------

options(Seurat.object.assay.version = "v3") # for backward compatibility with DoubletFinder

sc.sub <- readRDS("KRAS_model/objects/stratified.rds")

all.markers <- FindAllMarkers(sc.sub, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
table(all.markers$cluster)
tail(all.markers)

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "KRAS_model/results/stratified.markers.reclustering.csv")
#write.xlsx(all.markers, file="KRAS_model/results/stratified.markers.reclustering.xlsx")
#all.markers <- read.csv(file = "KRAS_model/results/stratified.markers.reclustering.csv", row.names = 1)

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

jpeg("KRAS_model/results/stratified.heatmap.jpeg", quality = 100, width = 1200, height = 1500)
plot_heatmap(dataset = sc.sub,  
             markers = unique(sel.markers$gene),
             anno_var = c("seurat_clusters","sample"),
             anno_colors = list(
               alphabet2(11),
               rainbow(5)),
#               alphabet(8)),
             row_font_size = 10
)
dev.off()

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 2, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

jpeg("KRAS_model/results/stratified.featureplots.jpeg", quality = 100, width = 1500, height = 1500)
FeaturePlot(sc.sub, features = sel.markers, ncol = 5)
dev.off()




# 

cells.9 <- WhichCells(sc.sub, idents = 9)
load("sc.all.celltypes.RData")
p1 <- DimPlot(cc.int.sct, cells.highlight = cells.9)
p2 <- DimPlot(cc.int.sct)
p2+p1
# Remove cluster 9 (rectum)
# This cluster presents several markers of rectum/glandular epithelium like Krt8, Krt19, 
# Mucine expression (Muc3, Muc13) and other secretory associated genes.
# Importantly, the different specific keratin of stratified epithelium (like Krt6, Krt1, Krt10, Krt14, Krt5 or Krt17) are not expressed in this cluster.

# not removed at this step, but mostly classified later as "unlabeled"
#sc.sub <- subset(sc.sub, ident = '9', invert = T)
#saveRDS(sc.sub, file = "stratified.rds")


# relabeling --------------------------------------------------------------

rm(list=ls())

sc.sub <- readRDS("KRAS_model/objects/stratified.rds")
Idents(sc.sub) <- sc.sub$seurat_clusters

# old stratified object labeled by Chloe (4/8/23)
sc.sub2 <- readRDS("KRAS_model/objects/old.stratified.rds")
head(sc.sub2)

Idents(sc.sub2)
sc.sub2$major.labels <- Idents(sc.sub2)
sc.sub$labels_Chloe <- "unlabeled"
Idents(sc.sub) <- sc.sub$labels_Chloe

for(i in 1:length(levels(sc.sub2$major.labels))){
  
  sub1 <- subset(sc.sub2, subset = major.labels == levels(sc.sub2$major.labels)[i])
  int1 <- intersect(colnames(sc.sub), colnames(sub1))
  sc.sub <- SetIdent(sc.sub, cells = int1, value = levels(sc.sub2$major.labels)[i])
  
  rm(sub1, int1)
}

table(Idents(sc.sub)) # 413 cells unlabeled
sc.sub$labels_Chloe <- Idents(sc.sub)
head(sc.sub)

d1 <- DimPlot(sc.sub, label = T, repel = T, label.size = 5) + ggtitle("Chloe's labels")
d2 <- DimPlot(sc.sub, label = T, group.by = "seurat_clusters", label.size = 6)
d1+d2

DimPlot(sc.sub, label = T, repel = T, label.size = 5)
DimPlot(sc.sub, label = F, repel = T, label.size = 5, split.by = "sample")



saveRDS(sc.sub, file = "KRAS_model/objects/stratified.rds")



# cell type proportions ---------------------------------------------------

sc.sub <- readRDS("KRAS_model/objects/stratified.rds")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:20]
condition.colors = c("brown","orange","yellow","pink")

DimPlot(sc.sub, label = F, repel = T, cols = cluster.colors)
DimPlot(sc.sub, label = F, repel = T, split.by = "sample", cols = cluster.colors)

tab2 <- table(sc.sub$sample, Idents(sc.sub))
write.csv(tab2, file = "KRAS_model/results/stratified/stratified.new.cluster.proportions.csv")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("New_clusters", "Sample","Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = New_clusters)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("New Cluster Proportions") +
  theme_ipsum() +
  xlab("")

p2


# Differential expression analysis ----

rm(list=ls())

output <- "KRAS_model/results/stratified/"

sc <- readRDS(file = "KRAS_model/objects/stratified.rds")
head(sc)
Idents(sc) <- sc$labels_Chloe
sc <- subset(sc, idents = "unlabeled", invert = T)
Idents(sc) <- sc$sample
DefaultAssay(sc) <- "integrated"
#DefaultAssay(sc) <- "SCT"
#sc <- PrepSCTFindMarkers(sc)


## global changes ----

comparisons <- list(Hyperplasia.vs.Normal = c("Hyperplasia","Normal"),
                    Dysplasia.vs.Normal = c("Dysplasia","Normal"),
                    Carcinoma.vs.Normal = c("Carcinoma","Normal"),
                    Dysplasia.vs.Hyperplasia = c("Dysplasia","Hyperplasia"),
                    Carcinoma.vs.Dysplasia = c("Carcinoma","Dysplasia"))

for(i in 1:length(comparisons)){
  DEGs <- FindMarkers(sc, ident.1 = comparisons[[i]][1], ident.2 = comparisons[[i]][2])
  DEGs <- DEGs[DEGs$p_val_adj < 0.05, ]
  write.csv(DEGs, file = paste0(output, names(comparisons)[i], ".csv"))
}


## by cell type ----

sc$celltype.condition <- paste(sc$labels_Chloe, sc$sample, sep="_")
Idents(sc) <- sc$celltype.condition
table(Idents(sc))

for(i in 1:length(comparisons)){
  for (level in levels(factor(sc$labels_Chloe))){
    try({
      ident1 <- paste0(level, "_", comparisons[[i]][1])
      ident2 <- paste0(level, "_", comparisons[[i]][2])
      degs <- FindMarkers(sc, 
                          ident.1 = ident1, 
                          ident.2 = ident2)
      write.csv(degs, file=paste0(output,"DEGs/DEGs_", names(comparisons)[i], "_", level,".csv"))
    })
    rm(ident1, ident2, degs)
  }
}




# ssGSEA ------------------------------------------------------------------

rm(list=ls())

sc <- readRDS(file = "KRAS_model/objects/stratified.rds")
head(sc)

library(escape); library(dittoSeq)
library(UCell); library(AUCell)

## global changes ----

#gene.sets <- getGeneSets(species = "Mus musculus", library = "H")
gene.sets <- getGeneSets(species = "Mus musculus", library = "C2")#, subcategory = "CP")
names(gene.sets)
length(names(gene.sets))
#gene.sets <- gene.sets[grep("BIOCARTA_", names(gene.sets))]
gene.sets <- gene.sets[grep("REACTOME_", names(gene.sets))]

ES <- enrichIt(obj = sc, 
               method = "ssGSEA", # "Ucell",
               gene.sets = gene.sets, 
               groups = 1000, cores = 12, 
               min.size = 50) # 20 for BioCarta and 50 for Reactome

head(ES)
sc2 <- AddMetaData(sc, ES)
head(sc2[[]])

ES2 <- data.frame(sc2[[]], Idents(sc2))
head(ES2)
colnames(ES2)[ncol(ES2)] <- "cluster"

PCA1 <- performPCA(enriched = ES2, gene.sets = names(gene.sets), groups = c("sample", "labels_Chloe"))
#pcaEnrichment(PCA1, PCx = "PC1", PCy = "PC2", contours = TRUE)
masterPCAPlot(ES2, gene.sets = names(gene.sets),
              PCx = "PC1", PCy = "PC2", top.contribution = 5) +
  theme(text=element_text(size=14),
        plot.title = element_text(size = 16)) +
#  ggtitle(paste0("PC plot of BioCarta gene sets"))
#ggsave(filename = paste0(output, "PCA_Biocarta.jpeg"))
  ggtitle(paste0("PC plot of Reactome gene sets"))
ggsave(filename = paste0(output, "PCA_Reactome.jpeg"))

enrich.scores <- getSignificance(ES2, 
                                 group = "sample", 
                                 gene.sets = names(gene.sets),
                                 fit = "ANOVA")
head(enrich.scores)
#write.csv(enrich.scores, file = "KRAS_model/results/stratified/BioCarta.csv")
write.csv(enrich.scores, file = "KRAS_model/results/stratified/Reactome.csv")



## new markers ----

sc <- readRDS(file = "KRAS_model/objects/stratified.rds")
DefaultAssay(sc) <- "integrated"
#sc <- PrepSCTFindMarkers(sc)

all.markers <- FindAllMarkers(sc, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
table(all.markers$cluster)
tail(all.markers)

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "KRAS_model/results/stratified/celltype.markers.csv")

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

jpeg("KRAS_model/results/stratified/celltypes.heatmap.jpeg", quality = 100, width = 1200, height = 1500)
plot_heatmap(dataset = sc,  
             markers = unique(sel.markers$gene),
             anno_var = c("labels_Chloe","sample"),
             anno_colors = list(
               alphabet2(11),
               rainbow(5)),
             #               alphabet(8)),
             row_font_size = 10
)
dev.off()

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 2, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

jpeg("KRAS_model/results/stratified/celltype.featureplots.jpeg", quality = 100, width = 1500, height = 1500)
FeaturePlot(sc, features = sel.markers, ncol = 5)
dev.off()




# end ---------------------------------------------------------------------
sessionInfo()
