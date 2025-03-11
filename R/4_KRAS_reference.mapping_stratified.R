
# Reference mapping

# https://satijalab.org/seurat/articles/integration_mapping
# https://satijalab.org/seurat/articles/covid_sctmapping
# https://satijalab.org/seurat/articles/multimodal_reference_mapping

# labels are tansfered from Atlas
# next, all samples are integrated, and
# using Atlas as "Normal" for differential expression


# settings ----------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/TZones/")

suppressPackageStartupMessages({
  library(Seurat); library(SeuratDisk)
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

ref <- readRDS("TZ_Atlas/data/integrated.normal.reference.KRAS.DKO.curated.rds")
head(ref)
ref$cell_type <- Idents(ref)
ref <- RunUMAP(ref, dims = 1:10, return.model=TRUE, reduction.name = "ref.umap")
DimPlot(ref)
DimPlot(ref, reduction = "ref.umap")


query <- readRDS(file = "KRAS_model/objects/stratified.rds")
head(query)
DimPlot(query)
query$precell_type <- Idents(query)



# label transfer from Atlas ----

anchors <- FindTransferAnchors(reference = ref, query = query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$cell_type,
                            dims = 1:30)
object <- AddMetaData(query, metadata = predictions)
head(object)

d1 <- DimPlot(object, group.by = "predicted.id", 
              repel = T,
              alpha = 0.5, label = TRUE, ncol = 1, label.size = 4) #+ NoLegend()
d2 <- DimPlot(object, group.by = "precell_type", 
              repel = T,
              alpha = 0.5, label = TRUE, ncol = 1, label.size = 4) #+ NoLegend()
d1+d2

table(object$precell_type)
table(object$predicted.id)
table(object$precell_type, object$predicted.id)

FeaturePlot(object, "prediction.score.AC.basal.proximal", cols = c("lightgrey", "darkred"))
plot_density(object, "prediction.score.AC.basal.proximal")


object$cell_type <- object$predicted.id
ref$precell_type <- ref$cell_type
object$condition <- object$sample
ref$condition <- "Atlas"


# integration with Atlas ----

DefaultAssay(ref) <- "SCT"

# remove contaminating rectum cells in query
Idents(object) <- object$cell_type
object <- subset(object, idents = c("Anal gland", "Fibroblasts", "Rectum", "Hair follicle"), invert = T)
table(object$predicted.id)

# remove cell types absent in query
ref <- subset(ref, idents = c("Anal gland", "Fibroblasts", "Rectum", "Hair follicle"), invert = T)
table(ref$cell_type)

ref.sub <- SplitObject(ref, split.by = "sample")
obj.sub <- SplitObject(object, split.by = "sample") # remove Normal sample, already in Atlas

reference.list <- c(ref.sub, obj.sub[-1])
names(reference.list)

reference.list <- lapply(X = reference.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = features)
reference.list <- lapply(X = reference.list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
int.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
table(int.sct$sample)
table(int.sct$condition)

int.sct <- RunPCA(int.sct) %>%
  FindNeighbors(dims = 1:10) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

head(int.sct)


saveRDS(int.sct, file = "KRAS_model/objects/stratified_Atlas.rds")




# inspection ----

sc <- readRDS(file = "KRAS_model/objects/stratified_Atlas.rds")
head(sc)

p1 <- DimPlot(sc, reduction = "umap", group.by = "sample")
p2 <- DimPlot(sc, reduction = "umap", group.by = "condition")
p3 <- DimPlot(sc, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  NoLegend()
p1 + p2 + p3

table(sc$condition)
table(sc$sample, sc$cell_type)

DimPlot(sc, split.by = "sample", group.by = "cell_type")
sc$condition <- factor(sc$condition, levels = c("Atlas","Hyperplasia","Dysplasia","Carcinoma"))
DimPlot(sc, split.by = "condition", group.by = "cell_type")

Idents(sc) <- sc$condition
cells.a <- sample(WhichCells(sc, idents = "Atlas"), 250)
cells.h <- sample(WhichCells(sc, idents = "Hyperplasia"), 250)
cells.d <- sample(WhichCells(sc, idents = "Dysplasia"), 250)
cells.c <- sample(WhichCells(sc, idents = "Carcinoma"), 250)

sc.subset <- subset(sc, cells = c(cells.a, cells.h, cells.d, cells.c))
table(sc.subset$condition)
DimPlot(sc.subset, split.by = "condition", group.by = "cell_type", alpha = 0.5, pt.size = 1.5)


saveRDS(sc, file = "KRAS_model/objects/stratified_Atlas.rds")



# Differential composition analysis ----

df_comp <- as.data.frame.matrix(table(sc$condition, sc$cell_type))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

# proportions

tab2 <- table(sc$condition, sc$cell_type)
write.csv(tab2, file = "KRAS_model/results/stratified/Atlas/cluster.proportions.csv")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("Cell_type", "Condition","Proportion")
head(tab2b)

ggplot(tab2b, aes(x = Condition, y = Proportion, fill = Cell_type)) +
  geom_bar(position = "fill", stat = "identity") +
  #  scale_fill_manual(values = cluster.colors) +
  ggtitle("Atlas-transfered labels") +
  theme_ipsum() +
  xlab("")





# Differential expression analysis ----

rm(list=ls())

sc <- readRDS(file = "KRAS_model/objects/stratified_Atlas.rds")
head(sc)
Idents(sc) <- sc$condition
DefaultAssay(sc) <- "integrated"
#DefaultAssay(sc) <- "SCT"
#sc <- PrepSCTFindMarkers(sc)


## global changes ----

comparisons <- list(Hyperplasia.vs.Atlas = c("Hyperplasia","Atlas"),
                    Dysplasia.vs.Atlas = c("Dysplasia","Atlas"),
                    Carcinoma.vs.Atlas = c("Carcinoma","Atlas"),
                    Dysplasia.vs.Hyperplasia = c("Dysplasia","Hyperplasia"),
                    Carcinoma.vs.Dysplasia = c("Carcinoma","Dysplasia"))

for(i in 1:length(comparisons)){
  DEGs <- FindMarkers(sc, ident.1 = comparisons[[i]][1], ident.2 = comparisons[[i]][2])
  DEGs <- DEGs[DEGs$p_val_adj < 0.05, ]
  write.csv(DEGs, file = paste0("KRAS_model/results/stratified/Atlas/", names(comparisons)[i], ".csv"))
}


## by cell type ----

output <- "KRAS_model/results/stratified/Atlas/DEGs/"
sc$celltype.condition <- paste(sc$cell_type, sc$condition, sep="_")
Idents(sc) <- sc$celltype.condition
table(Idents(sc))

for(i in 1:length(comparisons)){
  for (level in levels(factor(sc$cell_type))){
    try({
      ident1 <- paste0(level, "_", comparisons[[i]][1])
      ident2 <- paste0(level, "_", comparisons[[i]][2])
      degs <- FindMarkers(sc, 
                          ident.1 = ident1, 
                          ident.2 = ident2)
      write.csv(degs, file=paste0(output,"DEGs_", names(comparisons)[i], "_", level,".csv"))
    })
    rm(ident1, ident2, degs)
  }
}
  



# ssGSEA ------------------------------------------------------------------

rm(list=ls())

sc <- readRDS(file = "KRAS_model/objects/stratified_Atlas.rds")
head(sc)

library(escape); library(dittoSeq)
library(UCell); library(AUCell)

## global changes ----

#gene.sets <- getGeneSets(species = "Mus musculus", library = "H")
gene.sets <- getGeneSets(species = "Mus musculus", library = "C2")#, subcategory = "CP")
names(gene.sets)
length(names(gene.sets))
gene.sets <- gene.sets[grep("BIOCARTA_", names(gene.sets))]
#gene.sets <- gene.sets[grep("REACTOME_", names(gene.sets))]

ES <- enrichIt(obj = sc, 
               method = "ssGSEA", # "Ucell",
               gene.sets = gene.sets, 
               groups = 1000, cores = 12, 
               min.size = 20) # 

head(ES)
sc2 <- AddMetaData(sc, ES)
head(sc2[[]])

ES2 <- data.frame(sc2[[]], Idents(sc2))
head(ES2)
colnames(ES2)[ncol(ES2)] <- "cluster"

enrich.scores <- getSignificance(ES2, 
                          group = "condition", 
                          gene.sets = names(gene.sets),
                          fit = "ANOVA")
head(enrich.scores)
write.csv(enrich.scores, file = "KRAS_model/results/stratified/Atlas/BioCarta.csv")
#write.csv(enrich.scores, file = "KRAS_model/results/stratified/Atlas/Reactome.csv")



## by cell type ----

rm(list=ls())

sc <- readRDS(file = "KRAS_model/objects/stratified_Atlas.rds")
head(sc)

library(escape); library(dittoSeq)
library(UCell); library(AUCell)

output <- "KRAS_model/results/stratified/Atlas/pathways/"
Idents(sc) <- sc$cell_type
table(Idents(sc))

gene.sets <- getGeneSets(species = "Mus musculus", library = "C2")#, subcategory = "CP")
biocarta <- gene.sets[grep("BIOCARTA_", names(gene.sets))]
reactome <- gene.sets[grep("REACTOME_", names(gene.sets))]

for (level in levels(factor(sc$cell_type))){
  try({
    sc.subset <- subset(sc, idents = level)

    ES.biocarta <- enrichIt(obj = sc.subset, 
                   method = "ssGSEA", # "Ucell",
                   gene.sets = biocarta, 
                   groups = 1000, cores = 12, 
                   min.size = 20) #
    sc.subset.1 <- AddMetaData(sc.subset, ES.biocarta)
    ES.biocarta <- data.frame(sc.subset.1[[]], Idents(sc.subset.1))
    output.biocarta <- getSignificance(ES.biocarta, 
                          group = "condition", 
                          gene.sets = names(biocarta),
                          fit = "ANOVA")
    write.csv(output.biocarta, file = paste0(output,"BioCarta_", level,".csv"))
    
    PCA1 <- performPCA(enriched = ES.biocarta, gene.sets = names(biocarta), groups = c("condition", "cell_type"))
    masterPCAPlot(ES.biocarta, gene.sets = names(biocarta),
                  PCx = "PC1", PCy = "PC2", top.contribution = 5) +
      theme(text=element_text(size=14),
            plot.title = element_text(size = 16)) +
      ggtitle(paste0("PC plot of BioCarta gene sets in ", level))
    ggsave(filename = paste0(output, "PCA_Biocarta_", level, ".jpeg"))
    
    ES.reactome <- enrichIt(obj = sc.subset, 
                            method = "ssGSEA", # "Ucell",
                            gene.sets = reactome, 
                            groups = 1000, cores = 12, 
                            min.size = 50) #
    sc.subset.2 <- AddMetaData(sc.subset, ES.reactome)
    ES.reactome <- data.frame(sc.subset.2[[]], Idents(sc.subset.2))
    output.reactome <- getSignificance(ES.reactome, 
                                   group = "condition", 
                                   gene.sets = names(reactome),
                                   fit = "ANOVA")
    write.csv(output.reactome, file = paste0(output,"Reactome_", level,".csv"))
    
    PCA2 <- performPCA(enriched = ES.reactome, gene.sets = names(reactome), groups = c("condition", "cell_type"))
    masterPCAPlot(ES.reactome, gene.sets = names(reactome),
                  PCx = "PC1", PCy = "PC2", top.contribution = 5) +
      theme(text=element_text(size=14),
            plot.title = element_text(size = 16)) +
      ggtitle(paste0("PC plot of Reactome gene sets in ", level))
    ggsave(filename = paste0(output, "PCA_Reactome_", level, ".jpeg"))
    
  })
  
  rm(sc.subset, sc.subset.1, sc.subset.2, ES.biocarta, ES.reactome, PCA1, PCA2)

}





# end ---------------------------------------------------------------------

sessionInfo()



