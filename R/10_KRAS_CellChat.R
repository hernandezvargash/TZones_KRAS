
# KRAS model - Cell-Cell communication for major cell types

# summary -----------------------------------------------------------------

# Ligand-Receptor Interaction using CellChat: http://www.cellchat.org/
# https://www.nature.com/articles/s41467-021-21246-9


# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/TZones/")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(Nebulosa)
  library(xlsx)
  library(viridis)
  library(CellChat)
  library(patchwork)
  library(NMF)
  library(ggalluvial)
  library(pals)
  
})

set.seed(1)


# full integrated KRAS object ------------------------------------------------------------

sc.all <- readRDS(file = "KRAS_model/objects/integrated_annotated.rds")

head(sc.all[]) # 17721 cells
DefaultAssay(sc.all) <- "integrated"
Idents(sc.all)
DimPlot(sc.all)


# add labels of interest

stratified <- readRDS(file = "KRAS_model/objects/stratified_Atlas.rds")
stratified$cell_type <- as.factor(stratified$cell_type)
Idents(stratified) <- stratified$cell_type
DimPlot(stratified)
immune <- readRDS(file = "KRAS_model/objects/immune_Atlas.anno.rds")
Idents(immune) <- immune$new.cell.labels
DimPlot(immune)

for(i in 1:length(levels(stratified$cell_type))){
  
  sub1 <- subset(stratified, subset = cell_type == levels(stratified$cell_type)[i])
  int1 <- intersect(colnames(sc.all), colnames(sub1))
  sc.all <- SetIdent(sc.all, cells = int1, value = levels(stratified$cell_type)[i])
  
  rm(sub1, int1)
}

for(i in 1:length(levels(immune$new.cell.labels))){
  
  sub1 <- subset(immune, subset = new.cell.labels == levels(immune$new.cell.labels)[i])
  int1 <- intersect(colnames(sc.all), colnames(sub1))
  sc.all <- SetIdent(sc.all, cells = int1, value = levels(immune$new.cell.labels)[i])
  
  rm(sub1, int1)
}

table(Idents(sc.all))
DimPlot(sc.all)
sc.all$cellchat.label <- Idents(sc.all)

sc.all$condition <- sc.all$sample
sc.all$condition <- gsub("Hyperplasia.1", "Hyperplasia", sc.all$condition)
sc.all$condition <- gsub("Hyperplasia.2", "Hyperplasia", sc.all$condition)
sc.all$condition <- factor(sc.all$condition, levels = c("Normal", "Hyperplasia", "Dysplasia", "Carcinoma"))

table(sc.all$condition, sc.all$cellchat.label)

# group cell types to have enough cells for cellchat per condition

sc.all <- RenameIdents(sc.all,
                       "TZ suprabasal" = "TZ",
                       "TZ basal" = "TZ",
                       "AC suprabasal proximal" = "AC suprabasal",
                       "AC suprabasal distal" = "AC suprabasal",
                       "AC basal proximal" = "AC basal",
                       "AC basal distal" = "AC basal")
sc.all$cellchat.label <- Idents(sc.all)
table(sc.all$condition, sc.all$cellchat.label)

# some stratified were removed after stratified reclustering analysis
# also remove cell types with few cells

sc.all <- subset(sc.all, idents = c("B cells", "Stratified"), invert = T)
sc.all$cellchat.label <- Idents(sc.all)
table(sc.all$condition, sc.all$cellchat.label)

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:20]
condition.colors = c("pink","yellow","orange","brown")

p1 <- DimPlot(sc.all, pt.size = 0.3, cols = cluster.colors) + ggtitle("CellChat Cell Type")  +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(sc.all, pt.size = 0.2, cols = condition.colors, group.by = "condition") + ggtitle("Condition")
p1+p2

DimPlot(sc.all, pt.size = 0.5, cols = condition.colors, split.by = "condition", group.by = "condition") + ggtitle(NULL)
DimPlot(sc.all, pt.size = 0.5, cols = cluster.colors, split.by = "condition")


saveRDS(sc.all, file = "KRAS_model/objects/cellchat.rds")

DefaultAssay(sc.all) <- "RNA"



# CellChat major cell types ------------------------------------------

rm(list=ls())
sc.all <- readRDS(file = "KRAS_model/objects/cellchat.rds")
sc.all$cellchat.label <- Idents(sc.all)
table(sc.all$condition, sc.all$cellchat.label)
DefaultAssay(sc.all) <- "RNA"

cellchat <- createCellChat(object = sc.all, group.by = "cellchat.label")
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
#future::plan("multisession", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

future::plan("multisession", workers = 12)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10) # there are only 13 Neutrophils in Normal condition
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

jpeg("KRAS_model/results/cellchat/netVisual_circle_all.cells.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  jpeg(paste0("KRAS_model/results/cellchat/netVisual_circle", rownames(mat)[i], ".jpeg"), height = 900, width = 900, quality = 100)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

jpeg("KRAS_model/results/cellchat/netAnalysis_heatmap_all.cells.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 10, height = 25)
dev.off()


df.net <- subsetCommunication(cellchat)
head(df.net)
df.net.2 <- subsetCommunication(cellchat, slot.name = "netP") # signaling pathways
head(df.net.2)

write.csv(df.net, file = "KRAS_model/results/cellchat/cellchat.net.genes.csv")
write.csv(df.net.2, file = "KRAS_model/results/cellchat/cellchat.net.pathways.csv")


save(cellchat, file = "KRAS_model/objects/cellchat_sc.all.RData")





# CellChat by condition ------------------------------------------

rm(list=ls())
sc.all <- readRDS(file = "KRAS_model/objects/cellchat.rds")
sc.all$cellchat.label <- Idents(sc.all)
table(sc.all$condition, sc.all$cellchat.label)
DefaultAssay(sc.all) <- "RNA"

conditions <- levels(sc.all$condition)
table(sc.all$condition)

for(i in 1:length(conditions)){
  sct.temp <- subset(sc.all, subset = condition == conditions[i])
  DefaultAssay(sct.temp) <- "RNA"
  
  future::plan("multisession", workers = 1)
  cellchat <- createCellChat(object = sct.temp, group.by = "cellchat.label")
  cellchat@DB <- CellChatDB.mouse
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  future::plan("multisession", workers = 12)
  options(future.globals.maxSize = 8000 * 1024^2)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  
  jpeg(paste0("KRAS_model/results/cellchat/netVisual_circle_", conditions[i], ".jpeg"), height = 900, width = 900, quality = 100)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  mat <- cellchat@net$weight
  
  for (j in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[j, ] <- mat[j, ]
    
    jpeg(paste0("KRAS_model/results/cellchat/netVisual_circle_", conditions[i], "_", rownames(mat)[j], ".jpeg"), height = 900, width = 900, quality = 100)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j])
    dev.off()
  }
  
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  jpeg(paste0("KRAS_model/results/cellchat/netAnalysis_heatmap_", conditions[i], ".jpeg"), height = 1161, width = 1065, quality = 100)
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 10, height = 25) +
    netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 10, height = 25)
  dev.off()
  
  df.net <- subsetCommunication(cellchat)
  df.net.2 <- subsetCommunication(cellchat, slot.name = "netP") # signaling pathways

  write.csv(df.net, file = paste0("KRAS_model/results/cellchat/cellchat.net.genes_", conditions[i], ".csv"))
  write.csv(df.net.2, file = paste0("KRAS_model/results/cellchat/cellchat.net.pathways_", conditions[i], ".csv"))

  save(cellchat, file = paste0("KRAS_model/objects/cellchat_", conditions[i], ".RData"))
  
  rm(sct.temp, cellchat, mat, mat2, df.net, df.net.2)
  
}


# Normal vs All: major clusters ----------------------------

rm(list=ls())


load("~/Dropbox/BioInfo/Lab/TZones/KRAS_model/objects/cellchat/cellchat_Normal.RData")
cellchat.1 <- cellchat
load("~/Dropbox/BioInfo/Lab/TZones/KRAS_model/objects/cellchat/cellchat_Hyperplasia.RData")
cellchat.2 <- cellchat
load("~/Dropbox/BioInfo/Lab/TZones/KRAS_model/objects/cellchat/cellchat_Dysplasia.RData")
cellchat.3 <- cellchat
load("~/Dropbox/BioInfo/Lab/TZones/KRAS_model/objects/cellchat/cellchat_Carcinoma.RData")
cellchat.4 <- cellchat

object.list <- list(Normal = cellchat.1, Hyperplasia = cellchat.2, Dysplasia = cellchat.3, Carcinoma = cellchat.4)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Compare the number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4), measure = "count",
                           size.text = 16, 
                           title.name = "Number of Interactions")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4), measure = "weight",
                           size.text = 16, 
                           title.name = "Interaction Strength")

gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", comparison = c(1, 2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(1, 2))
# The width of edges represent the relative number of interactions or interaction strength. 
# Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", comparison = c(1, 3))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(1, 3),)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", comparison = c(1, 4))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(1, 4))

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", comparison = c(2, 3))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(2, 3),)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", comparison = c(3, 4))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(3, 4))

gg1 <- netVisual_heatmap(cellchat, measure = "count", comparison = c(1, 2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(1, 2))
gg1 + gg2
gg1 <- netVisual_heatmap(cellchat, measure = "count", comparison = c(1, 3))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(1, 3))
gg1 + gg2
gg1 <- netVisual_heatmap(cellchat, measure = "count", comparison = c(1, 4))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(1, 4))
gg1 + gg2

gg1 <- netVisual_heatmap(cellchat, measure = "count", comparison = c(2, 3))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(2, 3))
gg1 + gg2
gg1 <- netVisual_heatmap(cellchat, measure = "count", comparison = c(3, 4))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(3, 4))
gg1 + gg2

# To better control the node size and edge weights of the inferred networks across different datasets, 
# we compute the maximum number of cells per cell group and the maximum number of interactions 
# (or interaction weights) across all datasets.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,4), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


# Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

for(i in levels(cellchat@idents$joint)){
    gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = i, comparison = c(1, 2))
    gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = i, comparison = c(1, 3))
    gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = i, comparison = c(1, 4))
    gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = i, comparison = c(2, 3))
    gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = i, comparison = c(3, 4))
    patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4,gg5))
    ggsave(filename = paste0("KRAS_model/results/cellchat/signaling.changes_", i, ".jpeg"))
}



# Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", min_dist = 0.2) # Sensible values are in the range 0.001 to 0.5
cellchat <- netClustering(cellchat, type = "functional", k=10)
netVisual_embeddingPairwise(cellchat, type = "functional", 
                                   dot.size = c(2, 10), label.size = 3, 
                                   title = "Functional Similarity All")


comparison.list <- list(HvsN = c(1,2),
                        DvsN = c(1,3),
                        CvsN = c(1,4),
                        DvsH = c(2,3),
                        CvsD = c(3,4))
names(comparison.list)
for(i in 1:length(comparison.list)){
  cellchat.temp <- computeNetSimilarityPairwise(cellchat, type = "functional", comparison = comparison.list[[i]])
  cellchat.temp <- netEmbedding(cellchat.temp, type = "functional", comparison = comparison.list[[i]], min_dist = 0.2) # Sensible values are in the range 0.001 to 0.5
  cellchat.temp <- netClustering(cellchat.temp, type = "functional", k=10, comparison = comparison.list[[i]])
  gg1 <- netVisual_embeddingPairwise(cellchat.temp, 
                                     comparison = comparison.list[[i]],
                                     type = "functional", 
                                     dot.size = c(2, 10), label.size = 3, 
                                     title = paste0("Functional Similarity ", names(comparison.list)[i]))
  gg2 <- rankSimilarity(cellchat.temp, type = "functional", comparison1 = comp1)
  gg1+gg2
  ggsave(paste0("KRAS_model/results/cellchat/functional.similarity_", names(comparison.list)[i], ".jpeg"))
  rm(cellchat.temp, gg1, gg2)
}


# Compare the overall information flow of each signaling pathway
# 
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2,3,4))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,2,3,4))
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,2))
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,3))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,3))
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,4))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,4))
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(2,3))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(2,3))
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(3,4))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(3,4))
gg1 + gg2



# Compare outgoing (or incoming) signaling

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pathway.union <- union(pathway.union, object.list[[i+2]]@netP$pathways)
pathway.union <- union(pathway.union, object.list[[i+3]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 25)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25)
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25)
draw(ht1 + ht2 + ht3 + ht4, ht_gap = unit(2, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 25, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25, color.heatmap = "GnBu")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25, color.heatmap = "GnBu")
draw(ht1 + ht2 + ht3 + ht4, ht_gap = unit(2, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 25, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25, color.heatmap = "OrRd")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 25, color.heatmap = "OrRd")
draw(ht1 + ht2 + ht3 + ht4, ht_gap = unit(2, "cm"))


# Identify the upgulated and down-regulated signaling ligand-receptor pairs

netVisual_bubble(cellchat, 
                 sources.use = 5, 
                 targets.use = 1:3,  
                 comparison = c(1, 2, 3, 4), 
                 font.size = 16,
                 angle.x = 45)

netVisual_bubble(cellchat, 
                 sources.use = 6, 
                 targets.use = 1:3,  
                 comparison = c(1, 2, 3, 4), 
                 font.size = 16,
                 angle.x = 45)

netVisual_bubble(cellchat, 
                 sources.use = 7, 
                 targets.use = 1:3,  
                 comparison = c(1, 2, 3, 4), 
                 font.size = 16,
                 angle.x = 45)

netVisual_bubble(cellchat, 
                 sources.use = 8, 
                 targets.use = 1:3,  
                 comparison = c(1, 2, 3, 4), 
                 font.size = 16,
                 angle.x = 45)

netVisual_bubble(cellchat, 
                 sources.use = 1:3, 
                 targets.use = 5,  
                 comparison = c(1, 2, 3, 4), 
                 font.size = 16,
                 angle.x = 45)

netVisual_bubble(cellchat, 
                 sources.use = 1:3, 
                 targets.use = 6,  
                 comparison = c(1, 2, 3, 4), 
                 font.size = 16,
                 angle.x = 45)

netVisual_bubble(cellchat, 
                 sources.use = 1:3, 
                 targets.use = 7,  
                 comparison = c(1, 2, 3, 4), 
                 font.size = 16,
                 angle.x = 45)

netVisual_bubble(cellchat, 
                 sources.use = 1:3, 
                 targets.use = 8,  
                 comparison = c(1, 2, 3, 4), 
                 font.size = 16,
                 angle.x = 45)


# Identify dysfunctional signaling by using differential expression analysis

# The above method for identifying the upregulated and down-regulated signaling is performed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups.
# Alternative, we can identify the upregulated and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis. 
# Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, 
# and then obtain the upregulated and down-regulated signaling based on the fold change of ligands in the sender cells and receptors in the receiver cells. Such analysis can be done as follows.


## N vs H ----

object.list <- list(Normal = cellchat.1, Hyperplasia = cellchat.2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Hyperplasia"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Hyperplasia",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Normal",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 5:8, targets.use = 1:3, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 5:8, targets.use = 1:3, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1:3, targets.use = 5:8, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 1:3, targets.use = 5:8, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2


# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 5:8, targets.use = 1:3, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 5:8, targets.use = 1:3, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 1:3, targets.use = 5:8, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 1:3, targets.use = 5:8, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# visualize the enriched ligands in the first condition
net.down.2 <- net.down[is.finite(net.down$receptor.logFC), ]
computeEnrichmentScore(net.down.2, species = 'mouse')
# visualize the enriched ligands in the second condition
net.up.2 <- net.up[is.finite(net.up$receptor.logFC), ]
net.up.2 <- net.up.2[is.finite(net.up.2$ligand.logFC), ]
net.up.2 <- na.omit(net.up.2)
computeEnrichmentScore(net.up.2, species = 'mouse')

# Compare the signaling gene expression distribution between different datasets
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Normal", "Hyperplasia")) # set factor level
plotGeneExpression(cellchat, signaling = "TNF", split.by = "datasets", colors.ggplot = T)



## N vs D ----

object.list <- list(Normal = cellchat.1, Dysplasia = cellchat.3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pos.dataset = "Dysplasia"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Dysplasia",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Normal",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 5:8, targets.use = 1:3, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 5:8, targets.use = 1:3, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1:3, targets.use = 5:8, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 1:3, targets.use = 5:8, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 5:8, targets.use = 1:3, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 5:8, targets.use = 1:3, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 1:3, targets.use = 5:8, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 1:3, targets.use = 5:8, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

net.down.2 <- net.down[is.finite(net.down$receptor.logFC), ]
computeEnrichmentScore(net.down.2, species = 'mouse')
net.up.2 <- net.up[is.finite(net.up$receptor.logFC), ]
net.up.2 <- net.up.2[is.finite(net.up.2$ligand.logFC), ]
#net.up.2 <- na.omit(net.up.2)
computeEnrichmentScore(net.up.2, species = 'mouse')

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Normal", "Dysplasia")) # set factor level
plotGeneExpression(cellchat, signaling = "GRN", split.by = "datasets", colors.ggplot = T)



## N vs C ----

object.list <- list(Normal = cellchat.1, Carcinoma = cellchat.4)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pos.dataset = "Carcinoma"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Carcinoma",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Normal",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 5:8, targets.use = 1:3, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 5:8, targets.use = 1:3, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1:3, targets.use = 5:8, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 1:3, targets.use = 5:8, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 5:8, targets.use = 1:3, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 5:8, targets.use = 1:3, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 1:3, targets.use = 5:8, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 1:3, targets.use = 5:8, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

net.down.2 <- net.down[is.finite(net.down$receptor.logFC), ]
computeEnrichmentScore(net.down.2, species = 'mouse')
net.up.2 <- net.up[is.finite(net.up$receptor.logFC), ]
net.up.2 <- net.up.2[is.finite(net.up.2$ligand.logFC), ]
#net.up.2 <- na.omit(net.up.2)
computeEnrichmentScore(net.up.2, species = 'mouse')

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Normal", "Carcinoma")) # set factor level
plotGeneExpression(cellchat, signaling = "SELL", split.by = "datasets", colors.ggplot = T)





# end ---------------------------------------------------------------------

