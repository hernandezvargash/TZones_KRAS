
# KRAS model - Cell-Cell communication for cell types of interest
# stratified labels vs selected immune subtypes

# summary -----------------------------------------------------------------

# Ligand-Receptor Interaction using CellChat: http://www.cellchat.org/
# https://www.nature.com/articles/s41467-021-21246-9


# libraries ---------------------------------------------------------------

rm(list=ls())
gc()

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


# prepare integrated KRAS object ------------------------------------------------------------

sc.all <- readRDS(file = "KRAS_model/objects/integrated_annotated.rds")

head(sc.all[]) # 17721 cells
#DefaultAssay(sc.all) <- "integrated"
Idents(sc.all)
DimPlot(sc.all)


# add labels of interest
# (objects prepared by Chloe)
filenames <- list.files("KRAS_model/objects/objects_for_cellchat", full.names = T)

# load the objects in a list
objects <- sapply(filenames, function(x) mget(load(x)), simplify = TRUE) 
names(objects) <- str_remove(basename(filenames), ".RData")
lapply(objects, ncol)

labels_to_extract <- c("general_label","cluster_v1","labels_Chloe","clusters_v2")
names(labels_to_extract) <- names(objects)

# add the labels to the main object
for(o in 1:length(objects)){
  
  object.temp <- objects[[o]]
  object.temp@meta.data[, labels_to_extract[o]] <- factor(object.temp@meta.data[, labels_to_extract[o]])
  Idents(object.temp) <- object.temp@meta.data[, labels_to_extract[o]]
  
  sc.all <- AddMetaData(sc.all, "NA", col.name = labels_to_extract[o])
  Idents(sc.all) <- sc.all@meta.data[, ncol(sc.all@meta.data)]
  
  for(i in 1:length(levels(object.temp@meta.data[, labels_to_extract[o]]))){
    
    sub1 <- subset(object.temp, idents = levels(object.temp@meta.data[, labels_to_extract[o]])[i])
    int1 <- intersect(colnames(sc.all), colnames(sub1))
    sc.all <- SetIdent(sc.all, cells = int1, value = levels(object.temp@meta.data[, labels_to_extract[o]])[i])
    
    rm(sub1, int1)
  }
  
  sc.all@meta.data[, labels_to_extract[o]] <- Idents(sc.all)
  
  rm(object.temp)
  
}
head(sc.all)
table(sc.all$general_label)
table(sc.all$cluster_v1)
table(sc.all$labels_Chloe)
table(sc.all$clusters_v2)

Idents(sc.all) <- sc.all$major.labels
p1 <- DimPlot(sc.all, label = T, label.box = T, repel = T)
p1 + DimPlot(sc.all, group.by = "general_label")
p1 + DimPlot(sc.all, group.by = "cluster_v1")
p1 + DimPlot(sc.all, group.by = "labels_Chloe")
p1 + DimPlot(sc.all, group.by = "clusters_v2")

table(Idents(sc.all))
table(sc.all$sample)
sc.all$condition <- sc.all$sample
sc.all$condition <- gsub("Hyperplasia.1", "Hyperplasia", sc.all$condition)
sc.all$condition <- gsub("Hyperplasia.2", "Hyperplasia", sc.all$condition)
sc.all$condition <- factor(sc.all$condition, levels = c("Normal", "Hyperplasia", "Dysplasia", "Carcinoma"))
table(sc.all$condition)

# save the object
saveRDS(sc.all, file = "KRAS_model/objects/integrated_annotated.rds")


# prepare cellchat objects ----

sc.all <- readRDS(file = "KRAS_model/objects/integrated_annotated.rds")
head(sc.all)

table(sc.all$labels_Chloe)
table(sc.all$labels_Chloe, sc.all$condition)
Idents(sc.all) <- sc.all$labels_Chloe # from the stratified object

immune.vars <- c("general_label", "cluster_v1", "clusters_v2")
cellchat.list <- list()
for(i in 1:length(immune.vars)){
  sc.temp <- sc.all
  meta <- sc.temp@meta.data
  meta <- meta[meta[, immune.vars[1]] != "NA", ]
  meta[, immune.vars[i]] <- droplevels(meta[, immune.vars[i]])
  for(j in 1:length(levels(meta[, immune.vars[i]]))){
    sub1 <- meta[meta[, immune.vars[i]] == levels(meta[, immune.vars[i]])[j], ]
    int1 <- intersect(colnames(sc.temp), rownames(sub1))
    sc.temp <- SetIdent(sc.temp, cells = int1, value = levels(meta[, immune.vars[i]])[j])
    rm(sub1, int1)
  }
  cellchat.list[[paste0(immune.vars[i])]] <- sc.temp
  rm(sc.temp)  
}
lapply(cellchat.list, function(x) table(Idents(x)))
lapply(cellchat.list, ncol)
cellchat.list <- lapply(cellchat.list, function(x) subset(x, idents = "NA", invert = T))
lapply(cellchat.list, function(x) table(Idents(x), x$condition))

cellchat.dir <- "KRAS_model/cellchat/"

DimPlot(cellchat.list[[1]], label = T, label.box = T, repel = T)
ggsave(paste0(cellchat.dir, "immune_cellchat_labels.jpeg"), width = 12, height = 8)
DimPlot(cellchat.list[[2]], label = T, label.box = T, repel = T)
ggsave(paste0(cellchat.dir, "PMN_cellchat_labels.jpeg"), width = 12, height = 8)
DimPlot(cellchat.list[[3]], label = T, label.box = T, repel = T)
ggsave(paste0(cellchat.dir, "ILC_cellchat_labels.jpeg"), width = 12, height = 8)

# sample 50 cells per cell type and condition within each object
list.cells <- list()
for(l in 1:length(cellchat.list)){
  list.cells.temp.1 <- list()
  obj <- cellchat.list[[l]]
  obj$cellchat.cat <- as.factor(Idents(obj))
  cell.types <- levels(Idents(obj))
  condition <- levels(obj$condition)
  for(i in 1:length(condition)){
    list.cells.temp.2 <- list()
    Idents(obj) <- obj$condition
    temp1 <- WhichCells(obj, idents = condition[i])
    for(j in 1:length(cell.types)){
      Idents(obj) <- obj$cellchat.cat
      temp2 <- WhichCells(obj, idents = cell.types[j])
      temp3 <- intersect(temp1, temp2)
      if(length(temp3) < 50){
        temp3 <- sample(temp3, size=length(temp3))
      } else {
        temp3 <- sample(temp3, size=50)
      }
      list.cells.temp.2[[j]] <- temp3
      rm(temp2, temp3)
    }
    list.cells.temp.1[[i]] <- unlist(list.cells.temp.2)
    rm(temp1)
  }
  list.cells[[l]] <- unlist(list.cells.temp.1)
  obj <- obj[, unlist(list.cells[[l]])]
  cellchat.list[[l]] <- obj
  rm(list.cells.temp, obj, cell.types)
}

lapply(cellchat.list, function(x) table(Idents(x)))
lapply(cellchat.list, ncol)
lapply(cellchat.list, function(x) table(Idents(x), x$condition))

# save the object
saveRDS(cellchat.list, file = "KRAS_model/objects/cellchat_list.rds")



# CellChat integrated ------------------------------------------

cellchat.list <- readRDS(file = "KRAS_model/objects/cellchat_list.rds")
names(cellchat.list) <- c("Immune","PMNs","ILCs")

for(i in 1:length(cellchat.list)){
  sc.temp <- cellchat.list[[i]]
  DefaultAssay(sc.temp) <- "RNA"
  output.folder <- paste0("KRAS_model/cellchat/All/", names(cellchat.list)[i], "/")
  sc.temp$cellchat.cat <- droplevels(sc.temp$cellchat.cat)
  Idents(sc.temp) <- sc.temp$cellchat.cat
  cellchat <- createCellChat(object = sc.temp, group.by = "cellchat.cat", assay = "RNA")
  cellchat@DB <- CellChatDB.mouse
  cellchat <- subsetData(cellchat)
  
  future::plan("multisession", workers = 12)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  options(future.globals.maxSize = 8000 * 1024^2)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  jpeg(paste0(output.folder, "netVisual_circle_all.cells.jpeg"), height = 900, width = 900, quality = 100)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  mat <- cellchat@net$weight
  for (j in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[j, ] <- mat[j, ]
    jpeg(paste0(output.folder, "netVisual_circle_", make.names(rownames(mat))[j], ".jpeg"), height = 900, width = 900, quality = 100)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j])
    dev.off()
  }
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  pdf(paste0(output.folder, "netAnalysis_heatmap.pdf"), height = 14, width = 10)
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 8, height = 25) +
    netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 8, height = 25)
  dev.off()
  df.net <- subsetCommunication(cellchat)
  head(df.net)
  df.net.2 <- subsetCommunication(cellchat, slot.name = "netP") # signaling pathways
  head(df.net.2)
  max.y <- max(c(table(df.net$target), table(df.net$source)))
  jpeg(paste0(output.folder, "interactions_barplot.jpeg"), height = 900, width = 500, quality = 100)
  par(mfrow = (c(2,1)), mar=c(12,5,2,2))
  barplot(table(df.net$source), col = scPalette(length(levels(cellchat@idents))), 
          las = 2, ylim = c(0, max.y),
          main = "Outgoing interactions per cell type")
  barplot(table(df.net$target), col = scPalette(length(levels(cellchat@idents))), las = 2,
          las = 2, ylim = c(0, max.y),
          main = "Incoming interactions per cell type")
  dev.off()
  netAnalysis_signalingRole_scatter(cellchat)
  ggsave(paste0(output.folder, "netAnalysis_signalingRole_scatter.jpeg"))
  
  
  write.csv(df.net, file = paste0(output.folder, "cellchat.net.1.csv"))
  write.csv(df.net.2, file = paste0(output.folder, "cellchat.net.2.csv"))
  
  save(cellchat, file = paste0(output.folder, "cellchat", names(cellchat.list)[i], ".RData"))
  
}





# CellChat by condition ------------------------------------------

cellchat.list <- readRDS(file = "KRAS_model/objects/cellchat_list.rds")
names(cellchat.list) <- c("Immune","PMNs","ILCs")

for(i in 1:length(cellchat.list)){
  sc.temp.2 <- cellchat.list[[i]]
  DefaultAssay(sc.temp.2) <- "RNA"
  for(k in 1:length(levels(sc.temp.2$condition))){
    output.folder <- paste0("KRAS_model/cellchat/", levels(sc.temp.2$condition)[k], "/", names(cellchat.list)[i], "/")
    sc.temp <- subset(sc.temp.2, condition == levels(sc.temp.2$condition)[k])
    sc.temp$cellchat.cat <- droplevels(sc.temp$cellchat.cat)
    Idents(sc.temp) <- sc.temp$cellchat.cat
    cellchat <- createCellChat(object = sc.temp, group.by = "cellchat.cat", assay = "RNA")
    cellchat@DB <- CellChatDB.mouse
    cellchat <- subsetData(cellchat)
    
    future::plan("multisession", workers = 12)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    options(future.globals.maxSize = 8000 * 1024^2)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    if (dim(cellchat@data)[2] > 50){
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      groupSize <- as.numeric(table(cellchat@idents))
      jpeg(paste0(output.folder, "netVisual_circle_all.cells.jpeg"), height = 900, width = 900, quality = 100)
      netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
      dev.off()
      mat <- cellchat@net$weight
      for (j in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[j, ] <- mat[j, ]
        jpeg(paste0(output.folder, "netVisual_circle_", make.names(rownames(mat))[j], ".jpeg"), height = 900, width = 900, quality = 100)
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j])
        dev.off()
        rm(mat2)
      }
      cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  pdf(paste0(output.folder, "netAnalysis_heatmap.pdf"), height = 14, width = 10)
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 8, height = 25) +
    netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 8, height = 25)
  dev.off()
      df.net <- subsetCommunication(cellchat)
      head(df.net)
      df.net.2 <- subsetCommunication(cellchat, slot.name = "netP") # signaling pathways
      head(df.net.2)
      max.y <- max(c(table(df.net$target), table(df.net$source)))
      jpeg(paste0(output.folder, "interactions_barplot.jpeg"), height = 900, width = 500, quality = 100)
      par(mfrow = (c(2,1)), mar=c(12,5,2,2))
      barplot(table(df.net$source), col = scPalette(length(levels(cellchat@idents))), 
              las = 2, ylim = c(0, max.y),
              main = "Outgoing interactions per cell type")
      barplot(table(df.net$target), col = scPalette(length(levels(cellchat@idents))), las = 2,
              las = 2, ylim = c(0, max.y),
              main = "Incoming interactions per cell type")
      dev.off()
      netAnalysis_signalingRole_scatter(cellchat)
      ggsave(paste0(output.folder, "netAnalysis_signalingRole_scatter.jpeg"))
      
      
      write.csv(df.net, file = paste0(output.folder, "cellchat.net.1.csv"))
      write.csv(df.net.2, file = paste0(output.folder, "cellchat.net.2.csv"))
      
      save(cellchat, file = paste0(output.folder, "cellchat", names(cellchat.list)[i], ".RData"))
      
    }
    
    rm(cellchat, df.net, df.net.2, sc.temp)
  }
  rm(sc.temp.2)
}



# L-R pairs of interest ------------------------------------------------

#CXCL5-CXCR2 Basal AC proximal ----> immune atlas
#CXCL5-CXCR2 supra Basal AC proximal ----> immune atlas
#CXCL7-CXCR2 Basal AC proximal ----> immune atlas
#CXCL7-CXCR2 supra Basal AC proximal ----> immune atlas

#CXCL5-CXCR2 Basal AC proximal ----> PMNs
#CXCL5-CXCR2 supra Basal AC proximal ----> PMNs
#CXCL7-CXCR2 Basal AC proximal ----> PMNs
#CXCL7-CXCR2 supra Basal AC proximal ----> PMNs

output.dir <- "KRAS_model/cellchat/CXCL_pathway/"

dataset <- "All"
cell.types <- c("Immune", "PMNs")
cellchat.of.interest <- paste0(dataset, ".", cell.types)

i = 1
i = 2

cellchat.files <- list.files("KRAS_model/cellchat/", full.names = T, recursive = T, pattern = ".RData")
load(cellchat.files[grep(paste0(dataset,"/",cell.types[i]), cellchat.files)])
levels(cellchat@idents)

pathways.show <- c("CXCL") 
# Immune
sources.use <- c(13,15)
targets.use <- 1:6
sources.use <- 7:16
targets.use <- 2
# PMNs
sources.use <- c(12,14)
targets.use <- 1:5

# contribution or each L-R pair within selected pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
ggsave(paste0(output.dir, "netAnalysis_contribution_", cellchat.of.interest[i], ".jpeg"), height = 6, width = 10)
par(mfrow=c(1,1), mar=c(0,5,0,0))
# cell communication with all cell types
jpeg(paste0(output.dir, "netVisual_heatmap_", cellchat.of.interest[i], ".jpeg"), height = 400, width = 600, quality = 100)
netVisual_heatmap(cellchat, signaling = pathways.show, 
#                  sources.use = sources.use, targets.use = targets.use,
                  color.heatmap = "Reds")
dev.off()
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", title.space = 1) # signaling.name = "aggregated CXCL L-R pairs"
ggsave(paste0(output.dir, "netVisual_aggregate_", cellchat.of.interest[i], ".jpeg"), height = 6, width = 10)
# chord diagram with participating cell types
jpeg(paste0(output.dir, "netVisual_chord_", cellchat.of.interest[i], ".jpeg"), height = 600, width = 600, quality = 100)
netVisual_chord_gene(cellchat, sources.use = c(13,15), targets.use = 1:6, 
                     slot.name = "net",
                     directional = 1,
                     title.name = "CXCL L-R chord diagram",
                     legend.pos.x = -10,
                     legend.pos.y = 30,
                     lab.cex = 1.2,
                     signaling = pathways.show)
dev.off()
# bubble heatmap of all significant interactions
jpeg(paste0(output.dir, "netVisual_bubble_outgoing_", cellchat.of.interest[i], ".jpeg"), height = 400, width = 600, quality = 100)
netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, remove.isolate = FALSE, 
                 font.size = 12, font.size.title = 14, angle.x = 45, dot.size.min = 3, dot.size.max = 6,
                 signaling = pathways.show)
dev.off()
netVisual_bubble(cellchat, sources.use = targets.use, targets.use = sources.use, remove.isolate = FALSE, 
                 font.size = 12, font.size.title = 14, angle.x = 45, dot.size.min = 3, dot.size.max = 6,
                 signaling = pathways.show)
# gene expression of all cell types
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")
ggsave(paste0(output.dir, "plotGeneExpression_", cellchat.of.interest[i], ".jpeg"), height = 6, width = 10)
# Heatmap showing the centrality scores/importance of cell groups as senders, receivers, mediators and influencers 
# in a single intercellular communication network
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 18, height = 5, 
                                  font.size.title = 16, font.size = 14)
ggsave(paste0(output.dir, "netAnalysis_signalingRole_network_", cellchat.of.interest[i], ".jpeg"), height = 6, width = 10)


pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
# CXCL5_CXCR2
LR.show <- pairLR.CXCL[2,] # 1 for immune, 2 for PMNs
jpeg(paste0(output.dir, "netVisual_", LR.show, "_circle_", cellchat.of.interest[i], ".jpeg"), height = 600, width = 600, quality = 100)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = NULL, layout = "circle")
dev.off()
jpeg(paste0(output.dir, "netVisual_", LR.show, "_chord_", cellchat.of.interest[i], ".jpeg"), height = 700, width = 700, quality = 100)
netVisual_individual(cellchat, signaling = pathways.show, height = 5, pairLR.use = LR.show, vertex.receiver = NULL, layout = "chord")
dev.off()
# PPBP_CXCR2 (CXCL7)
LR.show <- pairLR.CXCL[5,] # 4 for Immune, 5 for PMNs
jpeg(paste0(output.dir, "netVisual_", LR.show, "_circle_", cellchat.of.interest[i], ".jpeg"), height = 600, width = 600, quality = 100)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = NULL, layout = "circle")
dev.off()
jpeg(paste0(output.dir, "netVisual_", LR.show, "_chord_", cellchat.of.interest[i], ".jpeg"), height = 700, width = 700, quality = 100)
netVisual_individual(cellchat, signaling = pathways.show, height = 5, pairLR.use = LR.show, vertex.receiver = NULL, layout = "chord")
dev.off()


# plots for figures ------------------------------------------------

output.dir <- "KRAS_model/cellchat/figures/"

dataset <- c("Normal", "Hyperplasia", "Dysplasia", "Carcinoma")
cell.types <- c("PMNs")
cellchat.of.interest <- paste0(dataset, ".", cell.types)

cellchat.files <- list.files("KRAS_model/cellchat/", full.names = T, recursive = T, pattern = ".RData")

for(d in 2:length(dataset)){
  load(cellchat.files[grep(paste0(dataset[d],"/",cell.types), cellchat.files)])
  groupSize <- as.numeric(table(cellchat@idents))
  mat <- cellchat@net$weight
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#  mat2["AC suprabasal proximal", ] <- mat["AC suprabasal proximal", ]
  mat2["AC basal proximal", ] <- mat["AC basal proximal", ]
  #netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F)
  rownames(mat2) <- rep("  ", nrow(mat2))
  colnames(mat2) <- rep("  ", nrow(mat2))
  par(mfrow=c(1,1), mar=c(0,0,0,0))
  netVisual_circle(mat2, #cellchat@net$weight,
                   remove.isolate = F,
                   #                 sources.use = "AC suprabasal proximal",
                   title.name = NULL, 
                   shape = "circle",
                   margin = c(0.2,0,0.4,0), arrow.size = 0.5,
                   #                 text.x = 0, text.y = 1.5,
                   vertex.label.cex = 1.2, alpha.edge = 0.8,
                   vertex.weight = groupSize, weight.scale = T)
  text(0,1.4, dataset[d] , cex = 1.8)
}



load(cellchat.files[grep(paste0(dataset[2],"/",cell.types), cellchat.files)])
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F)
mat <- cellchat@net$weight
rownames(mat)
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2["AC suprabasal proximal", ] <- mat["AC suprabasal proximal", ]
#netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F)
rownames(mat2) <- rep("  ", nrow(mat2))
colnames(mat2) <- rep("  ", nrow(mat2))
par(mfrow=c(1,1), mar=c(0,0,0,0))
netVisual_circle(mat2, #cellchat@net$weight,
                 remove.isolate = F,
#                 sources.use = "AC suprabasal proximal",
                 title.name = NULL, 
                 shape = "circle",
                 margin = c(0.2,0,0.4,0), arrow.size = 0.5,
#                 text.x = 0, text.y = 1.5,
                 vertex.label.cex = 1.2, alpha.edge = 0.8,
                 vertex.weight = groupSize, weight.scale = T)
text(0,1.5,"Hyperplasia", cex = 1.8)





dataset <- "All"
cell.types <- c("Immune")
cellchat.of.interest <- paste0(dataset, ".", cell.types)
cellchat.files <- list.files("KRAS_model/cellchat/", full.names = T, recursive = T, pattern = ".RData")
load(cellchat.files[grep(paste0(dataset,"/",cell.types), cellchat.files)])
levels(cellchat@idents)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, 
                                         signaling = "CXCL",
                                         pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, 
                                         signaling = "CXCL",
                                         pattern = "incoming")
ht1 + ht2

netAnalysis_signalingRole_heatmap(cellchat, 
                                  signaling = "CXCL",
                                  pattern = "all")

# end ---------------------------------------------------------------------

