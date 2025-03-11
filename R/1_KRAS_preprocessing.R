
# KRAS model - Preprocessing

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

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/TZones/")

suppressPackageStartupMessages({
  library(Seurat);  library(sctransform);  library(SeuratDisk); library(SeuratWrappers); library(SingleCellExperiment)
  library(dplyr);  library(stringr) # library(xlsx); 
  library(ggplot2); library(cowplot); library(viridis); library(gridExtra); library(RColorBrewer); library(patchwork); library(hrbrthemes); library(pals)
  library(SingleR); library(scRNAseq); library(scibet)
  library(slingshot, quietly = T)
  library(mclust, quietly = T); library(celldex); library(KernSmooth); library(dendextend); library(circlize)
  library(enrichR); library(clusterProfiler); library(enrichplot); library(DOSE); library(pathview); library(ReactomePA)
  library(ggpubr); library(eulerr); library(alluvial)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db); library(biomaRt)
  library(Scillus)
  library(Tempora)
  library(infercnv)
  library(SCENIC)
  library(snow)
  library(DoubletFinder)
})

set.seed(1)


# loading from loom files  -----------------------------------------------------------------

loom.files <- list.files(path = "shared/KRAS/data/loom", recursive = T, full.names = T, pattern = ".loom")

# two hyperplasia samples will be merged later
loom.files <- loom.files[c(4,3,2,1,5)]
#names <- str_sub(loom.files, 54, -6)

names <- c("Normal","Hyperplasia.1","Hyperplasia.2","Dysplasia","Carcinoma")

list.samples <- list()

for(i in 1:length(loom.files)){
  list.samples[[i]] <- ReadVelocity(file = loom.files[i])
  list.samples[[i]] <- as.Seurat(list.samples[[i]])
  list.samples[[i]][["RNA"]] <- list.samples[[i]][["spliced"]]
  DefaultAssay(list.samples[[i]]) <- "RNA"
  list.samples[[i]]$sample <- names[[i]]
}

names(list.samples) <- names

saveRDS(list.samples, file = "KRAS_model/objects/list.samples.rds")



# loading from barcode matrices (don't run, only for testing) -----------------------------------------------------------------

dirs <- list.dirs(path = "shared/KRAS/data/matrices", recursive = F, full.names = T)
dirs <- dirs[c(6,4,3,1,5)]

#names <- str_sub(dirs, 58, -1)
names <- c("Normal","Hyperplasia.1","Hyperplasia.2","Dysplasia","Carcinoma")

list.samples <- list()

for(i in 1:length(dirs)){
  list.samples[[i]] <- Read10X(data.dir = dirs[i])
  list.samples[[i]] <- CreateSeuratObject(list.samples[[i]])
  list.samples[[i]]$sample <- names[[i]]
}

names(list.samples) <- names



# initial filtering -------------------------------------------------------------

lapply(list.samples, dim)

object.list <- lapply(X = list.samples, FUN = PercentageFeatureSet, pattern = "^mt-", col.name = "percent.mt")
object.list <- lapply(X = object.list, FUN = PercentageFeatureSet, pattern = "^Rp[sl][[:digit:]]", col.name = "percent.ribo")
object.list <- lapply(X = object.list, FUN = subset, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)

lapply(object.list, dim)

head(object.list[[5]])



# doublet estimation -------------------------------------------------------


for(i in 1:length(object.list)){
  
  temp1 <- SCTransform(object.list[[i]])
  temp1 <- RunPCA(temp1)
  temp1 <- RunUMAP(temp1, dims = 1:10)
  
  ## pK Identification
  sweep.res <- paramSweep_v3(temp1, PCs = 1:20, sct = T) # using SCTransform normalization
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line() + ggtitle(paste0(names(object.list)[i]))
  
  pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric))
  pK <- as.numeric(as.character(pK$pK))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- temp1@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(temp1@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # run doubletFinder 
  temp1 <- doubletFinder_v3(temp1, 
                            PCs = 1:20, 
                            pN = 0.25, 
                            pK = pK, 
                            nExp = nExp_poi.adj,
                            reuse.pANN = FALSE, sct = T) # SCT was used for normalization
  names(temp1@meta.data)[ncol(temp1@meta.data)] <- "DF.classifications"
  DimPlot(temp1, reduction = 'umap', group.by = "DF.classifications") +
    ggtitle(paste0(names(object.list)[i]))
  table(temp1$DF.classifications)
  
  object.list[[i]] <- temp1
  
  rm(temp1, nExp_poi.adj, pK)
  
}


lapply(object.list, function(x) table(x@meta.data$DF.classifications))
doublet.plot <- function(x){
  DimPlot(x, group.by = colnames(x@meta.data)[ncol(x@meta.data)]) + ggtitle(x@meta.data$sample)
}
plist <- lapply(object.list, doublet.plot)
do.call(grid.arrange, plist)

saveRDS(object.list, file = "KRAS_model/objects/list.samples.prefiltered.doublets.rds")




# integration -------------------------------------------------------------


object.list <- readRDS("KRAS_model/objects/list.samples.prefiltered.doublets.rds")

lapply(object.list, dim)
object.list <- lapply(X = object.list, FUN = subset, subset = DF.classifications == "Singlet")

object.list <- lapply(X = object.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = features)
object.list <- lapply(X = object.list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = object.list, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
int.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
int.sct$sample <- factor(int.sct$sample, levels = c("Normal","Hyperplasia.1","Hyperplasia.2","Dysplasia","Carcinoma"))
table(int.sct$sample)
# Normal Hyperplasia.1 Hyperplasia.2     Dysplasia     Carcinoma 
# 2972          1150          1865          5020          6714 

int.sct <- RunPCA(int.sct, verbose = FALSE)
int.sct <- FindNeighbors(int.sct, dims = 1:10)
int.sct <- RunUMAP(int.sct, reduction = "pca", dims = 1:30)
int.sct <- FindClusters(int.sct, resolution = 0.5)


VlnPlot(int.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)
VlnPlot(int.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), group.by = "sample", ncol = 4)

DimPlot(int.sct, label = T)

p1 <- DimPlot(int.sct)
p2 <- DimPlot(int.sct, group.by = "sample")
p1+p2

DimPlot(int.sct, split.by = "sample")

head(int.sct[[]])
table(int.sct$orig.ident)
table(int.sct$sample)



saveRDS(int.sct, file = "KRAS_model/objects/integrated.rds")





# end ---------------------------------------------------------------------


