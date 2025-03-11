
# KRAS model - CNV analyses

# summary -----------------------------------------------------------------

# https://github.com/broadinstitute/infercnv
# inferCNV used as recommended one sample at a time 
# (after split of the integrated and annotated object)
# and using stroma and immune cells as "normal" reference


# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/TZones/")

suppressPackageStartupMessages({
  library(Seurat);  library(sctransform);  library(SeuratDisk); library(SeuratWrappers); library(SingleCellExperiment)
  library(dplyr); library(xlsx);  library(stringr)
  library(ggplot2); library(cowplot); library(viridis); library(gridExtra); library(RColorBrewer); library(patchwork); library(hrbrthemes); library(pals)
  library(SingleR); library(scRNAseq); library(scibet)
  library(slingshot, quietly = T)
  library(mclust, quietly = T); library(celldex); library(KernSmooth); library(dendextend); library(circlize); library(Nebulosa)
  library(enrichR); library(clusterProfiler); library(enrichplot); library(DOSE); library(pathview); library(ReactomePA)
  library(ggpubr); library(eulerr); library(alluvial)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db); library(biomaRt)
  library(Scillus); library(Nebulosa)
  library(Tempora)
  library(infercnv)
  library(SCENIC)
  library(snow)
  
})

#options(Seurat.object.assay.version = "v3")

set.seed(5)


# CNV analysis ------------------------------------------------------------

sc <- readRDS(file = "KRAS_model/objects/integrated_annotated.rds")
DimPlot(sc, split.by = "sample")
table(sc$sample)
table(sc$major.labels, sc$sample)
# only 1 Anal gland cell in Hyperplasia.2, needs to be removed
DefaultAssay(sc) <- "RNA"
sc.split <- SplitObject(sc, split.by = "sample")
names(sc.split)
head(sc.split$Normal)[]
sc.sub <- sc.split$Hyperplasia.2
sc.sub <- subset(sc, idents = "Hyperplasia.2")
sc.sub <- subset(sc.sub, idents = "Anal gland", invert = T)
sc.split$Hyperplasia.2 <- sc.sub

for(i in 1:length(sc.split)){
  
  sc.temp <- sc.split[[i]]
  counts.temp <- as.matrix(sc.temp[["RNA"]]$counts)
  annots.temp <- sc.temp[[]]
  annots.temp <- annots.temp[, c("major.labels"), drop = F]

  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  order.temp <- getBM(attributes=c('mgi_symbol', 'chromosome_name', 'start_position', 'end_position', 'strand'),
                      filters=c('mgi_symbol'), values=list(rownames(counts.temp)), mart=mouse)
  
  rownames(order.temp) <- make.names(order.temp$mgi_symbol, unique = T)
  order.temp <- order.temp[, -1]

  infercnv_object <- infercnv::CreateInfercnvObject(raw_counts_matrix = counts.temp, 
                                                    gene_order_file = order.temp,
                                                    annotations_file = annots.temp,
                                                    ref_group_names=c("Myeloid", "T cells", "B cells", "Stromal"))
  options(scipen = 100)
  out_dir = paste0("~/Dropbox/BioInfo/Lab/TZones/KRAS_model/results/CNV/", names(sc.split)[i])
  infercnv_obj_default = infercnv::run(
    infercnv_object,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=out_dir,
    cluster_by_groups=TRUE, 
    analysis_mode = "samples", # using subclusters gives a "ScaleData" error
    tumor_subcluster_partition_method = "leiden",
    num_threads = 12,
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=T,
    no_prelim_plot=F,
    png_res=300
  )
  
  # add info to seurat object
  sc.temp = infercnv::add_to_seurat(infercnv_output_path = out_dir,
                                     seurat_obj = sc.temp, # optional
                                     top_n=10)
  
  png(paste0(out_dir, "/top_loss_1.png"))
  FeaturePlot(sc.temp, reduction="umap", features="top_loss_1") + ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))
  dev.off()
  png(paste0(out_dir, "/top_dupli_1.png"))
  FeaturePlot(sc.temp, reduction="umap", features="top_dupli_1") + ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))
  dev.off()
  
  sc.split[[i]] <- sc.temp
  
  rm(sc.temp, counts.temp, order.temp, annots.temp, out_dir)
  
}

f1 <- FeaturePlot(sc.split$Normal, reduction="umap", features="top_loss_1") + ggtitle("Normal: top loss 1")
f2 <- FeaturePlot(sc.split$Hyperplasia.1, reduction="umap", features="top_loss_1") + ggtitle("Hyperplasia.1: top loss 1")
f3 <- FeaturePlot(sc.split$Hyperplasia.2, reduction="umap", features="top_loss_1") + ggtitle("Hyperplasia.2: top loss 1")
f4 <- FeaturePlot(sc.split$Dysplasia, reduction="umap", features="top_loss_1") + ggtitle("Dysplasia: top loss 1")
f5 <- FeaturePlot(sc.split$Carcinoma, reduction="umap", features="top_loss_1") + ggtitle("Carcinoma: top loss 1")
f1+f2+f3+f4+f5
f6 <- FeaturePlot(sc.split$Normal, reduction="umap", features="top_dupli_1") + ggtitle("Normal: top dupli 1")
f7 <- FeaturePlot(sc.split$Hyperplasia.1, reduction="umap", features="top_dupli_1") + ggtitle("Hyperplasia.1: top dupli 1")
f8 <- FeaturePlot(sc.split$Hyperplasia.2, reduction="umap", features="top_dupli_1") + ggtitle("Hyperplasia.2: top dupli 1")
f9 <- FeaturePlot(sc.split$Dysplasia, reduction="umap", features="top_dupli_1") + ggtitle("Dysplasia: top dupli 1")
f10 <- FeaturePlot(sc.split$Carcinoma, reduction="umap", features="top_dupli_1") + ggtitle("Carcinoma: top dupli 1")
grid.arrange(f1,f6,f2,f7,f3,f8,f4,f9,f5,f10, ncol = 2)


saveRDS(sc.split, file = "KRAS_model/objects/KRAS_CNV_list.rds")

sc.merge <- merge(sc.split[[1]], c(sc.split[[2]], sc.split[[3]], sc.split[[4]], sc.split[[5]]))
head(sc.merge)
table(sc.merge$sample)
meta <- sc.merge@meta.data
head(meta)

write.csv(meta, file = "KRAS_model/objects/CNV_metadata.csv")



# annotate stratified with CNV data ----

rm(list=ls())

cnv <- read.csv("KRAS_model/objects/CNV_metadata.csv", row.names = 1)
head(cnv)

sc <- readRDS(file = "KRAS_model/objects/stratified.rds")
setdiff(rownames(sc@meta.data), rownames(cnv))
setdiff(rownames(cnv), rownames(sc@meta.data))
int1 <- intersect(rownames(sc@meta.data), rownames(cnv))
cnv <- cnv[int1, ]
all(rownames(sc@meta.data) == rownames(cnv))
colnames(cnv)
sc <- AddMetaData(sc, cnv[, 30:239])
head(sc)

DimPlot(sc, label = T)
FeaturePlot(sc, reduction="umap", features="top_loss_1", split.by = "sample")
p1 <- plot_density(sc, "top_loss_1")
table(sc$top_loss_1)
DimPlot(sc, group.by="top_loss_2", pt.size=0.5 )
p2 <- plot_density(sc, "top_loss_2")

p3 <- plot_density(sc, "top_dupli_2")
p4 <- plot_density(sc, "top_dupli_7")
p5 <- plot_density(sc, "top_dupli_9")
p6 <- plot_density(sc, "top_dupli_10")

grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 2)

plot_density(sc, "proportion_cnv_9")
FeaturePlot(sc, "proportion_cnv_9", split.by = "sample")

FeaturePlot(sc, "proportion_scaled_dupli_14", split.by = "sample")



# test in integrated object ----

sc <- readRDS(file = "KRAS_model/objects/stratified_Atlas.rds")
DimPlot(sc)
DimPlot(sc, split.by = "condition")
table(sc$condition)
table(sc$cell_type, sc$condition)
DefaultAssay(sc) <- "RNA"

counts.temp <- as.matrix(sc[["RNA"]]$counts)
annots.temp <- sc[[]]
annots.temp <- annots.temp[, c("condition","cell_type"), drop = F]
head(annots.temp)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
order.temp <- getBM(attributes=c('mgi_symbol', 'chromosome_name', 'start_position', 'end_position', 'strand'),
                    filters=c('mgi_symbol'), values=list(rownames(counts.temp)), mart=mouse)

rownames(order.temp) <- make.names(order.temp$mgi_symbol, unique = T)
order.temp <- order.temp[, -1]
head(order.temp)

infercnv_object <- infercnv::CreateInfercnvObject(raw_counts_matrix = counts.temp, 
                                                  gene_order_file = order.temp,
                                                  annotations_file = annots.temp,
                                                  ref_group_names=c("Atlas"))
#options(scipen = 100)
out_dir = paste0("~/Dropbox/BioInfo/Lab/TZones/KRAS_model/results/CNV/Atlas")
infercnv_obj_default = infercnv::run(
  infercnv_object,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=out_dir,
  cluster_by_groups=TRUE, 
  analysis_mode = "samples", # using subclusters gives a "ScaleData" error
  tumor_subcluster_partition_method = "leiden",
  num_threads = 12,
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=T,
  no_prelim_plot=F,
  png_res=300
)

# add info to seurat object
sc = infercnv::add_to_seurat(infercnv_output_path = out_dir,
                                  seurat_obj = sc, # optional
                                  top_n=10)
head(sc)

DimPlot(sc, group.by = "cell_type")
FeaturePlot(sc, reduction="umap", features="top_loss_1") + ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))
FeaturePlot(sc, reduction="umap", features="top_loss_1") + ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))
FeaturePlot(sc, reduction="umap", features="top_dupli_1") + ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))
FeaturePlot(sc, reduction="umap", features="has_loss_13") + ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))

FeaturePlot(sc, "top_dupli_1", split.by = "condition")
DimPlot(sc, group.by = "top_dupli_1", split.by = "condition")

colnames(sc@meta.data)
gg.list <- list()
for(i in 253:271){
  gg.list[[i-252]] <- plot_density(sc, colnames(sc@meta.data)[i])
}
gg.list[[20]] <- DimPlot(sc, group.by = "cell_type", label = T, label.size = 2) + NoLegend()
do.call("grid.arrange", c(gg.list, ncol = 5)) 

gg.list <- list()
for(i in 253:271){
  gg.list[[i-252]] <- DimPlot(sc, group.by = colnames(sc@meta.data)[i], split.by = "condition")
}
do.call("grid.arrange", c(gg.list[1:4], ncol = 1)) 
do.call("grid.arrange", c(gg.list[5:9], ncol = 1)) 
do.call("grid.arrange", c(gg.list[10:14], ncol = 1)) 
do.call("grid.arrange", c(gg.list[15:19], ncol = 1)) 



saveRDS(sc, file = "KRAS_model/objects/Atlas_stratified_KRAS_CNV.rds")



# end ---------------------------------------------------------------------

sessionInfo()
