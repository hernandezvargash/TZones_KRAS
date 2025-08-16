# Revisions received on 22/05/25
# Deadline: end of June
# New revisions: 11/08/25
# Deadline: 25 September


# libraries ---------------------------------------------------------------

rm(list=ls())
gc()

setwd("~/Dropbox/BioInfo/Lab/TZones/")

suppressPackageStartupMessages({
  library(Seurat);  library(sctransform);  library(SeuratDisk); #library(SeuratWrappers); library(SingleCellExperiment)
  library(dplyr);  library(stringr) #library(xlsx); 
  library(ggplot2); library(cowplot); library(viridis); library(gridExtra); library(RColorBrewer); library(patchwork); library(hrbrthemes); library(pals)
  library(SingleR); library(scRNAseq); library(scibet)
  library(slingshot, quietly = T)
  library(mclust, quietly = T); library(celldex); library(KernSmooth); library(dendextend); library(circlize)
  library(enrichR); library(clusterProfiler); library(enrichplot); library(DOSE); library(pathview); library(ReactomePA)
  library(ggpubr); library(eulerr); library(alluvial); library(scales)
  library(escape)
  library(Nebulosa)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db); library(biomaRt)
  library(snow)
  library(speckle) # propeller
})

set.seed(1)

#results.dir <- "manuscript/KRAS/revisions_June_2025/"
results.dir <- "manuscript/KRAS/revisions_Aug_2025/"
data.dir <- "KRAS_model/objects/"


# all cells ----

int.sct <- readRDS(file = paste0(data.dir, "integrated_annotated.rds"))
head(int.sct[[]])
DimPlot(int.sct, label = T, label.box = T) + NoLegend()

# identify missing cell type
#colonocyte-like are positive for EpCam, Keratin 8, Keratin 18, Keratin 20a, Cxcr2, cd24a, Fcgr3, S100a8/9, Muc3, Car4, Saa1
FeaturePlot(int.sct, features = c("Epcam", "Krt8", "Krt18", "Krt20a", "Cxcr2", "Cd24a", "Fcgr3", "S100a8", "S100a9", "Muc3", "Car4", "Saa1"), ncol = 4)
# these cells were identified in the immune subset, as seurat cluster 5
# and having mixed annotations Neutrophils/Cancer cells
# to.remove <- WhichCells(immune, expression = new.cell.labels == "Neutrophils" & customclassif == "Cancer cells") # 262 cells
DimPlot(int.sct, label = T, label.box = T, cells.highlight = to.remove) + NoLegend()
int.sct <- SetIdent(int.sct, cells = to.remove, value = "Colonocyte-like")
DimPlot(int.sct, label = T, label.box = T, repel = T) + NoLegend()


# rename identities

Idents(int.sct)
int.sct <- RenameIdents(int.sct, 
                            `Myeloid` = "Immune cells",
                            `Stratified` = "Stratified epithelium",
                            `Rectum` = "Glandular epithelium",
                            `T cells` = "Immune cells",
                            `Anal gland` = "Anal gland",
                            `Hair follicle` = "Hair follicle",
                            `B cells` = "Immune cells",
                            `Stromal` = "Mesenchymal cells"
                           )
int.sct$cell_type <- Idents(int.sct)

DimPlot(int.sct, label = F, label.box = T, repel = T) +
  ggtitle("Total Cells")
ggsave(filename = paste0(results.dir, "all.cells.identities.png"), width = 7, height = 5)

# replace condition categories in case of alphabetical ordering of the conditions
Idents(int.sct) <- int.sct$condition
int.sct <- RenameIdents(int.sct, 
                           `Normal` = "A.Normal",
                           `Hyperplasia` = "B.Hyperplasia",
                           `Dysplasia` = "C.Dysplasia",
                           `Carcinoma` = "D.Carcinoma"
)

# convert main variables to vectors to avoid numeric conversion
int.sct$orig.ident <- as.character(int.sct$orig.ident)
int.sct$sample <- as.character(int.sct$sample)
int.sct$condition <- as.character(Idents(int.sct))
int.sct$major.labels <- as.character(int.sct$major.labels)
int.sct$cell_type <- as.character(int.sct$cell_type)


SaveH5Seurat(int.sct, filename = paste0(data.dir, "all.cells.h5Seurat"), overwrite = TRUE)
Convert(paste0(data.dir, "all.cells.h5Seurat"), dest = "h5ad", overwrite = TRUE)


saveRDS(int.sct, file = paste0(data.dir, "all.cells_manuscript.rds"))


# Aug 2025
# DimPlot with all cells, split by sample

int.sct <- readRDS(file = paste0(data.dir, "all.cells_manuscript.rds"))
head(int.sct[[]])

levels(int.sct$cell_type)
Idents(int.sct) <- factor(int.sct$cell_type,
                          levels = c("Mesenchymal cells", "Stratified epithelium", "Immune cells", 
                                     "Glandular epithelium", "Anal gland", "Hair follicle", "Colonocyte-like"))

DimPlot(int.sct, label = F, label.box = T, repel = T) +
  ggtitle("Total Cells")
ggsave(filename = paste0(results.dir, "non.cherry.picked.umap.2A.png"), width = 7, height = 5)

DimPlot(int.sct, 
        split.by = "condition",
        label = F, label.box = T, repel = T) +
  ggtitle("Total Cells")
ggsave(filename = paste0(results.dir, "non.cherry.picked.umap.2A.split.png"), width = 18, height = 5)

table(int.sct$condition)
# A.Normal B.Hyperplasia   C.Dysplasia   D.Carcinoma 
# 2972          3015          5020          6714 


## proportions ----

table(int.sct$cell_type)
Idents(int.sct) <- int.sct$cell_type

type.colors = alphabet(length(levels(as.factor(int.sct$cell_type))))
show_col(type.colors)
names(type.colors) <- levels(as.factor(int.sct$cell_type))

# Plot cell type proportions
tab1 <- table(int.sct$condition, int.sct$cell_type)
tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Cell_Type", "Condition","Proportion")
head(tab1)
ggplot(tab1, aes(x = Condition, y = Proportion, fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = type.colors) +
  ggtitle("Cell Type Proportions") +
  theme_ipsum() +
  xlab("") #+ RotatedAxis()
ggsave(paste0(results.dir, "celltype.proportions_all.cells.png"), width=7, height=5, dpi=300)

# when having replicates (only hyperplasia), we can use propeller
meta <- int.sct@meta.data
table(meta$cell_type, meta$sample)
prop.results <- propeller(clusters = meta$cell_type, sample = meta$sample, 
                          group = meta$condition)

plotCellTypeProps(clusters=meta$cell_type, sample=meta$condition) +
  scale_fill_manual(values = type.colors)




# stratified ----

#stratified <- readRDS(file = paste0(data.dir, "objects/stratified.rds"))
stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))

head(stratified[[]])
Idents(stratified) <- stratified$cell_type
DimPlot(stratified, label = T, label.box = T) + NoLegend()

## Rename identities ----

Idents(stratified)
stratified <- RenameIdents(stratified, 
                            `IFN responsive cells` = "Reactive cells",
#                            `AC suprabasal distal` = "",
#                            `AC basal proximal` = "",
                            `Cancerous TZ` = "Cancerous cells",
#                            `AC suprabasal proximal` = "",
#                            `AC proliferative basal proximal` = "",
#                            `TZ` = "",
                            `EMT cells` = "EMT-like cells",
                            `Wound Healing` = "Wounded cells",
                            `unlabeled` = "AC basal distal"
                           )
stratified$cell_type <- Idents(stratified)

# replace condition categories in case of alphabetical ordering of the conditions
Idents(stratified) <- stratified$condition
stratified <- RenameIdents(stratified, 
                            `Normal` = "A.Normal",
                            `Hyperplasia` = "B.Hyperplasia",
                            `Dysplasia` = "C.Dysplasia",
                            `Carcinoma` = "D.Carcinoma"
                           )

# convert main variables to vectors to avoid numeric conversion
# https://github.com/satijalab/seurat/issues/1508
stratified$orig.ident <- as.character(stratified$orig.ident)
stratified$sample <- as.character(stratified$sample)
stratified$condition <- as.character(Idents(stratified))
stratified$cell_type <- as.character(stratified$cell_type)

SaveH5Seurat(stratified, filename = paste0(data.dir, "stratified.h5Seurat"), overwrite = TRUE)
Convert(paste0(data.dir, "stratified.h5Seurat"), dest = "h5ad", overwrite = TRUE)

saveRDS(stratified, file = paste0(data.dir, "stratified_manuscript.rds"))


DimPlot(stratified, label = F, label.box = T, repel = T) + 
  ggtitle("Stratified cells")
ggsave(filename = paste0(results.dir, "stratified.identities.png"), width = 7, height = 4)


# August 2025

stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))

Idents(stratified) <- stratified$cell_type

DimPlot(stratified, label = F, label.box = T, repel = T) + 
  ggtitle("Stratified cells")
ggsave(filename = paste0(results.dir, "non.cherry.picked.umap.2B.png"), width = 8, height = 5)

DimPlot(stratified, 
        pt.size = 1.5,
        split.by = "condition",
        label = F, label.box = T, repel = T) +
  ggtitle("Stratified cells")
ggsave(filename = paste0(results.dir, "non.cherry.picked.umap.S3G.split.png"), width = 18, height = 5)

table(stratified$condition)
# A.Normal B.Hyperplasia   C.Dysplasia   D.Carcinoma 
#  151           304          2162          2661 



## proportions ----

table(stratified$cell_type)
Idents(stratified) <- stratified$cell_type

type.colors = alphabet(length(levels(as.factor(stratified$cell_type))))
show_col(type.colors)
names(type.colors) <- levels(as.factor(stratified$cell_type))

# Plot cell type proportions
tab1 <- table(stratified$condition, stratified$cell_type)
tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Cell_Type", "Condition","Proportion")
head(tab1)
ggplot(tab1, aes(x = Condition, y = Proportion, fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = type.colors) +
  ggtitle("Cell Type Proportions") +
  theme_ipsum() +
  xlab("") #+ RotatedAxis()
ggsave(paste0(results.dir, "celltype.proportions_stratified.png"), width=7, height=5, dpi=300)

# when having replicates (only hyperplasia), we can use propeller
meta <- stratified@meta.data
table(meta$cell_type, meta$sample)
table(meta$cell_type, meta$orig.ident)
prop.results <- propeller(clusters = meta$cell_type, sample = meta$orig.ident, 
                          group = meta$condition)

plotCellTypeProps(clusters=meta$cell_type, sample=meta$condition) +
  scale_fill_manual(values = type.colors)


## Scores ----

library(GO.db)
library(biomaRt)
library(Nebulosa)

stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
head(stratified[[]])
Idents(stratified) <- stratified$cell_type
DimPlot(stratified, label = T, label.box = T) + NoLegend()

DefaultAssay(stratified) <- "SCT"
#DefaultAssay(stratified) <- "RNA"


### IL17 ----
library(GO.db)

go_terms <- as.list(GOTERM)
keyword <- "17"
matches <- sapply(go_terms, function(term) {
  grepl(keyword, Term(term), ignore.case = TRUE)
})
matched_terms <- go_terms[matches]

# View matching GO terms
for (term in matched_terms) {
  cat(paste0(GOID(term), ": ", Term(term), "\n"))
}

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
go_id <- "GO:0097398"  # cellular response to interleukin-17
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = go_id,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "GO"
)
head(genes_in_go)
genes_in_go$SYMBOL
# "Socs3"    "Cxcl1"    "Cxcl10"   "Il1b"     "Il6"      "Nfkb1"    "Ccl1"     "Stat3"    
# "Nfkbiz"   "Nfkbiz"   "Traf3ip2" "Srsf1"    "Dcp1b"    "Mir409"   "Mir1896" 
int1 <- intersect(genes_in_go$SYMBOL, rownames(stratified))
# "Socs3"    "Cxcl1"    "Cxcl10"   "Il1b"     "Il6"      "Nfkb1"    "Stat3"    "Nfkbiz"   "Traf3ip2" "Srsf1"    "Dcp1b" 

stratified <- AddModuleScore(stratified, 
                                          features = list(genes_in_go$SYMBOL), 
                                          name = paste0("GO_","0097398"))
head(stratified)
# following level order as in Figure 2G
stratified$cell_type <- factor(stratified$cell_type, 
                                      levels = c("TZ","AC basal proximal","AC proliferative basal proximal",
                                                 "AC suprabasal proximal","AC basal distal","AC suprabasal distal",
                                                 "EMT-like cells","Wounded cells","Reactive cells",
                                                 "Cancerous cells"))
VlnPlot(stratified, paste0("GO_","0097398","1"), group.by = "cell_type", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") +
  ggtitle("GO:0097398 - Cellular response to interleukin-17")
#labs(title = "GO:0097398", subtitle = "Cellular response to interleukin-17")
ggsave(filename = paste0(results.dir, "violin_GO_0097398_stratified.png"), width = 12, height = 8)

FeaturePlot(stratified, genes_in_go$SYMBOL)
# stronger signal in cancerous cells: Socs3, Cxcl1, Nfkbiz
plot_density(stratified, paste0("GO_","0097398","1")) +
  labs(title = "GO:0097398", subtitle = "Cellular response to interleukin-17")
ggsave(filename = paste0(results.dir, "density_GO_0097398_stratified.png"), width = 6, height = 5)
plot_density(stratified, int1)
ggsave(filename = paste0(results.dir, "density_GO_0097398_genes_stratified.png"), width = 20, height = 15)
VlnPlot(stratified, c("Il17a","Junb"), group.by = "cell_type")

VlnPlot(stratified, paste0("GO_","0097398","1"), 
        group.by = "condition", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(0,1.4) +
  labs(title = "GO:0097398", subtitle = "Cellular response to interleukin-17") +
  stat_compare_means(label = "p.format", 
                 comparisons = list(c("A.Normal", "B.Hyperplasia"),
                                    c("B.Hyperplasia", "C.Dysplasia"),
                                    c("C.Dysplasia", "D.Carcinoma")))
ggsave(filename = paste0(results.dir, "violin_GO_0097398_condition_stratified.png"), width = 7, height = 6)


### TNF ----

go_terms <- as.list(GOTERM)
keyword <- "tumor necrosis factor"
matches <- sapply(go_terms, function(term) {
  grepl(keyword, Term(term), ignore.case = TRUE)
})
matched_terms <- go_terms[matches]
# View matching GO terms
for (term in matched_terms) {
  cat(paste0(GOID(term), ": ", Term(term), "\n"))
}
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
go_id <- "GO:0005031"  # tumor necrosis factor receptor activity
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = go_id,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "GO"
)
head(genes_in_go)
genes_in_go$SYMBOL
# [1] "Fas"       "Fas"       "Tnfrsf11a" "Tnfrsf11a" "Tnfrsf18"  "Tnfrsf1a"  
# "Tnfrsf1a"  "Tnfrsf1a"  "Tnfrsf1a"  "Tnfrsf1a"  "Tnfrsf1a"  "Tnfrsf1a"  
# "Tnfrsf1b"  "Tnfrsf1b"  "Tnfrsf1b"  "Tnfrsf4"   "Tnfrsf4"  "Eda2r"  
int1 <- intersect(genes_in_go$SYMBOL, rownames(stratified))
#  "Fas"       "Tnfrsf11a" "Tnfrsf18"  "Tnfrsf1a"  "Tnfrsf1b"  "Tnfrsf4"   "Eda2r"

stratified <- AddModuleScore(stratified, 
                             features = list(genes_in_go$SYMBOL), 
                             name = paste0("GO_","0005031"))
head(stratified)
VlnPlot(stratified, paste0("GO_","0005031","1"), group.by = "cell_type", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") +
  ggtitle("GO:0005031 - tumor necrosis factor receptor activity")
ggsave(filename = paste0(results.dir, "violin_GO_0005031_stratified.png"), width = 12, height = 8)
VlnPlot(stratified, paste0("GO_","0005031","1"), 
        group.by = "condition", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(0,0.7) +
  labs(title = "GO:0005031", subtitle = "tumor necrosis factor receptor activity") +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("A.Normal", "B.Hyperplasia"),
                                        c("B.Hyperplasia", "C.Dysplasia"),
                                        c("C.Dysplasia", "D.Carcinoma")))
ggsave(filename = paste0(results.dir, "violin_GO_0005031_condition_stratified.png"), width = 7, height = 6)



### IL6 ----

go_terms <- as.list(GOTERM)
keyword <- "interleukin-6"
matches <- sapply(go_terms, function(term) {
  grepl(keyword, Term(term), ignore.case = TRUE)
})
matched_terms <- go_terms[matches]
# View matching GO terms
for (term in matched_terms) {
  cat(paste0(GOID(term), ": ", Term(term), "\n"))
}
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
go_id <- "GO:0070102"  # interleukin-6-mediated signaling pathway
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = go_id,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "GO"
)
head(genes_in_go)
genes_in_go$SYMBOL
# "Cebpa" "Cebpa" "Fer"   "Fer"   "Gab1"  "Il6"   "Il6"   
# "Il6"   "Il6ra" "Il6ra" "Il6ra" "Il6st" "Il6st" "Jak1"  
# "Smad4" "Spi1"  "Src"   "Src"   "Stat3" "Stat3" "Ctr9"  
# "Yap1"  "Tut4"  "St18"  "St18"  
int1 <- intersect(genes_in_go$SYMBOL, rownames(stratified))
#  "Cebpa" "Fer"   "Gab1"  "Il6"   "Il6ra" "Il6st" "Jak1" 
#   "Smad4" "Spi1"  "Src"   "Stat3" "Ctr9"  "Yap1"  "Tut4" 

stratified <- AddModuleScore(stratified, 
                             features = list(genes_in_go$SYMBOL), 
                             name = paste0("GO_","0070102"))
#head(stratified)
VlnPlot(stratified, paste0("GO_","0070102","1"), group.by = "cell_type", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") +
  ggtitle("GO:0070102 - interleukin-6-mediated signaling pathway")
ggsave(filename = paste0(results.dir, "violin_GO_0070102_stratified.png"), width = 12, height = 8)
VlnPlot(stratified, paste0("GO_","0070102","1"), 
        group.by = "condition", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(0,0.8) +
  labs(title = "GO:0070102", subtitle = "interleukin-6-mediated signaling pathway") +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("A.Normal", "B.Hyperplasia"),
                                        c("B.Hyperplasia", "C.Dysplasia"),
                                        c("C.Dysplasia", "D.Carcinoma")))
ggsave(filename = paste0(results.dir, "violin_GO_0070102_condition_stratified.png"), width = 7, height = 6)




### TGFb ----

go_terms <- as.list(GOTERM)
keyword <- "transforming growth factor beta"
matches <- sapply(go_terms, function(term) {
  grepl(keyword, Term(term), ignore.case = TRUE)
})
matched_terms <- go_terms[matches]
# View matching GO terms
for (term in matched_terms) {
  cat(paste0(GOID(term), ": ", Term(term), "\n"))
}
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
go_id <- "GO:0005024"  # transforming growth factor beta receptor activity
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = go_id,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "GO"
)
head(genes_in_go)
genes_in_go$SYMBOL
# "Acvrl1"  "Bmpr1b"  "Bmpr2"   "Tgfbr1"  "Tgfbr1"  "Tgfbr1"  "Tgfbr2"  "Tgfbr2"  
# "Tgfbr2"  "Tgfbr3"  "Tgfbr3"  "Tgfbr3"  "Tgfbr3"  "Amhr2"   "Amhr2"   "Tgfbr3l"
int1 <- intersect(genes_in_go$SYMBOL, rownames(stratified))
#  "Acvrl1"  "Bmpr2"   "Tgfbr1"  "Tgfbr2"  "Tgfbr3"  "Amhr2"   "Tgfbr3l"

stratified <- AddModuleScore(stratified, 
                             features = list(genes_in_go$SYMBOL), 
                             name = paste0("GO_","0005024"))
#head(stratified)
VlnPlot(stratified, paste0("GO_","0005024","1"), group.by = "cell_type", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") +
  ggtitle("GO:0005024 - transforming growth factor beta receptor activity")
ggsave(filename = paste0(results.dir, "violin_GO_0005024_stratified.png"), width = 12, height = 8)
VlnPlot(stratified, paste0("GO_","0005024","1"), 
        group.by = "condition", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(0,0.8) +
  labs(title = "GO:0005024", subtitle = "transforming growth factor beta receptor activity") +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("A.Normal", "B.Hyperplasia"),
                                        c("B.Hyperplasia", "C.Dysplasia"),
                                        c("C.Dysplasia", "D.Carcinoma")))
ggsave(filename = paste0(results.dir, "violin_GO_0005024_condition_stratified.png"), width = 7, height = 6)


saveRDS(stratified, file = paste0(data.dir, "stratified_manuscript.rds"))


## DEGs ----

stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
head(stratified[[]])
table(stratified$cell_type)
Idents(stratified) <- stratified$cell_type
DimPlot(stratified, label = T, label.box = T) + NoLegend()

### conserved markers ----

stratified <- PrepSCTFindMarkers(stratified)
#DefaultAssay(stratified) <- "SCT"
# use more selective markers to distinguish relevant cell types
celltypes.of.interest <- c("EMT-like cells", "Wounded cells", "Reactive cells")
i = 2
sel.celltype <- celltypes.of.interest[i]
con.markers <- FindConservedMarkers(stratified,
                                    only.pos = F, 
                                    min.pct = 0.2, 
                                    logfc.threshold = 0.5,
                                    ident.1 = sel.celltype,
                                    grouping.var = "sample")
con.markers <- con.markers[con.markers$max_pval < 0.05, ]
head(con.markers)
table(con.markers$Carcinoma_avg_log2FC > 0)
write.csv(con.markers, file = paste0(results.dir, "conserved.markers_", sel.celltype, ".csv"))



markers <- FindAllMarkers(stratified, 
                          only.pos = T, 
                          min.pct = 0.2, 
                          logfc.threshold = 0.5,
                          assay = "SCT")
head(markers)
markers <- markers[markers$p_val_adj < 0.05, ]
table(markers$cluster)
write.csv(markers, file = paste0(results.dir, "stratified.markers.csv"))


# split tables 

markers <- read.csv(file = paste0(results.dir, "stratified.markers.csv"), row.names = 1)
table(markers$cluster)

for(i in unique(markers$cluster)) {
  markers.sub <- markers[markers$cluster == i, ]
  write.csv(markers.sub, file = paste0(results.dir, "stratified.markers.", i, ".csv"))
}




### Cancerous cells ----

markers.cancerous <- markers[markers$cluster == "Cancerous cells", ] # 492
write.csv(markers.cancerous, file = paste0(results.dir, "stratified.markers.cancerous.csv"))

saveRDS(stratified, file = paste0(data.dir, "stratified_manuscript.rds"))


# plot scores from previous section (as in Figure 2g)

cancerous <- subset(stratified, cell_type == "Cancerous cells")
table(cancerous$sample)
head(cancerous)
Idents(cancerous) <- cancerous$sample
cancerous <- subset(cancerous, idents = c("Dysplasia","Carcinoma"))
cancerous$sample <- factor(cancerous$sample, 
                                        levels = c("Dysplasia", "Carcinoma"))

VlnPlot(cancerous, paste0("GO_","0097398","1"), 
        group.by = "sample", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(0,1) +
  labs(title = "GO:0097398", subtitle = "Cellular response to interleukin-17") +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("Dysplasia", "Carcinoma")))
ggsave(filename = paste0(results.dir, "violin_GO_0097398_condition_cancerous.png"), width = 6, height = 6)

VlnPlot(cancerous, paste0("GO_","0005031","1"), 
        group.by = "sample", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
#  xlab("") + ylim(0,5) +
  labs(title = "GO:0005031", subtitle = "tumor necrosis factor receptor activity") +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("Dysplasia", "Carcinoma")))
ggsave(filename = paste0(results.dir, "violin_GO_0005031_condition_cancerous.png"), width = 7, height = 6)

VlnPlot(cancerous, paste0("GO_","0070102","1"), 
        group.by = "sample", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(0,0.8) +
  labs(title = "GO:0070102", subtitle = "interleukin-6-mediated signaling pathway") +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("Dysplasia", "Carcinoma")))
ggsave(filename = paste0(results.dir, "violin_GO_0070102_condition_cancerous.png"), width = 7, height = 6)

VlnPlot(cancerous, paste0("GO_","0005024","1"), 
        group.by = "sample", sort = F) +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(0,0.8) +
  labs(title = "GO:0005024", subtitle = "transforming growth factor beta receptor activity") +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("Dysplasia", "Carcinoma")))
ggsave(filename = paste0(results.dir, "violin_GO_0005024_condition_cancerous.png"), width = 7, height = 6)



### cell subtype signatures ----

# to distinguish EMT-like, Reactive and Wounded cells
# EMT is the top hallmark in EMT-like cluster markers
# to identify top markers, we used conserved markers (above)
# and selected top pathways using EnrichR
# wounded cells only had negative markers with this strategy,
# so we used a published wound healing signature (https://www.nature.com/articles/ncomms14684)

reactive.genes <- list(
  Glutathione = c("Gclc","Gsta4","Gsta3"),
  Oxidations = c("Aldh3a1","Gclc","Gsta4","Gsta3","Fmo2"),
  Chemical_stress = c("Cat","Gclc","Gsta3"),
  Cornification = c("Krt13","Krt12","Dsc2"),
  Epithelial_differentiation = c("Krtdap","Crabp2","Fabp5","Krt14","Emp1","Krt10",
                                 "Lgals3","Krt17","Edf1","Tagln2","Sult2b1")
  )
EMT.like.genes <- list(
  EMT = c("Serpine2","Mgp","Snai2","Fbln2"),
  ECM = c("Tgfb2","Ltbp4","Pdgfa","Agrn","Fbln2")
  )
wound_healing_signature <- list(
  phase_I = c("Il24", "Il33", "S100a8", "S100a9"), # Inflammatory phase
  phase_II = c("Itga5", "Itga6", "Plau", "Plaur", "Ephb2", "Efnb1", "Cxcr4", 
               "C5ar1", "Myh9", "Procr", "Wnt5a", "Elk3", "Hmga2", "Inhba", 
               "Pcdh7", "Pcdhb19", "Pcdhga", "Fn1", "Lama3", "Lamb3", "Lamc2",
               "Cdsn", "Fscn1", "Cald1", "Nav2", "Fmnl2", 
               "Gjb6","Cx30", "Gjb2","Cx26", "Myo1b", 
               "Myo5b", "Myh9", 
               "Tpm1", "Tpm2",
               "Tubb2a", "Tubb3","Tubb6",
               "Gprc5a","E2f7","Fgf18","Aim2","TnC",
               "Krt17","Krt6","Fos","Junb","Egr1","Egr2",
               "E2f7","Fgf18","Dsc2","Areg","Ereg",
               "Emb","Epgn","Fscn1",
               "Emb","Sprr1b","Sprr2h",
               "Flrt2","Myo1b"), # Regenerative phase
  phase_III = c("Mmp9", "Mmp13", "Mmp1b", 
                "Timp3", "Slc2a1", "Hk2", "Pgk1") # Remodeling phase
)

sel.gene.sets <- list(
  'Reactive genes' = unique(unlist(reactive.genes)),
  'EMT-like genes' = unique(unlist(EMT.like.genes)),
  'Wound healing genes' = unique(unlist(wound_healing_signature))
)

names(sel.gene.sets)

stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
head(stratified[[]])
stratified$sample <- factor(stratified$sample, levels = c("Normal","Hyperplasia","Dysplasia","Carcinoma"))

stratified <- AddModuleScore(stratified, sel.gene.sets, name = names(sel.gene.sets))
head(stratified)

Idents(stratified) <- stratified$cell_type
stratified.subset <- subset(stratified, idents = c("EMT-like cells", "Wounded cells", "Reactive cells"))

sel.genes <- c("Krt13","Gsta4","Lgals3","Tagln2","Sult2b1",
               "Serpine2","Mgp","Snai2","Pdgfa","Fbln2",
               "Itga6","Myh9","Junb","Egr1","Slc2a1")

DotPlot(stratified,
        idents = c("EMT-like cells", "Wounded cells", "Reactive cells"),
        cols = c("blue","red"),
        features = sel.genes) + 
  RotatedAxis() + xlab("") + ylab("")
ggsave(filename = paste0(results.dir, "dotplot_stratified_subtypes.png"), 
         width = 10, height = 6, dpi = 300)

i=1
VlnPlot(stratified.subset, paste0(names(sel.gene.sets)[i],i), 
        group.by = "cell_type", sort = F) + NoLegend() +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") +
  ggtitle(names(sel.gene.sets)[i]) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text=element_text(size=10))
ggsave(filename = paste0(results.dir, "violin_", 
                           names(sel.gene.sets)[i], "_stratified_subtypes.png"), 
         width = 4, height = 5)
#plot_density(stratified, paste0(names(sel.gene.sets)[i],i), method = "wkde")
VlnPlot(stratified, paste0(names(sel.gene.sets)[i],i), 
#        split.by = "cell_type",
        group.by = "sample", sort = F) + NoLegend() +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") +
  ylim(0,1) +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("Normal", "Hyperplasia"),
                                        c("Hyperplasia", "Dysplasia"),
                                        c("Dysplasia", "Carcinoma"))) +
  ggtitle(names(sel.gene.sets)[i]) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text=element_text(size=10))
ggsave(filename = paste0(results.dir, "violin_", 
                           names(sel.gene.sets)[i], "_stratified_subtypes_condition.png"), 
         width = 6, height = 6)


i=1
v1 <- VlnPlot(stratified.subset, paste0(names(sel.gene.sets)[i],i), 
        group.by = "cell_type", sort = F) + NoLegend() +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(-0.5,2) +
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("EMT-like cells", "Wounded cells"),
                                        c("Wounded cells", "Reactive cells"),
                                        c("EMT-like cells", "Reactive cells"))) +
  ggtitle(names(sel.gene.sets)[i]) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text=element_text(size=10))
i=2
v2 <- VlnPlot(stratified.subset, paste0(names(sel.gene.sets)[i],i), 
              group.by = "cell_type", sort = F) + NoLegend() +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(-0.6,2) +
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("EMT-like cells", "Wounded cells"),
                                        c("Wounded cells", "Reactive cells"),
                                        c("EMT-like cells", "Reactive cells"))) +
  ggtitle(names(sel.gene.sets)[i]) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text=element_text(size=10))
i=3
v3 <- VlnPlot(stratified.subset, paste0(names(sel.gene.sets)[i],i), 
              group.by = "cell_type", sort = F) + NoLegend() +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") + ylim(-0.2,0.45) +
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("EMT-like cells", "Wounded cells"),
                                        c("Wounded cells", "Reactive cells"),
                                        c("EMT-like cells", "Reactive cells"))) +
  ggtitle(names(sel.gene.sets)[i]) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text=element_text(size=10))

ga1 <- grid.arrange(v1,v2,v3, ncol = 3)
ggsave(filename = paste0(results.dir, "violin_stratified_subtypes.png"), 
       ga1,
         width = 12, height = 5, dpi = 300)





saveRDS(stratified, file = paste0(data.dir, "stratified_manuscript.rds"))


#### wound healing ----
# in addition to global wound healing above,
# wound healing by phase

wound_healing_signature <- list(
  WH.phase_I = c("Il24", "Il33", "S100a8", "S100a9"), # Inflammatory phase
  WH.phase_II = c("Itga5", "Itga6", "Plau", "Plaur", "Ephb2", "Efnb1", "Cxcr4", 
               "C5ar1", "Myh9", "Procr", "Wnt5a", "Elk3", "Hmga2", "Inhba", 
               "Pcdh7", "Pcdhb19", "Pcdhga", "Fn1", "Lama3", "Lamb3", "Lamc2",
               "Cdsn", "Fscn1", "Cald1", "Nav2", "Fmnl2", 
               "Gjb6","Cx30", "Gjb2","Cx26", "Myo1b", 
               "Myo5b", "Myh9", 
               "Tpm1", "Tpm2",
               "Tubb2a", "Tubb3","Tubb6",
               "Gprc5a","E2f7","Fgf18","Aim2","TnC",
               "Krt17","Krt6","Fos","Junb","Egr1","Egr2",
               "E2f7","Fgf18","Dsc2","Areg","Ereg",
               "Emb","Epgn","Fscn1",
               "Emb","Sprr1b","Sprr2h",
               "Flrt2","Myo1b"), # Regenerative phase
  WH.phase_III = c("Mmp9", "Mmp13", "Mmp1b", 
                "Timp3", "Slc2a1", "Hk2", "Pgk1") # Remodeling phase
)

names(wound_healing_signature)
stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
head(stratified[[]])
stratified <- AddModuleScore(stratified, wound_healing_signature, name = names(wound_healing_signature))
head(stratified)
i=3
VlnPlot(stratified, paste0(names(wound_healing_signature)[i],i), 
        #        split.by = "cell_type",
        group.by = "sample", sort = F) + NoLegend() +
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white") +
  xlab("") +
  ylim(-0.3,1.5) +
  stat_compare_means(label = "p.format", 
                     comparisons = list(c("Normal", "Hyperplasia"),
                                        c("Hyperplasia", "Dysplasia"),
                                        c("Dysplasia", "Carcinoma"))) +
  ggtitle(names(wound_healing_signature)[i]) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text=element_text(size=10))
ggsave(filename = paste0(results.dir, "violin_wound.healing_phase.", 
                         names(wound_healing_signature)[i], "_stratified_subtypes_condition.png"), 
       width = 6, height = 6)



## subcluster stratified ----
# to better discriminate cancerous, EMT-like, wounded, and reactive

obj <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
table(obj$cell_type)
Idents(obj) <- obj$cell_type
obj <- subset(obj, idents = c("Cancerous cells", "Wounded cells",
                              "EMT-like cells", "Reactive cells"))
DimPlot(obj, label = T, label.box = T) + NoLegend()
obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj)
ElbowPlot(obj)
obj <- FindNeighbors(obj, dims = 1:10, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(obj, label = T, label.box = T) + NoLegend()
DimPlot(obj, group.by = "cell_type", label = T, label.box = T) + NoLegend()

# while reactive and cancerous cells remain distinct,
# EMT-like and wounded cells are rather mixed

DimPlot(obj, label = T, label.box = T,
        split.by = "cell_type") + NoLegend()
DimPlot(obj, split.by = "condition")
DimPlot(obj,
        group.by = "cell_type",
        split.by = "condition")

markers <- FindAllMarkers(obj, 
                          only.pos = T, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25,
                          assay = "SCT")
table(markers$cluster)
for(i in unique(markers$cluster)) {
  markers.sub <- markers[markers$cluster == i, ]
  write.csv(markers.sub, file = paste0(results.dir, "subclustering_stratified/markers.", i, ".csv"))
}




# immune ----

#immune <- readRDS(file = paste0(data.dir, "immune.rds"))
immune <- readRDS(file = paste0(data.dir, "immune_Atlas.anno.rds"))
head(immune[[]])

Idents(immune) <- immune$new.cell.labels
Idents(immune) <- immune$customclassif
Idents(immune) <- immune$seurat_clusters
DimPlot(immune, label = T, label.box = T) + NoLegend()

# I don't have the last version of the object
# but "new.cell.labels" seem to match the last version of categories
# some changes need to be made
# remove "cancer cells" from neutrophil cluster
# add Naive B cells from "customclassif" or "seurat_clusters"
# group gd and T into "T cells"

# is the colonocyte-like cluster here?
colo.like.genes <- c("Epcam", "Krt8", "Krt18", "Krt20a", "Cxcr2", 
                "Cd24a", "Fcgr3", "S100a8", "S100a9", "Muc3", "Car4", "Saa1")
FeaturePlot(immune, 
            features = colo.like.genes, 
            ncol = 3)
immune <- AddModuleScore(immune, list(colo.like.genes))
VlnPlot(immune, features = "Cluster1", group.by = "seurat_clusters")
ident.5 <- WhichCells(immune, idents = "5") # 443

Bcells <- WhichCells(immune, expression = seurat_clusters == "6") # 301 cells
to.remove <- WhichCells(immune, expression = new.cell.labels == "Neutrophils" & customclassif == "Cancer cells") # 262 cells

immune <- subset(immune, cells = to.remove, invert = T)
Idents(immune) <- immune$new.cell.labels
immune <- SetIdent(immune, cells = Bcells, value = "Naive B cells")
immune <- RenameIdents(immune, 
                            `gd T cells` = "T cells",
                            `T cells` = "T cells"
                           )

# replace condition categories in case of alphabetical ordering of the conditions
table(immune$sample)
immune$condition <- immune$sample
Idents(immune) <- immune$condition
immune <- RenameIdents(immune, 
                           `Normal` = "A.Normal",
                           `Hyperplasia` = "B.Hyperplasia",
                           `Dysplasia` = "C.Dysplasia",
                           `Carcinoma` = "D.Carcinoma"
)

# convert main variables to vectors to avoid numeric conversion
# https://github.com/satijalab/seurat/issues/1508
immune$orig.ident <- as.character(immune$orig.ident)
immune$sample <- as.character(immune$sample)
immune$condition <- as.character(Idents(immune))
immune$cell_type <- as.character(immune$cell_type)

SaveH5Seurat(immune, filename = paste0(data.dir, "immune.h5Seurat"), overwrite = TRUE)
Convert(paste0(data.dir, "immune.h5Seurat"), dest = "h5ad", overwrite = TRUE)

saveRDS(immune, file = paste0(data.dir, "immune_manuscript.rds"))


Idents(immune) <- immune$cell_type
DimPlot(immune, label = F, label.box = T, repel = T) + 
  ggtitle("Immune cells")
ggsave(filename = paste0(results.dir, "immune.identities.png"), width = 7, height = 6)


# August 2025

immune <- readRDS(file = paste0(data.dir, "immune_manuscript.rds"))

Idents(immune) <- immune$cell_type

DimPlot(immune, label = F, label.box = T, repel = T) + 
  ggtitle("Immune cells")
ggsave(filename = paste0(results.dir, "non.cherry.picked.umap.3A.png"), width = 7, height = 6)

DimPlot(immune, 
        split.by = "condition",
        label = F, label.box = T, repel = T) +
  ggtitle("Immune cells")
ggsave(filename = paste0(results.dir, "non.cherry.picked.umap.S6C.split.png"), width = 16, height = 5)

table(immune$condition)
# A.Normal B.Hyperplasia   C.Dysplasia   D.Carcinoma 
# 1318           492          1710          1954




## proportions ----

table(immune$cell_type)
#Idents(immune) <- immune$cell_type

type.colors = alphabet(length(levels(as.factor(immune$cell_type))))
show_col(type.colors)
names(type.colors) <- levels(as.factor(immune$cell_type))

# Plot cell type proportions
tab1 <- table(immune$condition, immune$cell_type)
tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Cell_Type", "Condition","Proportion")
head(tab1)
ggplot(tab1, aes(x = Condition, y = Proportion, fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = type.colors) +
  ggtitle("Cell Type Proportions") +
  theme_ipsum() +
  xlab("") #+ RotatedAxis()
ggsave(paste0(results.dir, "celltype.proportions_immune.png"), width=7, height=5, dpi=300)

# when having replicates (only hyperplasia), we can use propeller
meta <- immune@meta.data
table(meta$cell_type, meta$sample)
table(meta$cell_type, meta$orig.ident)
prop.results <- propeller(clusters = meta$cell_type, sample = meta$orig.ident, 
                          group = meta$condition)

plotCellTypeProps(clusters=meta$cell_type, sample=meta$condition) +
  scale_fill_manual(values = type.colors)



table(immune$orig.ident, immune$condition)

props <- getTransformedProps(immune$cell_type, immune$orig.ident, transform="logit")
barplot(props$Proportions, col = c("orange","purple","dark green"),legend=TRUE, 
        ylab="Proportions")
names(props)
props$TransformedProps

group <- c("C","B","B","A","D")
design <- model.matrix(~ 0 + group)
design

propeller.anova(prop.list=props, design=design, coef = c(1,2,3,4), 
                robust=TRUE, trend=FALSE, sort=TRUE)
library(limma)
mycontr <- makeContrasts(groupB-groupA, levels=design)
mycontr <- makeContrasts(groupC-groupB, levels=design)
mycontr <- makeContrasts(groupD-groupC, levels=design)
mycontr <- makeContrasts(groupC-groupA, levels=design)
mycontr <- makeContrasts(groupD-groupA, levels=design)
propeller.ttest(props, 
                design, 
                contrasts = mycontr, 
                robust=TRUE, trend=FALSE, 
                sort=TRUE)

progress <- c(3,2,2,1,4)
des.progress <- model.matrix(~progress)
des.progress

fit <- lmFit(props$TransformedProps,des.progress)
fit <- eBayes(fit, robust=TRUE)
topTable(fit)


## other cytokines ----

library(GO.db)
library(biomaRt)
library(Nebulosa)

immune <- readRDS(file = paste0(data.dir, "immune_manuscript.rds"))
head(immune[[]])
Idents(immune) <- immune$cell_type
DimPlot(immune, label = T, label.box = T) + NoLegend()

VlnPlot(immune, c("Il6","Il17a","Tnf","Tgfb1"), group.by = "cell_type", ncol = 2) &
  xlab("")
ggsave(filename = paste0(results.dir, "violin_IL6_IL17_TNF_TGFB1_immune.png"), width = 8, height = 7)

# use format as in Figure 3b
immune <- RenameIdents(immune, 
                            `T cells` = "T cells",
                            `Neutrophils` = "Neutrophils",
                            `Myeloid cells` = "Myeloid",
                            `ILCs` = "ILC",
                            `Naive B cells` = "Naive B",
                            `B cells` = "B cells"
                           )
immune$cell_type <- Idents(immune)
levels(as.factor(immune$cell_type))

DefaultAssay(immune) <- "SCT"
#DefaultAssay(immune) <- "integrated"
#DefaultAssay(immune) <- "RNA"

VlnPlot(immune, c("Il6","Tnf","Tgfb1"),
        pt.size = 0.1,
        group.by = "cell_type", ncol = 3) &
  xlab("") &
  geom_boxplot(aes(ymin = -Inf, ymax = Inf), 
               width = 0.3, alpha = 0.2, color = "white")
ggsave(filename = paste0(results.dir, "violin_IL6_TNF_TGFB1_immune.png"), width = 12, height = 4)





# T cells ----

load(file = paste0(data.dir, "objects_for_cellchat/t_ILC.RData"))
t_ILC
DimPlot(t_ILC)



saveRDS(t_ILC, file = paste0(data.dir, "t_ILC_manuscript.rds"))

 
# compositional analyses ----

## sscomp ----
# https://github.com/MangiolaLaboratory/sccomp

library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(tidyr)


### stratified ----

# using group as continuous
# similar results in pairwise categorial comparisons
stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
head(stratified)
stratified$sample <- factor(stratified$sample, levels = c("Normal","Hyperplasia","Dysplasia","Carcinoma"))

DimPlot(stratified, group.by = "cell_type")
DimPlot(stratified, #pt.size = 0.2,
        split.by = "sample",
        group.by = "cell_type") +
  ggtitle("")
ggsave(filename = paste0(results.dir, "stratified.identities_split.by.group.png"), width = 12, height = 4)
DimPlot(stratified, #pt.size = 0.2,
        group.by = "sample") +
  ggtitle("")
ggsave(filename = paste0(results.dir, "stratified.identities_grouped.by.stage.png"), width = 7, height = 6)

Idents(stratified) <- stratified$sample
stratified <- RenameIdents(stratified, 
                            `Normal` = "1",
                            `Hyperplasia` = "2",
                            `Dysplasia` = "3",
                            `Carcinoma` = "4"
)

stratified$type.continuous <- Idents(stratified)
stratified$type.continuous <- as.numeric(stratified$type.continuous)
table(stratified$sample)
table(stratified$type.continuous)

sccomp_result = 
  stratified |>
  sccomp_estimate( 
#    formula_composition = ~ 0 + sample, 
    formula_composition = ~ type.continuous, 
    sample = "orig.ident", 
    cell_group = "cell_type", 
    cores = 8
  ) 
#sccomp_test(sccomp_result, contrasts = c("sampleCarcinoma - sampleNormal"))
sccomp_result_outliers <- sccomp_remove_outliers(sccomp_result)
#sccomp_result_test <- sccomp_test(sccomp_result, contrasts = "type.continuous")
sccomp_result_test <- sccomp_test(sccomp_result_outliers, contrasts = "type.continuous")

plot(sccomp_result_test)
ggsave(filename = paste0(results.dir, "sccomp_stratified.png"), width = 5, height = 5)
#sccomp_boxplot(sccomp_result_test, factor = "sample")


# Dimplot with same cell number

table(stratified$sample)
normal.cells <- WhichCells(stratified, expression = sample == "Normal") # 151
hyper.cells <-  sample(WhichCells(stratified, expression = sample == "Hyperplasia"), length(normal.cells))
dys.cells <-  sample(WhichCells(stratified, expression = sample == "Dysplasia"), length(normal.cells))
cancer.cells <-  sample(WhichCells(stratified, expression = sample == "Carcinoma"), length(normal.cells))

stratified.151 <- subset(stratified, cells = c(normal.cells,
                                               hyper.cells,
                                               dys.cells,
                                               cancer.cells))
DimPlot(stratified.151, pt.size = 2, alpha = 0.5,
        split.by = "sample",
        group.by = "cell_type") +
  ggtitle("")
ggsave(filename = paste0(results.dir, "stratified.151.identities_split.by.group.png"), width = 12, height = 4)
DimPlot(stratified.151, pt.size = 2, alpha = 0.7,
        group.by = "sample") +
  ggtitle("")
ggsave(filename = paste0(results.dir, "stratified.151.identities_grouped.by.stage.png"), width = 7, height = 6)



saveRDS(stratified, file = paste0(data.dir, "stratified_manuscript.rds"))


### immune ----

immune <- readRDS(file = paste0(data.dir, "immune_manuscript.rds"))
immune$sample <- factor(immune$sample, levels = c("Normal","Hyperplasia","Dysplasia","Carcinoma"))
DimPlot(immune)

Idents(immune) <- factor(immune$cell_type, 
                         levels = c("T cells","Neutrophils","Myeloid cells","ILCs","Naive B cells","B cells"))
immune$cell_type <- Idents(immune)
DimPlot(immune, #pt.size = 0.2,
        split.by = "sample",
        group.by = "cell_type") +
  ggtitle("")
ggsave(filename = paste0(results.dir, "immune.identities_split.by.group.png"), width = 12, height = 4)
DimPlot(immune, #pt.size = 0.2,
        group.by = "sample") +
  ggtitle("")
ggsave(filename = paste0(results.dir, "immune.identities_grouped.by.stage.png"), width = 7, height = 6)


Idents(immune) <- immune$sample
immune <- RenameIdents(immune, 
                           `Normal` = "1",
                           `Hyperplasia` = "2",
                           `Dysplasia` = "3",
                           `Carcinoma` = "4"
)

immune$type.continuous <- Idents(immune)
immune$type.continuous <- as.numeric(immune$type.continuous)
table(immune$sample)
table(immune$type.continuous)

sccomp_result = 
  immune |>
  sccomp_estimate( 
    #    formula_composition = ~ 0 + sample, 
    formula_composition = ~ type.continuous, 
    sample = "orig.ident", 
    cell_group = "cell_type", 
    cores = 8
  ) 
#sccomp_test(sccomp_result, contrasts = c("sampleCarcinoma - sampleNormal"))
sccomp_result_outliers <- sccomp_remove_outliers(sccomp_result)
#sccomp_result_test <- sccomp_test(sccomp_result, contrasts = "type.continuous")
sccomp_result_test <- sccomp_test(sccomp_result_outliers, contrasts = "type.continuous")

plot(sccomp_result_test)
ggsave(filename = paste0(results.dir, "sccomp_immune.png"), width = 5, height = 5)
#sccomp_boxplot(sccomp_result_test, factor = "sample")



# Dimplot with same cell number

table(immune$sample)
hyper.cells <-  WhichCells(immune, expression = sample == "Hyperplasia") # 492
normal.cells <- sample(WhichCells(immune, expression = sample == "Normal"), length(hyper.cells))
dys.cells <-  sample(WhichCells(immune, expression = sample == "Dysplasia"), length(hyper.cells))
cancer.cells <-  sample(WhichCells(immune, expression = sample == "Carcinoma"), length(hyper.cells))

immune.492 <- subset(immune, cells = c(normal.cells,
                                               hyper.cells,
                                               dys.cells,
                                               cancer.cells))
DimPlot(immune.492, pt.size = 1, alpha = 0.5,
        split.by = "sample",
        group.by = "cell_type") +
  ggtitle("")
ggsave(filename = paste0(results.dir, "immune.492.identities_split.by.group.png"), width = 12, height = 4)
DimPlot(immune.492, pt.size = 1, alpha = 0.7,
        group.by = "sample") +
  ggtitle("")
ggsave(filename = paste0(results.dir, "immune.492.identities_grouped.by.stage.png"), width = 7, height = 6)


saveRDS(immune, file = paste0(data.dir, "immune_manuscript.rds"))


## miloR ----
# as suggested by the August reviews
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7617075/
# https://github.com/MarioniLab/miloR?tab=readme-ov-file

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)


### stratified ----

stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
head(stratified)
stratified$sample <- factor(stratified$sample, levels = c("Normal","Hyperplasia","Dysplasia","Carcinoma"))

DimPlot(stratified, group.by = "cell_type")
DimPlot(stratified, #pt.size = 0.2,
        split.by = "sample",
        group.by = "cell_type") +
  ggtitle("")

pca_embeddings <- stratified[["pca"]]@cell.embeddings
umap_embeddings <- stratified[["umap"]]@cell.embeddings
sce <- as.SingleCellExperiment(stratified)
reducedDim(sce, "PCA", withDimnames = TRUE) <- pca_embeddings
reducedDim(sce, "UMAP", withDimnames = TRUE) <- umap_embeddings
reducedDimNames(sce)

#sce <- sce[,apply(reducedDim(sce, "PCA"), 1, function(x) !all(is.na(x)))]
#sce <- runUMAP(sce, dimred = "PCA", name = 'umap')
#plotReducedDim(sce, colour_by="condition", dimred = "umap") 

milo <- Milo(sce)
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d=30, 
                   refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(milo)

milo <- countCells(milo, 
                   meta.data = as.data.frame(colData(milo)), 
                   sample="orig.ident")
head(nhoodCounts(milo))

milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")

design <- data.frame(colData(milo))[,c("orig.ident", "condition", "type.continuous")]
head(design)
# changing to numeric doesn't seem to change the result:
design$type.continuous <- as.numeric(factor(design$condition))
table(design$type.continuous)
table(design$condition)
design <- distinct(design)
rownames(design) <- design$orig.ident
design

da_results <- testNhoods(milo, 
#                         design = ~condition, 
                         design = ~type.continuous, 
#                         design = ~as.numeric(type.continuous), 
                         design.df=design[colnames(nhoodCounts(milo)), ])
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

milo <- buildNhoodGraph(milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = "UMAP", colour_by="condition", text_by = "cell_type", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
ggsave(filename = paste0(results.dir, "milo_stratified_umap_nhoods.png"), 
       width = 12, height = 6, dpi = 300)

da_results <- annotateNhoods(milo, da_results, coldata_col = "cell_type")
head(da_results)

ggplot(da_results, aes(cell_type_fraction)) + geom_histogram(bins=50)

da_results$celltype <- ifelse(da_results$cell_type_fraction < 0.7, "Mixed", da_results$cell_type)
da_results$celltype <- factor(da_results$celltype, 
                              levels = c("Mixed",
                                         "AC suprabasal distal",
                                         "TZ",  
                                         "Wounded cells",
                                         "AC basal distal",
                                         "AC suprabasal proximal",
                                         "AC proliferative basal proximal",
                                         "Reactive cells",
                                         "EMT-like cells",
                                         "AC basal proximal", 
                                         "Cancerous cells"
                              ))
plotDAbeeswarm(da_results, group.by = "celltype")
ggsave(filename = paste0(results.dir, "milo_stratified_beeswarm.png"), 
       width = 9, height = 10, dpi = 300)

write.csv(da_results, file = paste0(results.dir, "milo_stratified_results.csv"))




### immune ----

immune <- readRDS(file = paste0(data.dir, "immune_manuscript.rds"))
DimPlot(immune, #pt.size = 0.2,
        split.by = "sample",
        group.by = "cell_type") +
  ggtitle("")

pca_embeddings <- immune[["pca"]]@cell.embeddings
umap_embeddings <- immune[["umap"]]@cell.embeddings
sce <- as.SingleCellExperiment(immune)
reducedDim(sce, "PCA", withDimnames = TRUE) <- pca_embeddings
reducedDim(sce, "UMAP", withDimnames = TRUE) <- umap_embeddings
reducedDimNames(sce)

milo <- Milo(sce)
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d=30, 
                   refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(milo)

milo <- countCells(milo, 
                   meta.data = as.data.frame(colData(milo)), 
                   sample="orig.ident")
head(nhoodCounts(milo))

milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")

design <- data.frame(colData(milo))[,c("orig.ident", "condition", "type.continuous")]
head(design)
table(design$type.continuous)
table(design$condition)
design <- distinct(design)
rownames(design) <- design$orig.ident
design

da_results <- testNhoods(milo, 
                         design = ~type.continuous, 
                         design.df=design[colnames(nhoodCounts(milo)), ])
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

milo <- buildNhoodGraph(milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = "UMAP", colour_by="condition", text_by = "cell_type", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
ggsave(filename = paste0(results.dir, "milo_immune_umap_nhoods.png"), 
       width = 12, height = 6, dpi = 300)

da_results <- annotateNhoods(milo, da_results, coldata_col = "cell_type")
head(da_results)

ggplot(da_results, aes(cell_type_fraction)) + geom_histogram(bins=50)

da_results$celltype <- ifelse(da_results$cell_type_fraction < 0.7, "Mixed", da_results$cell_type)
da_results$celltype <- factor(da_results$celltype, 
                              levels = c("Mixed",
                                         "Naive B cells",
                                         "Myeloid cells",
                                         "ILCs",
                                         "T cells",
                                         "B cells",
                                         "Neutrophils"
                              ))
plotDAbeeswarm(da_results, group.by = "celltype")
ggsave(filename = paste0(results.dir, "milo_immune_beeswarm.png"), 
       width = 7, height = 7, dpi = 300)

write.csv(da_results, file = paste0(results.dir, "milo_immune_results.csv"))




# Trajectory ------------------------------------------------

library(tidyr)
library(slingshot)
library(Lamian)
library(condiments)
library(tradeSeq)
library(UpSetR)
library(SingleCellExperiment)
library(scales)
library(ggplot2)
library(Seurat)
library(pheatmap)

prepare_data_for_trajectory <- function(seurat_obj, time_var = "timepoint", 
                                        celltype_var = "cell_type") {
  
  # Convert Seurat to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj)
  
  # Add metadata
  colData(sce)$timepoint <- seurat_obj@meta.data[[time_var]]
  colData(sce)$cell_type <- seurat_obj@meta.data[[celltype_var]]
  
  # Add dimensionality reductions
  reducedDim(sce, "PCA") <- Embeddings(seurat_obj, "pca")
  reducedDim(sce, "UMAP") <- Embeddings(seurat_obj, "umap")
  
  return(sce)
}

stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
head(stratified[[]])
table(stratified$condition)
table(stratified$cell_type, stratified$condition)

DimPlot(stratified)
levels(Idents(stratified))

cluster.colors <- hue_pal()(10)
names(cluster.colors) <- levels(stratified$cell_type)
DimPlot(stratified, group.by = "cell_type", cols = cluster.colors)

sce <- prepare_data_for_trajectory(stratified, "condition", "cell_type")
sce <- slingshot(sce, 
                 clusterLabels = "cell_type", 
                 reducedDim = 'UMAP',
                 start.clus = "TZ", end.clus = NULL
)

par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(reducedDims(sce)$UMAP, main = paste0("Trajectories from TZ Root Cluster"), 
     pch=20, asp = 1, cex = 0.5,
     col = cluster.colors[colData(sce)$cell_type])
lines(SlingshotDataSet(sce), lwd=2, col= "black")

SlingshotDataSet(sce)
# Lineage1: TZ  AC basal distal  Wounded cells  Reactive cells  AC proliferative basal proximal  AC basal proximal  
# Lineage2: TZ  AC basal distal  Wounded cells  Reactive cells  AC suprabasal proximal  
# Lineage3: TZ  AC basal distal  Wounded cells  EMT-like cells  Cancerous cells  
# Lineage4: TZ  AC basal distal  AC suprabasal distal  

png(filename = paste0(results.dir, "trajectories/Pseudotime_Lineages.jpeg"), width = 1200, height = 400)
par(mfrow=c(1,4), mar=c(2,2,2,2))
plot(reducedDims(sce)$UMAP, 
     main = paste0("Lineage 1 from TZ Root Cluster"), 
     pch=20, asp = 1, cex = 0.5, col = cluster.colors[colData(sce)$cell_type])
lines(SlingshotDataSet(sce)@curves[[1]], lwd=2, col= "black")
plot(reducedDims(sce)$UMAP, 
     main = paste0("Lineage 2 from TZ Root Cluster"), 
     pch=20, asp = 1, cex = 0.5, col = cluster.colors[colData(sce)$cell_type])
lines(SlingshotDataSet(sce)@curves[[2]], lwd=2, col= "black")
plot(reducedDims(sce)$UMAP, 
     main = paste0("Lineage 3 from TZ Root Cluster"), 
     pch=20, asp = 1, cex = 0.5, col = cluster.colors[colData(sce)$cell_type])
lines(SlingshotDataSet(sce)@curves[[3]], lwd=2, col= "black")
plot(reducedDims(sce)$UMAP, 
     main = paste0("Lineage 4 from TZ Root Cluster"), 
     pch=20, asp = 1, cex = 0.5, col = cluster.colors[colData(sce)$cell_type])
lines(SlingshotDataSet(sce)@curves[[4]], lwd=2, col= "black")
dev.off()

png(filename = paste0(results.dir, "trajectories/Pseudotime_Lineages.2.jpeg"), width = 1200, height = 400)
par(mfrow=c(1,4), mar=c(2,2,2,2))
umap <- reducedDims(sce)$UMAP
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(nrow(umap))
lineages <- grep("slingPseudotime", colnames(colData(sce)), value = T)
for(j in 1:length(lineages)) {
  plotcol <- colors[cut(colData(sce)[, lineages[j]], breaks=nrow(umap))]
  plot(reducedDims(sce)$UMAP, 
       main = paste0("Lineage ", j, " from TZ Cluster"),
       col = plotcol, 
       pch=16, 
       asp = 1, 
       xlab = 'UMAP 1', 
       ylab = 'UMAP 2')
  lines(SlingshotDataSet(sce)@curves[[j]], lwd=2, col='black')
}
dev.off()

save(sce, file = paste0(data.dir, "slingshot_root.TZ.RData"))


## Condiments ----
# https://hectorrdb.github.io/condiments/index.html
# https://github.com/HectorRDB/condiments/
# https://www.nature.com/articles/s41467-024-44823-0#code-availability

load(file = paste0(data.dir, "slingshot_root.TZ.RData"))


#### focus on linege 3 ----
# that leads from TZ to cancerous cells

#### Trajectory Inference ----

scores <- imbalance_score(
  Object = reducedDims(sce)$UMAP, 
  condition = colData(sce)$condition)

df <- data.frame(
  umapcca_1 = reducedDims(sce)$UMAP[, 1],
  umapcca_2 = reducedDims(sce)$UMAP[, 2],
  condition = colData(sce)$condition,
  scores = scores$scores,
  scaled_scores = scores$scaled_scores,
  cl = colData(sce)$cell_type
)
head(df)

ggplot(df, aes(x = umapcca_1, y = umapcca_2, col = cl)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cluster.colors) +
  ggplot(df, aes(x = umapcca_1, y = umapcca_2, col = cl)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cluster.colors) +
  geom_path(data =  slingCurves(sce, as.df = TRUE) %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 0.5)
ggsave(paste0(results.dir, "Slingshot_trajectories_root4.jpeg"), width = 8, height = 4, dpi = 300)

# tested directly with the sce object (took 4 min)
topologyTest(sce, conditions = colData(sce)$condition, rep = 100, parallel = T)
# method thresh statistic p.value
# 1 Classifier   0.01  0.882609       0


#### Differential progression ----
## Kolmogorov-Smirnov Test (TradeSeq)
ks.test(slingPseudotime(sce)[colData(sce)$condition == "A.Normal", 3],
        slingPseudotime(sce)[colData(sce)$condition == "B.Hyperplasia", 3]) # p-value = 0.00236
ks.test(slingPseudotime(sce)[colData(sce)$condition == "B.Hyperplasia", 3],
        slingPseudotime(sce)[colData(sce)$condition == "C.Dysplasia", 3]) # p-value = 4.618e-10
ks.test(slingPseudotime(sce)[colData(sce)$condition == "C.Dysplasia", 3],
        slingPseudotime(sce)[colData(sce)$condition == "D.Carcinoma", 3]) # p-value < 2.2e-16

psts <- slingPseudotime(sce) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         condition = df$condition) %>%
  pivot_longer(starts_with("Lineage"), values_to = "pseudotime", names_to = "lineages")

ggplot(psts, aes(x = pseudotime, fill = condition)) +
  geom_density(alpha = .5) +
#  scale_fill_brewer(type = "qual") +
  facet_wrap(~lineages) +
  theme(legend.position = "bottom")
ggsave(paste0(results.dir, "trajectories/Pseudotime_density.png"), width = 7, height = 4)

# plot only lineage 3 (leading to epithelial cluster 10)
psts <- slingPseudotime(sce) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         condition = df$condition) %>%
  pivot_longer(starts_with("Lineage"), values_to = "pseudotime", names_to = "lineages")

psts.3 <- psts[psts$lineages=="Lineage3",]
psts.3 <- na.omit(psts.3)
ggplot(psts.3, aes(x = pseudotime, fill = condition)) +
  geom_density(alpha = .5) +
  #  scale_fill_manual(values = condition.colors) +
  #  scale_fill_brewer(type = "qual") +
  ggtitle("Differential progression\nLineage 3")
ggsave(paste0(results.dir, "trajectories/Pseudotime_density_Lineage.3.png"), width = 8, height = 4)


# progressionTest(condiments) is similar to KS test (Tradeseq)
# Test whether or not the pseudotime distribution are identical within lineages between conditions
prog_res <- progressionTest(sce, conditions = df$condition, global = TRUE, lineages = TRUE)
knitr::kable(prog_res)
#|lineage | statistic|   p.value|
#  |:-------|---------:|---------:|
#  |All     | 0.3066667| 0.0376242|
#  |1       | 0.3175862| 0.0566055|
#  |2       | 0.3337500| 0.0182185|
#  |3       | 0.4515385| 0.0000061|
#  |4       | 0.3337500| 0.0253150|


#### Differential differentiation ----

df$weight_1 <- slingCurveWeights(sce, as.probs = TRUE)[, 1]
df$weight_2 <- slingCurveWeights(sce, as.probs = TRUE)[, 2]
df$weight_3 <- slingCurveWeights(sce, as.probs = TRUE)[, 3]
df$weight_4 <- slingCurveWeights(sce, as.probs = TRUE)[, 4]
head(df)
gg1 <- ggplot(df, aes(x = weight_1, fill = condition)) +
  geom_density(alpha = .5) +
#  scale_fill_brewer(type = "qual") +
  labs(x = "Curve weight for lineage 1")
gg2 <- ggplot(df, aes(x = weight_2, fill = condition)) +
  geom_density(alpha = .5) +
#  scale_fill_brewer(type = "qual") +
  labs(x = "Curve weight for lineage 2")
gg3 <- ggplot(df, aes(x = weight_3, fill = condition)) +
  geom_density(alpha = .5) +
#  scale_fill_brewer(type = "qual") +
  labs(x = "Curve weight for lineage 3")
gg4 <- ggplot(df, aes(x = weight_4, fill = condition)) +
  geom_density(alpha = .5) +
#  scale_fill_brewer(type = "qual") +
  labs(x = "Curve weight for lineage 4")
ga1 <- grid.arrange(gg1, gg2, gg3, gg4, ncol = 2)
ggsave(paste0(results.dir, "trajectories/Curve_weight_density.png"), ga1, width = 8, height = 4)

# Test whether or not the cell repartition between lineages is independent of the conditions
# default, the classifier test (Lopez-Paz and Oquab 2016).
dif_res <- differentiationTest(sce, condition = df$condition, global = T, pairwise = TRUE)
knitr::kable(dif_res)
#|pair | statistic|   p.value|
#  |:----|---------:|---------:|
#  |All  | 0.3788889| 0.0000535|
#  |1vs2 | 0.3757143| 0.0005193|
#  |1vs3 | 0.3983333| 0.0002141|
#  |1vs4 | 0.3057895| 0.0580569|
#  |2vs3 | 0.4142424| 0.0000203|
#  |2vs4 | 0.3619512| 0.0006536|
#  |3vs4 | 0.2828571| 0.1888455|


#### Differential expression along lineages ----
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#differential-expression-1

#load(file = paste0(data.dir, "slingshot_root4.RData"))

# select number of knots
head(slingPseudotime(sce))
head(slingCurves(sce))
head(slingCurveWeights(sce))

icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
                   pseudotime = slingPseudotime(sce),
                   cellWeights = slingCurveWeights(sce),
                   conditions = factor(colData(sce)$condition),
                   nGenes = 100,
                   k = 3:7,
                   parallel = T, BPPARAM = BiocParallel::MulticoreParam(12)
)

# select most variable genes
#obj <- readRDS(file = "data/synovial/objects/mesenchymal_tumor.rds")
var_genes <- VariableFeatures(FindVariableFeatures(stratified, nfeatures = 1000)) # 
var_genes <- var_genes[!startsWith(var_genes, "mt-")] # 994
#rm(obj)

# this may take hours:
sce <- fitGAM(sce, 
              conditions = factor(colData(sce)$condition),
              genes = var_genes,
              parallel = T, BPPARAM = BiocParallel::MulticoreParam(12),
              nknots = 6) # based on the knot plots
mean(rowData(sce)$tradeSeq$converged)

save(sce, file = paste0(data.dir, "sce_fitGAM_root.TZ.RData"))


# 

load(file = paste0(data.dir, "sce_fitGAM_root.TZ.RData"))

# test whether average gene expression is associated with pseudotime.
rowData(sce)$assocRes <- associationTest(sce, global = T, lineages = TRUE, l2fc = log2(2))
head(na.omit(rowData(sce)$assocRes))
dim(na.omit(rowData(sce)$assocRes)) # 486 52

assocRes <- na.omit(rowData(sce)$assocRes)
head(assocRes)
table(assocRes$pvalue <= 0.05) # 235
hist(assocRes$meanLogFC)
table(assocRes$meanLogFC > 1) # 350
table(assocRes$pvalue <= 0.05 & assocRes$meanLogFC > 1) # 211

A.Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage3_conditionA.Normal, "fdr") <= 0.05)
] # 
B.Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage3_conditionB.Hyperplasia, "fdr") <= 0.05)
] # 
C.Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage3_conditionC.Dysplasia, "fdr") <= 0.05)
] # 
D.Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage3_conditionD.Carcinoma, "fdr") <= 0.05)
] # 
UpSetR::upset(fromList(list(Normal = A.Genes, 
                            Hyperplasia = B.Genes,
                            Dysplasia = C.Genes,
                            Carcinoma = D.Genes)),
              main.bar.color = "skyblue", sets.bar.color = "skyblue", 
              sets.x.label = "No. of DEGs", text.scale = 2)



write.csv(assocRes, file = paste0(results.dir, "trajectories/assocRes_root.TZ.csv"))

# DE genes in lineage 3
sel.genes <- c(A.Genes, B.Genes, C.Genes, D.Genes)
sel.genes <- unique(sel.genes) # 385

assocRes.lin3 <- assocRes[sel.genes, ]
colnames(assocRes.lin3)
assocRes.lin3 <- assocRes.lin3[, c(1:3,28:39,52)]
assocRes.lin3 <- assocRes.lin3[order(assocRes.lin3$pvalue_lineage3_conditionD.Carcinoma), ]
head(assocRes.lin3)
write.csv(assocRes.lin3, file = paste0(results.dir, "trajectories/assocRes_lineage3_root.TZ.csv"))


## Visualization of DE genes ----

plotSmoothers(sce, assays(sce)$counts, gene = "Dbi", alpha = 1, border = TRUE) + ggtitle("CDH1")

unique.D.Genes <- setdiff(D.Genes, c(A.Genes, B.Genes, C.Genes))
write.csv(unique.D.Genes, file = paste0(results.dir, "trajectories/unique.D.Genes_lineage3_root.TZ.csv"))
unique.C.D.Genes <- setdiff(c(C.Genes, D.Genes), c(A.Genes, B.Genes))
write.csv(unique.C.D.Genes, file = paste0(results.dir, "trajectories/unique.C.D.Genes_lineage3_root.TZ.csv"))

sel.genes <-  rownames(assocRes)[
  assocRes$pvalue_lineage3_conditionD.Carcinoma <= 0.01 & 
    assocRes$meanLogFC > 2
] # 

sel.genes <- rownames(assocRes.lin3)[1:50]


# due to reported inconsistencies
# checked this https://github.com/statOmics/tradeSeq/issues/230
# and this: https://github.com/statOmics/tradeSeq/issues/232
# and this:  https://github.com/statOmics/tradeSeq/issues/237
# I chose to use tidy=T in the heatmap

nPoints <- 100
#yhatSmooth <- predictSmooth(sce, gene = sel.genes, nPoints = nPoints, tidy = TRUE) # used for labeled heatmap
yhatSmooth <- predictSmooth(sce, gene = D.Genes, nPoints = nPoints, tidy = TRUE)
# Examine the structure of tidy data
head(yhatSmooth)
str(yhatSmooth)
# Filter for lineage 3 only
yhatSmooth_lineage3 <- yhatSmooth[yhatSmooth$lineage == "3", ]
# Filter for carcinoma only
yhatSmooth_lineage3 <- yhatSmooth_lineage3[yhatSmooth_lineage3$condition == "D.Carcinoma", ]
head(yhatSmooth_lineage3)
# Convert tidy data back to matrix format
# Assuming columns are: gene, pseudotime, yhat, lineage, conditions
library(tidyr)
library(tibble)
yhatMatrix <- yhatSmooth_lineage3 %>%
  dplyr::select(gene, time, yhat, condition) %>%  # adjust column names as needed
  unite("condition_time", condition, time, sep = "_") %>%
  pivot_wider(names_from = condition_time, values_from = yhat) %>%
  column_to_rownames("gene") %>%
  as.matrix()
head(yhatMatrix)
# Apply scaling
#yhatSmoothScaled <- t(log1p(t(yhatMatrix)))
yhatSmoothScaled <- t(scale(t(yhatMatrix)))
# Create heatmaps
heatSmooth <- pheatmap(yhatSmoothScaled,
                       cluster_cols = FALSE,
                       show_rownames = F, 
                       show_colnames = FALSE)

# plot only carcinoma curve and cells assigned to lineage 3
curvesCols <- c(rep("transparent",11), "#440154FF", rep("transparent",4))
sel.smooth.gene <- "Pdgfa"
plotSmoothers(sce, assays(sce)$counts, gene = sel.smooth.gene, curvesCols = curvesCols,
              border = FALSE) +
  ggtitle(sel.smooth.gene) +
  theme(legend.position = "none") +
  ggplot2::scale_color_manual(values = curvesCols)

# IL17 activity GO term genes present in stratified dataset
# "Socs3"    "Cxcl1"    "Cxcl10"   "Il1b"     "Il6"      "Nfkb1"    "Stat3"    "Nfkbiz"   "Traf3ip2" "Srsf1"    "Dcp1b" 

sel.smooth.gene <- "Nfkbiz"
plotSmoothers(sce, assays(sce)$counts, gene = sel.smooth.gene, curvesCols = curvesCols,
              border = FALSE) +
  ggtitle(sel.smooth.gene) +
  theme(legend.position = "none") +
  ggplot2::scale_color_manual(values = curvesCols)



# Extract the hierarchical clustering object for rows
row_clusters <- heatSmooth$tree_row

# You can visualize the dendrogram if you want
# plot(row_clusters)
# Cut the dendrogram into 3 clusters
# Adjust 'k' based on the number of main clusters you observe in your heatmap
num_clusters <- 3 # Change this value based on your observation
gene_clusters <- cutree(row_clusters, k = num_clusters)
# 'gene_clusters' is a named vector where names are gene names and values are cluster assignments
head(gene_clusters)
table(gene_clusters)
write.csv(gene_clusters, file = paste0(results.dir, "trajectories/gene_clusters_lineage3_root.TZ.csv"))


## pathways ----

gene_clusters <- read.csv(file = paste0(results.dir, 
                                        "trajectories/gene_clusters_lineage3_root.TZ.csv"),
                          row.names = 1)
head(gene_clusters)

library(enrichR)
dbs <- listEnrichrDbs()
head(dbs)
grep("Reactome", dbs$libraryName, value = T)
#dbs <- "MSigDB_Hallmark_2020"
#dbs <- "GO_Biological_Process_2023"
dbs <- "Reactome_Pathways_2024"

select.cluster <- 1
enriched <- enrichr(rownames(gene_clusters)[gene_clusters$x == select.cluster], dbs)
head(enriched[[1]], 20)
table(enriched[[1]]["Adjusted.P.value"] < 0.05) #

#write.csv(enriched[[1]], file = paste0(results.dir, "enrichR_lineage_3_cluster_", select.cluster, "_root.TZ.csv"))
write.csv(enriched[[1]], file = paste0(results.dir, "Reactome_lineage_3_cluster_", select.cluster, "_root.TZ.csv"))

# remove redundant terms
enriched[[1]] <- enriched[[1]][!grepl("Developmental Biology", enriched[[1]]$Term), ]
enriched[[1]] <- enriched[[1]][!grepl("Immune System", enriched[[1]]$Term), ]

b1 <- plotEnrich(enriched[[1]], 
           showTerms = 10, 
           numChar = 40, 
           y = "Count", 
           orderBy = "P.value",
           title = paste0("Lineage 3 Cluster ", 
                          select.cluster, 
           #               "\nGO Biological Process")
           "\nReactome Pathways")
)
#ggsave(paste0(results.dir, "enrichR_barplot_cluster_", select.cluster, ".png"), 
#       width = 8, height = 6, dpi = 300)
ggsave(paste0(results.dir, "Reactome_barplot_cluster_", select.cluster, ".png"),
       b1,
       width = 8, height = 6, dpi = 300)

# make a dotplot
df <- enriched[[1]]
head(df)
df$GeneRatio <- sapply(df$Overlap, function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1] / parts[2]  # Divides numerator by denominator
})
df$GeneCount <- sapply(df$Overlap, function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1]  # take only first number
})
d1 <- ggplot2::ggplot(df[1:10, ], aes(x = reorder(Term, -P.value), y = -Adjusted.P.value)) +
  geom_point(aes(size = GeneRatio, color = -Adjusted.P.value)) +
  scale_size(range = c(2, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  ylim(-0.085, 0.005) +
  coord_flip() +
  labs(x = "Term", y = "Adjusted.P.value") +
  ggtitle(paste0("Lineage 3 Cluster ", select.cluster, "")) +
  theme(
    axis.text.y = element_text(size = 10),       # Term text
    axis.text.x = element_text(size = 10),       # Y-axis values
    axis.title = element_text(size = 12),        # Axis titles
    plot.title = element_text(size = 14),        # Plot title
    legend.title = element_text(size = 12),      # Legend titles
    legend.text = element_text(size = 10)        # Legend text
  )
#ggsave(paste0(results.dir, "enrichR_dotplot_cluster_", select.cluster, ".png"), 
#       width = 12, height = 6, dpi = 300)
ggsave(paste0(results.dir, "Reactome_dotplot_cluster_", select.cluster, ".png"), 
       d1,
       width = 8, height = 6, dpi = 300)



select.cluster <- 2
enriched <- enrichr(rownames(gene_clusters)[gene_clusters$x == select.cluster], dbs)
head(enriched[[1]], 20)
table(enriched[[1]]["Adjusted.P.value"] < 0.05) #

#write.csv(enriched[[1]], file = paste0(results.dir, "enrichR_lineage_3_cluster_", select.cluster, "_root.TZ.csv"))
write.csv(enriched[[1]], file = paste0(results.dir, "Reactome_lineage_3_cluster_", select.cluster, "_root.TZ.csv"))

# remove redundant terms
enriched[[1]] <- enriched[[1]][!grepl("Attenuation Phase", enriched[[1]]$Term), ]

b2 <- plotEnrich(enriched[[1]], 
                 showTerms = 10, 
                 numChar = 40, 
                 y = "Count", 
                 orderBy = "P.value",
                 title = paste0("Lineage 3 Cluster ", 
                                select.cluster, 
                                #               "\nGO Biological Process")
                                "\nReactome Pathways")
)
#ggsave(paste0(results.dir, "enrichR_barplot_cluster_", select.cluster, ".png"), 
#       width = 8, height = 6, dpi = 300)
ggsave(paste0(results.dir, "Reactome_barplot_cluster_", select.cluster, ".png"),
       b2,
       width = 8, height = 6, dpi = 300)

# make a dotplot
df <- enriched[[1]]
head(df)
df$GeneRatio <- sapply(df$Overlap, function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1] / parts[2]  # Divides numerator by denominator
})
df$GeneCount <- sapply(df$Overlap, function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1]  # take only first number
})
d2 <- ggplot2::ggplot(df[1:10, ], aes(x = reorder(Term, -P.value), y = -Adjusted.P.value)) +
  geom_point(aes(size = GeneRatio, color = -Adjusted.P.value)) +
  scale_size(range = c(2, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  ylim(-0.0025, 0.0005) +
  coord_flip() +
  labs(x = "Term", y = "Adjusted.P.value") +
  ggtitle(paste0("Lineage 3 Cluster ", select.cluster, "")) +
  theme(
    axis.text.y = element_text(size = 10),       # Term text
    axis.text.x = element_text(size = 10),       # Y-axis values
    axis.title = element_text(size = 12),        # Axis titles
    plot.title = element_text(size = 14),        # Plot title
    legend.title = element_text(size = 12),      # Legend titles
    legend.text = element_text(size = 10)        # Legend text
  )
#ggsave(paste0(results.dir, "enrichR_dotplot_cluster_", select.cluster, ".png"), 
#       width = 12, height = 6, dpi = 300)
ggsave(paste0(results.dir, "Reactome_dotplot_cluster_", select.cluster, ".png"), 
       d2,
       width = 8, height = 6, dpi = 300)


b1+b2
ggsave(paste0(results.dir, "Reactome_barplot_clusters_1_2.png"), 
       width = 14, height = 6, dpi = 300)

d1+d2
ggsave(paste0(results.dir, "Reactome_dotplot_clusters_1_2.png"), 
       width = 16, height = 6, dpi = 300)


## scores ----

# plot pathways scores along lineage 3

stratified <- readRDS(file = paste0(data.dir, "stratified_manuscript.rds"))
head(stratified)

#load(file = paste0(data.dir, "slingshot_root4.RData"))
SlingshotDataSet(sce)
head(sce)
colData(sce)

all(rownames(colData(sce))==colnames(stratified))

stratified$slingPseudotime_3 <- colData(sce)$slingPseudotime_3
head(stratified[[]])

FeaturePlot(stratified, "slingPseudotime_3") +
  ggtitle("Pseudotime on Lineage 3, Root in TZ cells")
ggsave(filename = paste0(results.dir, "trajectories/slingPseudotime_3.png"), width = 6, height = 6)

meta <- stratified@meta.data
head(meta)
# remove NAs in lineage 3
meta <- meta[!is.na(meta$slingPseudotime_3), ] # 7057
table(meta$condition)
# A.Normal B.Hyperplasia   C.Dysplasia   D.Carcinoma 
# 88           121          1264          1503 
meta <- meta[order(meta$slingPseudotime_3, decreasing = F), ]

obj2 <- subset(stratified, cells = rownames(meta))
VlnPlot(obj2, "GO_00973981", group.by = "sample") +
  geom_boxplot(width = 0.2, outlier.size = 0.5, outlier.colour = "black", alpha = 0.8) +
  ggtitle("IL17 Score on Lineage 3, Root in TZ cells") +
  xlab("") +
  ylim(-0.5, 1.5) +
#  scale_fill_manual(values = as.vector(alphabet())[c(8,2)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "wilcox.test", 
                   comparisons = list(c("Normal", "Hyperplasia"),
                                      c("Hyperplasia", "Dysplasia"),
                                      c("Dysplasia", "Carcinoma")),
                     label = "p.format", size = 5, label.x.npc = "center")
ggsave(filename = paste0(results.dir, "trajectories/GO_00973981_VlnPlot.png"), width = 7, height = 6)

head(stratified)
ggplot(meta, aes(x = slingPseudotime_3, y = GO_00973981)) +
  geom_point(alpha = 0.4, size = 0.2) + # Optional: make the raw points slightly transparent
  geom_line(alpha = 0.4, color = "white") + # Optional: make the raw lines slightly transparent
  geom_smooth(method = "loess", se = F, size = 2, alpha = 1) + # Adjust size for visibility
  labs(title = "IL17 Activity Score Over Pseudotime",
       subtitle = "Lineage 3, Root in TZ cells",
       x = "Pseudotime",
       y = "IL17 Score [GO:0097398]",
       color = "Condition") + # Customize the legend title
  theme_classic()
ggsave(filename = paste0(results.dir, "trajectories/GO_00973981_Pseudotime.png"), width = 5, height = 5)

meta_long <- meta %>%
  pivot_longer(cols = c(GO_00973981, GO_00050311, GO_00701021, GO_00050241),
               names_to = "GO_term", 
               values_to = "value")

ggplot(meta_long, aes(x = slingPseudotime_3, y = value, color = GO_term)) +
  geom_point(alpha = 0.4, size = 0.2) +
  geom_smooth(method = "loess", se = F, size = 2, alpha = 1) +
       xlab("Pseudotime") +
       ylab("Pathway Score") +
  scale_color_manual(values = c("GO_00973981" = "blue", 
                                "GO_00050311" = "red",
                                "GO_00701021" = "green", 
                                "GO_00050241" = "orange"),
                     labels = c("GO_00973981" = "IL17 activity",
                                "GO_00050311" = "TNF activity", 
                                "GO_00701021" = "IL6 activity",
                                "GO_00050241" = "TGFb activity"),
                     name = "GO Terms") +
  ggtitle("Pathway Activity Scores Over Pseudotime",
         subtitle = "Lineage 3: from TZ cells to cancerous cells") +
  theme_classic()
ggsave(filename = paste0(results.dir, "trajectories/Pathway_Activity_Scores_Pseudotime.png"), width = 6, height = 5)

# GO_0097398: cellular response to interleukin-17
# GO_0005031: tumor necrosis factor receptor activity
# GO_0070102: interleukin-6-mediated signaling pathway
# GO_0005024: transforming growth factor beta receptor activity




  


# end ---------------------------------------------------------------------
sessionInfo()
