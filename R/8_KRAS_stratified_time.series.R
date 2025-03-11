
# KRAS model - Time Series Trajectories on stratified cells

# summary -----------------------------------------------------------------

# https://github.com/BaderLab/Tempora
# https://www.baderlab.org/Software/Tempora

# browseVignettes(package = 'psupertime')
# https://github.com/wmacnair/psupertime
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9235474/


# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/TZones/")

suppressPackageStartupMessages({
  library(Seurat)
  library('psupertime')
  library('SingleCellExperiment')
})

set.seed(2)


# Trajectories ------------------------------------------------------------

## all stratified cells ----

filt.sct <- readRDS(file = "KRAS_model/objects/stratified.rds")
head(filt.sct[[]])
table(filt.sct$sample)
table(filt.sct$sample, filt.sct$labels_Chloe)
DimPlot(filt.sct)


# tested with SCT and integrated assays
# and done for all genes and TFs separately

DefaultAssay(filt.sct) <- "SCT"
sce <- as.SingleCellExperiment(filt.sct)
y           = sce$sample

psuper_obj  = psupertime(sce, y, sel_genes='hvg', penalization = "best")
psuper_obj  = psupertime(sce, y, sel_genes='tf_mouse', penalization = "best")

save(psuper_obj, file = "KRAS_model/objects/psupertime_stratified_SCT_hgv.RData")
save(psuper_obj, file = "KRAS_model/objects/psupertime_stratified_SCT_tf.RData")

DefaultAssay(filt.sct) <- "integrated"
sce <- as.SingleCellExperiment(filt.sct)
y           = sce$sample
psuper_obj  = psupertime(sce, y, sel_genes='hvg', penalization = "best")
psuper_obj  = psupertime(sce, y, sel_genes='tf_mouse', penalization = "best")

save(psuper_obj, file = "KRAS_model/objects/psupertime_stratified_integrated_hgv.RData")
save(psuper_obj, file = "KRAS_model/objects/psupertime_stratified_integrated_tf.RData")

test <-psuper_obj$beta_dt
head(test, 20)
write.csv(test, file = "KRAS_model/results/stratified/TimeSeries/SCT/hvg.correlations.csv")
write.csv(test, file = "KRAS_model/results/stratified/TimeSeries/SCT/tf.correlations.csv")
write.csv(test, file = "KRAS_model/results/stratified/TimeSeries/integrated/hvg.correlations.csv")
write.csv(test, file = "KRAS_model/results/stratified/TimeSeries/integrated/tf.correlations.csv")

# Model diagnostics
g       = plot_train_results(psuper_obj)
(g)
# psupertime ordering of cells
g       = plot_labels_over_psupertime(psuper_obj, label_name='sample')
(g)
# Genes identified by psupertime
g       = plot_identified_gene_coefficients(psuper_obj, n = 50)
(g)
g       = plot_identified_genes_over_psupertime(psuper_obj, n = 30, plot_ratio = 1.2, label_name='sample')
(g)
# psupertime as a classifier
g       = plot_predictions_against_classes(psuper_obj)
(g)


#

go_list 		= psupertime_go_analysis(psuper_obj, org_mapping='org.Mm.eg.db')

g 				= plot_profiles_of_gene_clusters(go_list, label_name='Donor age\n(years)')
(g)

# plot go terms
g 				= plot_go_results(go_list)
(g)

plot_heatmap_of_gene_clusters
plot_profiles_of_gene_clusters


## by cell type (Atlas) ----

rm(list=ls())

filt.sct <- readRDS(file = "KRAS_model/objects/stratified_Atlas.rds")
head(filt.sct[[]])
table(filt.sct$condition)
table(filt.sct$condition, filt.sct$cell_type)
DimPlot(filt.sct)

# re-labeling to ensure enough cells per group
Idents(filt.sct) <- filt.sct$cell_type
filt.sct <- RenameIdents(filt.sct, 
                         "AC basal distal" = "AC_basal",
                         "AC basal proximal" = "AC_basal",
                         "AC suprabasal distal" = "AC_basal",
                         "AC suprabasal proximal" = "AC_suprabasal",
                         "TZ basal" = "TZ",
                         "TZ suprabasal" = "TZ")
table(Idents(filt.sct))
# AC_basal AC_suprabasal            TZ 
# 4960          1534          1462
table(Idents(filt.sct), filt.sct$condition)
# Atlas Hyperplasia Dysplasia Carcinoma
# AC_basal       1901          73      1134      1852
# AC_suprabasal   449         161       459       465
# TZ              863          34       437       128

filt.sct$major.labels <- as.factor(Idents(filt.sct))
for(l in levels(filt.sct$major.labels)){
  
  sc.temp1 <- subset(filt.sct, idents = l)
  DefaultAssay(sc.temp1) <- "SCT"
  sc.temp2 <- subset(filt.sct, idents = l)
  DefaultAssay(sc.temp2) <- "integrated"
  
  sce1 <- as.SingleCellExperiment(sc.temp1)
  y1           = sce1$condition
  sce2 <- as.SingleCellExperiment(sc.temp2)
  y2           = sce2$condition
  
  psuper_obj1  = psupertime(sce1, y1, sel_genes='hvg', penalization = "best")
  psuper_obj2  = psupertime(sce1, y1, sel_genes='tf_mouse', penalization = "best")
  psuper_obj3  = psupertime(sce2, y2, sel_genes='hvg', penalization = "best")
  psuper_obj4  = psupertime(sce2, y2, sel_genes='tf_mouse', penalization = "best")
  
  super.list <- list(SCT_hgv = psuper_obj1,
                     SCT_tf = psuper_obj2,
                     INT_hgv = psuper_obj3,
                     INT_tf = psuper_obj4)
  
  for(i in 1:length(super.list)){
    
    super <- super.list[[i]]
    save(super, file = paste0("KRAS_model/objects/psupertime_stratified.Atlas_", l, "_", names(super.list)[i], ".RData"))
    write.csv(super$beta_dt, file = paste0("KRAS_model/results/stratified/TimeSeries/Atlas/", l, "_", names(super.list)[i], ".correlations.csv"))

    # Plots
    g       = plot_train_results(super)
    ggplot2::ggsave(filename = paste0("KRAS_model/results/stratified/TimeSeries/Atlas/", l, "_", names(super.list)[i], ".diagnostics.jpeg"), g, height=8, width=5)
    g       = plot_labels_over_psupertime(super, label_name='condition')
    ggplot2::ggsave(filename = paste0("KRAS_model/results/stratified/TimeSeries/Atlas/", l, "_", names(super.list)[i], ".cell_ordering.jpeg"), g, height=3, width=5)
    g       = plot_identified_gene_coefficients(super, n = 50)
    ggplot2::ggsave(filename = paste0("KRAS_model/results/stratified/TimeSeries/Atlas/", l, "_", names(super.list)[i], ".geneplot1.jpeg"), g, height=4, width=8)
    g       = plot_identified_genes_over_psupertime(super, n = 30, plot_ratio = 1.2, label_name='sample')
    ggplot2::ggsave(filename = paste0("KRAS_model/results/stratified/TimeSeries/Atlas/", l, "_", names(super.list)[i], ".geneplot2.jpeg"), g, height=12, width=14)
    g       = plot_predictions_against_classes(super)
    ggplot2::ggsave(filename = paste0("KRAS_model/results/stratified/TimeSeries/Atlas/", l, "_", names(super.list)[i], ".classifier.jpeg"), g, height=5, width=7)
    
  }
  
  rm(sc.temp1, sc.temp2, sce1, sce2, y1, y2, psuper_obj1, psuper_obj2, psuper_obj3, psuper_obj4, super.list, super, g)
  
}





# end ---------------------------------------------------------------------

