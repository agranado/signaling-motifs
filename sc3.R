#clustering of tabula muris data using SC3-R package (alternative method to Seurat)



library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)


setwd("/home/agranado/MEGA/Caltech/rnaseq/tabula-muris")
load("00_data_ingest/11_global_robj/FACS_all.Robj")
sce <- Convert(from=tiss_FACS, to="sce")

#maybe there is a bug in which Suerat conversion stores dgCMatrix objects into SCE which doesnt work
#this works : from https://github.com/hemberg-lab/SC3/issues/53
counts(sce) <- as.matrix(counts(sce))
#normcounts(sce) <- as.matrix(normcounts(sce))
logcounts(sce) <- as.matrix(logcounts(sce))
gc()
rowData(sce)$feature_symbol = rowData(sce)$gene
