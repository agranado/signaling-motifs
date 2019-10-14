


suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scMerge)
  library(scater)
  })

root.folder = "/home/agranado/"
# needed for SCT integration to run (otherwise it will crash )
options(future.globals.maxSize = 4000 * 1024^2)


# all cells (4 technologies/papers)
# meta data includes a tech field that we can use to separate them
# Also available at Google Drive CaltechProject/rna-seq/2019_postLM/RawData/Pancreas_ButlerCell2019
pancreas.data <- readRDS(file = paste(root.folder,"Downloads/pancreas_v3_files/pancreas_expression_matrix.rds",sep=""))
metadata <- readRDS(file = paste(root.folder,"Downloads/pancreas_v3_files/pancreas_metadata.rds",sep=""))

sce_pancreas = SingleCellExperiment(list(counts = as.matrix(pancreas)),colData = metadata)

# We need a logcounts assay
# we can compute it directly OR we can normalize with SCRAN (generates logcount assay automatically )
clusters <- quickCluster(sce_pancreas, min.size=100) # TAKES a long time
sce <- computeSumFactors(sce_pancreas, cluster=clusters)

#with tpm function from Fxns.R
log_counts = log2(1+tpm(as.matrix(pancreas.data)))
# OR
sce_pancreas<-normalize(sce_pancreas) #built in normalization for SCE
