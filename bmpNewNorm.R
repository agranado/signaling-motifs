#latest method (after loading the tiss object )

setwd("MEGA/Caltech/rnaseq/signaling-clusters/")
load("../tabula-muris/tiss_filteredJan10_2019.rdata")
source("bmp_profiles.R")
bmp.genes = pathway.genes("bmp")

library(Matrix)
#let's take the raw data to manual normalization then SC3
manual.norm<-tiss@raw.data
#take bmp genes to filter cells:
manual.norm.bmp<-manual.norm[bmp.genes,]
#at least twice the number of genes for nReads
count.cells = Matrix::colSums(manual.norm.bmp)
manual.norm.bmp = manual.norm.bmp[, count.cells>length(bmp.genes)*2]
#at least 3 genes have >0 expression values
genes.at.least.one = Matrix::colSums(manual.norm.bmp>0)
manual.norm.bmp = manual.norm.bmp[, genes.at.least.one>=3]
#NOW: let's take the cells in the new filtered matrix from the normalized data in tiss@data
#


#let's find which cells overlap with the tiss@data (we did take raw.data after all...)
manual.norm.bmp=manual.norm.bmp[,which(colnames(manual.norm.bmp) %in% colnames(tiss@data))]


#meta data:
meta.data = tiss@meta.data
meta.data  = meta.data[colnames(manual.norm.bmp),]


counts.bmp<-manual.norm.bmp #this was always the raw data
norm.bmp<-tiss@data[bmp.genes,colnames(manual.norm.bmp)]
dim(counts.bmp)
rm(manual.norm)
rm(manual.norm.bmp)

library(SC3)
library(scater)

sce<-SingleCellExperiment(assays = list( counts = as.matrix(counts.bmp),logcounts = as.matrix(norm.bmp)),colData = meta.data)
tiss<-Convert(from=sce,to="seurat")

colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads' # this is not UMI data so Seurat calculates only the number of reads

#make Seurat think we normalized using their method...
tiss <- NormalizeData(object = tiss, scale.factor = 1e6) #default normalization by Seurat

nGenes = Matrix::colSums(counts.bmp>0)

tiss@meta.data$nGene = nGene

nReads = Matrix::colSums(counts.bmp)

names(tiss@meta.data)[3] = "_nReads"

tiss@data = norm.bmp
#######
tiss <- ScaleData(object = tiss)

 tiss <- RunPCA(object = tiss, pc.genes =bmp.genes, do.print = FALSE, pcs.compute = 20)

 tiss <- ProjectPCA(object = tiss, do.print = FALSE)

 PCElbowPlot(object = tiss, num.pc = 20)

res.used = 1.2
tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:18,
                      resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc=T) #DONE

tiss <- RunTSNE(object = tiss, dims.use = 1:18, seed.use = 10, perplexity=30,
                                       check_duplicates = F)

x11();TSNEPlot(tiss)



save(tiss,file = "../datasets/TabulaMuris_bmp/bmp_clusteredOK_NoVarGenes_03142019.rda")
