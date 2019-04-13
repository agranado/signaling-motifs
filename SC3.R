#Analysis of Tabula Muris tissues using SC3


setwd("MEGA/Caltech/rnaseq/signaling-clusters/")
load("../tabula-muris/tiss_filteredJan10_2019.rdata")
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


###

tissue = "Large_Intestine"
meta.data %>% group_by(tissue) %>% summarize(count= n())
intes.index = meta.data$tissue=="Large_Intestine"


intes.norm  = norm.bmp[,intes.index]
intes.counts  = counts.bmp[,intes.index]
intes.meta = meta.data[intes.index,]


#correction
intes.counts = tiss@raw.data[bmp.genes,colnames(intes.norm)]
library(scater)
library(SC3)



#create SC3 object: you the both the normalized and raw counts
sce<-SingleCellExperiment(assays = list( counts = as.matrix(intes.counts),logcounts = as.matrix(intes.norm)),colData = intes.meta)

rowData(sce)$feature_symbol = rownames(intes.counts)
rownames(sce)<-rownames(intes.counts)

tiss<- Convert(from=sce,to="seurat")
colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads' # this is not UMI data so Seurat calculates only the number of reads

tiss <- NormalizeData(object = tiss, scale.factor = 1e6) #default normalization by Seurat
#rewrite the values for counts and nGenes according to the BMP
tiss@meta.data$nGene = Matrix::colSums(intes.counts>0)
names(tiss@meta.data)[3]<-"_nReads"


#convert to Seurat:
tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5,num.bin =10)
n.pcs = 6
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 10)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)
#test several resolution parameters:

tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
                     resolution = c(0.5,0.7,0.9,1.1,1.5,2,3), print.output = 0, save.SNN = TRUE,force.recalc=T) #DONE

#perplexity means how many close nieghbours each cell has
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs, seed.use = 10, perplexity=20,
                                     check_duplicates = F)

#analyse how clusters change with the resolution parameter 
# FROM: http://oshlacklab.com/combes-organoid-paper/04_Organoids_Clustering.html
x11();clustree(tiss)
 # # # # # # #
#UMAP: https://satijalab.org/seurat/conversion_vignette.html
p2 <- DimPlot(object = tiss, reduction.use = "umap", no.legend = TRUE, do.return = TRUE,
                vector.friendly = F, pt.size = 0.5) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
p2
#make it big
p2 <- DimPlot(object = tiss, reduction.use = "umap", no.legend = TRUE, do.return = TRUE,
              vector.friendly = F, pt.size = 3) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

p1<-RidgePlot(tiss,features.plot = bmp.genes, y.lab.rot=T,group.by ="ident",do.return=T)
#or
p1<-RidgePlot(tiss,features.plot = tiss@var.genes, y.lab.rot=T,group.by ="ident",do.return=T)
