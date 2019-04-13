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


#### additional filter:
# based on what we know, we can define which cells have active BMP signaling:

# Heidi:
# I think you need (at a minimum):
# * A Type I receptor
    type1.receptor=  c("Acvrl1","Acvr1","Bmpr1a","Acvr1b","Tgfbr1","Bmpr1b","Acvr1c")
# * A Type II receptor
    type2.receptor = c("Acvr2a","Acvr2b","Bmpr2","Tgfbr2")
# * A BMP R-Smad (so, Smad1 or 5 or 8)
    r_smads =  c("Smad1","Smad5","Smad8")
# * Smad4 (which complexes with R-Smad)
#     c("Smad4")

fil1 =  Matrix::colSums(counts.bmp[ rownames(counts.bmp) %in% type1.receptor ,])>0
fil2 =  Matrix::colSums(counts.bmp[ rownames(counts.bmp) %in% type2.receptor ,])>0
fil3 =  Matrix::colSums(counts.bmp[ rownames(counts.bmp) %in% r_smads ,])>0
fil4 =  sum(counts.bmp["Smad4",]>0)

bmp.active.cells = fil1 & fil2 & fil3 & fil4

counts.bmp =counts.bmp[,bmp.active.cells]
norm.bmp = norm.bmp[, bmp.active.cells]
meta.data = meta.data[bmp.active.cells,]

library(SC3)
library(scater)

sce<-SingleCellExperiment(assays = list( counts = as.matrix(counts.bmp),logcounts = as.matrix(norm.bmp)),colData = meta.data)
tiss<-Convert(from=sce,to="seurat")

colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads' # this is not UMI data so Seurat calculates only the number of reads

#make Seurat think we normalized using their method...
tiss <- NormalizeData(object = tiss, scale.factor = 1e6) #default normalization by Seurat

nGene = Matrix::colSums(counts.bmp>0)

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
tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:15,
                      resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc=T ) #, plot.SNN =T) #DONE

tiss<-SetIdent(tiss,ident=tiss@meta.data$res.1.2) #just in case

tiss <- RunTSNE(object = tiss, dims.use = 1:18, seed.use = 10, perplexity=30,
                                       check_duplicates = F)

x11();TSNEPlot(tiss)


#save heatmaps
cell.vars=c("ontology","tissue","id")
quant.vars =c("pct.exp","avg.log.exp","avg.exp")
for(i in 1:length(cell.vars)){
  for(j in 1:length(quant.vars)){

    if(cell.vars[i]=="ontology"){
      width = 10;height=15
    }else{
      width = 6; height =10
    }


    filename = paste("plots/bmp_PCA_allCells/heatmaps/clustering_",cell.vars[i] ,"_", quant.vars[j],"_.pdf",sep="")


                  heatmap.pipeline2(which.var =cell.vars[i],quant.var = quant.vars[j],filename = filename,
                  cluster.cols =F,fontsize = 5,width = width,height=height)
  }
}

save(tiss,file = "../datasets/TabulaMuris_bmp/bmp_clusteredOK_NoVarGenes_03142019.rda")
