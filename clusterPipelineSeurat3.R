#latest method (after loading the seurat.pathway object )

setwd("/home/agranado/MEGA/Caltech/rnaseq/signaling-clusters/")

#optimize for Seurat 3
load("../tabula-muris/tiss_filteredJan10_2019.rdata")
source("bmp_profiles.R")

library(Matrix)
library(SC3)
library(scater)
library(Seurat)

get.pathway.expression<-function( pathway.genes_,  min.genes.per.cell = 3,frac.express.components = 2){
  #let's take the raw data to manual normalization then SC3
  #manual.norm<-tiss@raw.data
  manual.norm<-tiss[["RNA"]]@counts
  #take pathway genes to filter cells:
  manual.norm.pathway<-manual.norm[pathway.genes_,]
  #at least twice the number of genes for nReads
  count.cells = Matrix::colSums(manual.norm.pathway)
  manual.norm.pathway = manual.norm.pathway[, count.cells>length(pathway.genes_)* frac.express.components]
  #at least 3 genes have >0 expression values
  genes.at.least.one = Matrix::colSums(manual.norm.pathway>0)
  manual.norm.pathway = manual.norm.pathway[, genes.at.least.one>= min.genes.per.cell]
  #NOW: let's take the cells in the new filtered matrix from the normalized data in seurat.pathway@data
  #


  #let's find which cells overlap with the seurat.pathway@data (we did take raw.data after all...)
  manual.norm.pathway=manual.norm.pathway[,which(colnames(manual.norm.pathway) %in% colnames(tiss[["RNA"]]@data))]


  #meta data:
  meta.data = tiss@meta.data
  meta.data  = meta.data[colnames(manual.norm.pathway),]


  counts.pathway<-manual.norm.pathway #this was always the raw data
  norm.pathway<-tiss[["RNA"]]@data [pathway.genes_,colnames(manual.norm.pathway)]
  dim(counts.pathway)
  rm(manual.norm)
  rm(manual.norm.pathway)

  return( list(counts.pathway, norm.pathway, meta.data ))
}

setup.seurat<-function(counts.pathway, norm.pathway, meta.data){
    #sce<-SingleCellExperiment(assays = list( counts = as.matrix(counts.pathway),logcounts = as.matrix(norm.pathway)),colData = meta.data)
    #seurat.pathway<-Convert(from=sce,to="seurat")
    seurat.pathway = CreateSeuratObject(counts = counts.pathway, project = "motifs",meta.data = meta.data)
    colnames(seurat.pathway@meta.data)[colnames(seurat.pathway@meta.data) == 'nUMI'] <- 'nReads' # this is not UMI data so Seurat calculates only the number of reads

    #make Seurat think we normalized using their method...
    seurat.pathway <- NormalizeData(object = seurat.pathway, scale.factor = 1e6) #default normalization by Seurat

    nGene = Matrix::colSums(counts.pathway>0)

    seurat.pathway@meta.data$nGene = nGene

    nReads = Matrix::colSums(counts.pathway)

    names(seurat.pathway@meta.data)[3] = "_nReads"

    seurat.pathway[["RNA"]]@data = norm.pathway
    return(seurat.pathway)
}


all.pathways = c("glucose","notch","shh","inositol","mtorc")
pca.manual = c(9,11,14,13,8)
for( p in 1:length(all.pathways)){

    which.pathway = all.pathways[p]
    pathway.genes_ = pathway.genes(which.pathway)
    #check for existing genes in the tiss object before retrieving data
    pathway.genes_ = pathway.genes_[pathway.genes_ %in% rownames(tiss[["RNA"]]@counts)]
    #function call
    res.list = get.pathway.expression( pathway.genes_)
    counts.pathway = res.list[[1]]
    norm.pathway = res.list[[2]]
    meta.data = res.list[[3]]

    #function call
    seurat.pathway = setup.seurat( counts.pathway, norm.pathway, meta.data) #works v3

    #######
     seurat.pathway <- ScaleData(object = seurat.pathway, features = rownames(seurat.pathway))
     pcs.compute = round(length(pathway.genes_) * 3/4)
     seurat.pathway <- RunPCA(object = seurat.pathway, features =pathway.genes_, do.print = FALSE, npcs = pcs.compute,maxit =10000) #no print

    # seurat.pathway <- ProjectPCA(object = seurat.pathway, do.print = FALSE) #no print

     a= ElbowPlot(object = seurat.pathway, ndims = pcs.compute)


     pca.used = pca.manual[p]

     seurat.pathway = JackStraw(seurat.pathway, num.replicate = 200, prop.freq = .1)
     seurat.pathway = ScoreJackStraw(seurat.pathway, dims = 1:10)


     filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_JackStraw_", toString(pcs.compute), "_.pdf",sep="")
     pdf(filename)
      JackStrawPlot(seurat.pathway, dims = 1:10)
     dev.off()


     #plot call
     filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_elbow_", toString(pcs.compute), "_.pdf",sep="")
     pdf(filename)
       plot(a)
     dev.off()


    # min.pcn = 4
    # diff(a$data$data.use[-c(1:min.pcn)])


    res.used = 0.9
    seurat.pathway = FindNeighbors(seurat.pathway, dims = 1:10) #for Seurat3

  #  seurat.pathway <- FindClusters(object = seurat.pathway, reduction.type = "pca", dims.use = 1:pca.used,
  #                        resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc=T ) #, plot.SNN =T) #DONE

    seurat.pathway = FindClusters(seurat.pathway, resolution = 0.5)
    #FindClusters now automatically sets the ident (and erases other idents)

  #  seurat.pathway<-SetIdent(seurat.pathway,ident=seurat.pathway@meta.data$res.0.9) #just in case
    #V2
    #seurat.pathway <- RunTSNE(object = seurat.pathway, dims.use = 1:pca.used, seed.use = 10, perplexity=30,
    #                                       check_duplicates = F)
    #v3:
    seurat.pathway<-RunTSNE(seurat.pathway, dims = 1:10, perplexity = 30)
    seurat.pathway = RunUMAP(seurat.pathway, dims = 1:10)


    #plot call :
    filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_TSNE_", toString(pcs.compute), "_.pdf",sep="")
    pdf(filename)
    # TSNEPlot(seurat.pathway)
    DimPlot(seurat.pathway,reduction="tsne")
    dev.off()


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


        filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_clustering_",cell.vars[i] ,"_", quant.vars[j],"_.pdf",sep="")


                      heatmap.pipeline2.v3(seurat.obj  = seurat.pathway, which.path = which.pathway, which.var =cell.vars[i],quant.var = quant.vars[j],filename = filename,
                      cluster.cols =T,fontsize = 5,width = width,height=height)
      }
    }

    save(seurat.pathway,file = paste("../datasets/TabulaMuris_bmp/", which.pathway, "_clusteredOK_NoVarGenes_04082019.rda"))
    rm(seurat.pathway)


}
