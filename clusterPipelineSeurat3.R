#latest method (after loading the seurat.pathway object )

setwd("/home/agranado/MEGA/Caltech/rnaseq/signaling-clusters/")

#optimize for Seurat 3
load("../tabula-muris/tiss_filteredJan10_2019.rdata")
tiss = UpdateSeuratObject(tiss)
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

all.pathways = list.files("../pathways/")
all.pathways = all.pathways[grep(".csv",all.pathways)]
all.pathways = gsub(patter = "\\.csv$","",all.pathways)


#main clustering, plotting , pca ALL pipeline
pca.manual = c(9,11,14,13,8)
plot.heatmaps = F
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
     pcs.compute = length(pathway.genes_) -1
     seurat.pathway <- RunPCA(object = seurat.pathway, features =pathway.genes_, do.print = FALSE, npcs = pcs.compute,maxit =10000) #no print

    # seurat.pathway <- ProjectPCA(object = seurat.pathway, do.print = FALSE) #no print



     pca.used = pca.manual[p]

     seurat.pathway = JackStraw(seurat.pathway, num.replicate = 200, prop.freq = .1)
     seurat.pathway = ScoreJackStraw(seurat.pathway, dims = 1:10)


     filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_JackStraw_", toString(pcs.compute), "_.pdf",sep="")
     pdf(filename)
      JackStrawPlot(seurat.pathway, dims = 1:10)
     dev.off()


     #plot call
     filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_elbow_", toString(pcs.compute), "_.pdf",sep="")
     a= ElbowPlot(object = seurat.pathway, ndims = pcs.compute)
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
    a = DimPlot(seurat.pathway,reduction="tsne")
    plot(a)
    dev.off()

    if(plot.heatmaps ==T){
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
    }

    save(seurat.pathway,file = paste("../datasets/TabulaMuris_bmp/", which.pathway, "_clusteredOK_NoVarGenes_04082019.rda"),sep="")
    rm(seurat.pathway)


}



make.plots=function(list.pathways, type ="tsne", file = "",seurat.pathway =NA){

  for(i in 1:length(list.pathways)){
    which.pathway = list.pathways[i]
    if(file==""){
      file = paste("../datasets/TabulaMuris_bmp/", which.pathway, "_clusteredOK_NoVarGenes_04082019.rda")
      load(file)
    }else{
      seurat.pathway = seurat.pathway
    }
      Sys.sleep(3)
    filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_TSNE_.pdf",sep="")
    pdf(filename)

        a=DimPlot(object =seurat.pathway,reduction ="tsne",do.return=T)
        plot(a)

        a = DimPlot(object =seurat.pathway,reduction ="umap",do.return = T)
        plot(a)

        a=DimPlot(object =seurat.pathway,reduction ="pca",do.return = T)
        plot(a)

    dev.off()
  }


}

make.random.pathways<-function(obj = tiss,path = "../pathways/random/", path.size = c(15,20,30,40,60,80,100), nrepeats = 100){

  all.genes = rownames(tiss[['RNA']]@data)
  for(j in 1:nrepeats){
    rep1 = lapply(path.size, sample, x=all.genes )
    for(i in 1:length(rep1)){
      write.csv(rep1[[i]], file =paste( path,"rand_",toString(length(rep1[[i]] )),"_",toString(runif(1)),"_.csv",sep="" ))
    }
  }
}


################################################
#loads one of the previously saved objects
do.pca<-function(which.pathway=""){
    if(exists("seurat.pathway")) rm(seurat.pathway)
    file = paste("../datasets/TabulaMuris_bmp/", which.pathway, "_clusteredOK_NoVarGenes_04082019.rda",sep="")
    load(file)
    pathway.genes_ = rownames(seurat.pathway)
    npcs = length(pathway.genes_)-1
    seurat.pathway <- RunPCA(object = seurat.pathway, features =pathway.genes_, do.print = FALSE, npcs = npcs,maxit =10000) #
    return(seurat.pathway@reductions$pca@stdev)
}
#takes the list and tiss object
do.pca.from.list<-function(which.pathway = ""){
  pathway.genes_ = pathway.genes(pathway = which.pathway)
  res.list = get.pathway.expression( pathway.genes_)
  counts.pathway = res.list[[1]]
  norm.pathway = res.list[[2]]
  meta.data = res.list[[3]]

  #function call
   seurat.pathway = setup.seurat( counts.pathway, norm.pathway, meta.data) #works v3


   seurat.pathway <- ScaleData(object = seurat.pathway, features = rownames(seurat.pathway))
   pcs.compute = length(pathway.genes_) -1
   seurat.pathway <- RunPCA(object = seurat.pathway, features =pathway.genes_, do.print = FALSE, npcs = pcs.compute,maxit =10000) #no print

   return(seruat.pathway )
}
#############################################
cl <- parallel::makeForkCluster(6)
doParallel::registerDoParallel(cl)

do.pca.all = function(pathway.list){

    results  = foreach(p = pathway.list) %dopar% do.pca(which.pathway =p)

}

    parallel::stopCluster(cl)
