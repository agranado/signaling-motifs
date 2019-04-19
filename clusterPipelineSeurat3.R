#latest method (after loading the seurat.pathway object )

#setwd("/home/agranado/MEGA/Caltech/rnaseq/signaling-clusters/")

#optimize for Seurat 3
load("../tabula-muris/tiss_filteredJan10_2019.rdata")

library(Matrix)
library(SC3)
library(scater)
library(Seurat)
library(data.table)

tiss = UpdateSeuratObject(tiss)


#get the annotated genes (manually curated lists)
pathway.genes<-function(pathway ="bmp"){
  bmp.receptors<-c("Bmpr1a","Bmpr1b","Acvr1","Acvrl1","Acvr1b","Tgfbr1","Acvr1c","Acvr2a","Acvr2b","Bmpr2","Tgfbr2")
  bmp.ligands<-c("Bmp2","Bmp3","Bmp4","Bmp5","Bmp6","Bmp7",
              "Bmp8a","Gdf9","Gdf10","Gdf15")
  bmp.smads<-c("Smad1" ,"Smad2" ,"Smad3", "Smad4", "Smad5", "Smad6", "Smad7", "Smad9")

  notch.all<-c(
    "Dll1","Dll3","Dll4","Dtx1","Jag1","Jag2","Adam10","Psen1","Psen2","Psenen","Notch1",
    "Notch2","Notch3","Notch4","Mfng","Rfng","Lfng")


    if(pathway =="bmp"){
      genes.plot = c(bmp.receptors,bmp.ligands,bmp.smads)
    }else if(length(grep("rand",pathway))>0){
          a = fread( paste( "../pathways/random/",pathway,sep=""))
          genes.plot = a$V2[-1]
    }else{ #real pathway

      a = fread( paste("../pathways/", pathway, ".csv",sep=""))
      genes.plot = a$To
    }

    # }else if(pathway=="notch"){
    #   genes.plot = notch.all
    # }else if(pathway =="")

    return (genes.plot)
}
#optional:

get.pathway.expression<-function( pathway.genes_,  min.frac.genes.expressed = 0.1, fold.nreads = 2){
  # fold.nreads :   min(nreads) > fold.nreads * length(pathway)
  #let's take the raw data to manual normalization then SC3
  #manual.norm<-tiss@raw.data
  manual.norm<-tiss[["RNA"]]@counts
  #take pathway genes to filter cells:
  manual.norm.pathway<-manual.norm[pathway.genes_,]
  #at least twice the number of genes for nReads
  count.cells = Matrix::colSums(manual.norm.pathway) #genes per cell
  #min nunmber genes per cell:
  manual.norm.pathway = manual.norm.pathway[, count.cells>length(pathway.genes_)* fold.nreads]
  #at least 3 genes have >0 expression values
  genes.at.least.one = Matrix::colSums(manual.norm.pathway>0)

  manual.norm.pathway = manual.norm.pathway[, genes.at.least.one>= min.genes.per.cell]
  #NOW: let's take the cells in the new filtered matrix from the normalized data in seurat.pathway@data
  #

  min.genes.per.cell = round(length(pathway.genes_) * min.frac.genes.expressed)
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


#all.pathways = c("glucose","notch","shh","inositol","mtorc")
# LOAD all pathway names
all.pathways = list.files("../pathways/")
all.pathways = all.pathways[grep(".csv",all.pathways)]
all.pathways = gsub(pattern = "\\.csv$","",all.pathways)

all.random = list.files("../pathways/random")
all.random = all.random[grep(".csv",all.random)]
all.random = gsub(pattern = "\\.csv$","",all.random)


#main clustering, plotting , pca ALL pipeline
#pca.manual = c(9,11,14,13,8)

#for( p in 1:length(all.pathways)){
full.pipeline<-function(which.pathway,plot.heatmaps = F,plot.tsne = F,  plot.elbow = T){
    #which.pathway = all.pathways[p]

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
     pcs.compute = 14
     pcs.compute = round(length(pathway.genes_) * 3/4)
     seurat.pathway <- RunPCA(object = seurat.pathway, features =pathway.genes_, do.print = FALSE, npcs = pcs.compute,maxit =10000) #no print

    # seurat.pathway <- ProjectPCA(object = seurat.pathway, do.print = FALSE) #no print



     #pca.used = pca.manual[p]

     # seurat.pathway = JackStraw(seurat.pathway, num.replicate = 200, prop.freq = .1)
     # seurat.pathway = ScoreJackStraw(seurat.pathway, dims = 1:pcs.compute)
     #
     #
     # filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_JackStraw_", toString(pcs.compute), "_.pdf",sep="")
     # pdf(filename)
     #  a = JackStrawPlot(seurat.pathway, dims = 1:pcs.compute)
     #  plot(a)
     # dev.off()


     #plot call
     if(plot.elbow == T){
         filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_elbow_", toString(pcs.compute), "_.pdf",sep="")
         a= ElbowPlot(object = seurat.pathway, ndims = pcs.compute)
         pdf(filename)
           plot(a)
         dev.off()
     }

    # min.pcn = 4
    # diff(a$data$data.use[-c(1:min.pcn)])


    res.used = 0.9
    seurat.pathway = FindNeighbors(seurat.pathway, dims = 1:pcs.compute) #for Seurat3

  #  seurat.pathway <- FindClusters(object = seurat.pathway, reduction.type = "pca", dims.use = 1:pca.used,
  #                        resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc=T ) #, plot.SNN =T) #DONE

    seurat.pathway = FindClusters(seurat.pathway, resolution = 0.5)
    #FindClusters now automatically sets the ident (and erases other idents)

  #  seurat.pathway<-SetIdent(seurat.pathway,ident=seurat.pathway@meta.data$res.0.9) #just in case
    #V2
    #seurat.pathway <- RunTSNE(object = seurat.pathway, dims.use = 1:pca.used, seed.use = 10, perplexity=30,
    #                                       check_duplicates = F)
    #v3:
    seurat.pathway<-RunTSNE(seurat.pathway, dims = 1:pcs.compute, perplexity = 30)
    seurat.pathway = RunUMAP(seurat.pathway, dims = 1:pcs.compute)

    if( plot.tsne ==T){
        #plot call :
        filename = paste("../results/tabulaMurisSeurat/", which.pathway, "_TSNE_", toString(pcs.compute), "_.pdf",sep="")
        pdf(filename)
        # TSNEPlot(seurat.pathway)
        a = DimPlot(seurat.pathway,reduction="tsne")
        plot(a)
        a = DimPlot(seurat.pathway,reduction="umap")
        plot(a)
        dev.off()
    }

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


    if(  length(grep(".*rand.*",which.pathway))>0 ){
    #Save random pathways in a different directory
      save(seurat.pathway,file = paste("../datasets/TabulaMuris_bmp/random/", which.pathway, "_clusteredOK_NoVarGenes_04082019.rda",sep=""))
    }else{
      save(seurat.pathway,file = paste("../datasets/TabulaMuris_bmp/", which.pathway, "_clusteredOK_NoVarGenes_04082019.rda",sep=""))
    }

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
do.pca.from.list<-function(which.pathway = "",maxit =10000){
  pathway.genes_ = pathway.genes(pathway = which.pathway)
  pathway.genes_ = pathway.genes_[pathway.genes_ %in% rownames(tiss[["RNA"]]@counts)]

  res.list = get.pathway.expression( pathway.genes_)
  counts.pathway = res.list[[1]]
  norm.pathway = res.list[[2]]
  meta.data = res.list[[3]]

  #function call
   seurat.pathway = setup.seurat( counts.pathway, norm.pathway, meta.data) #works v3


   seurat.pathway <- ScaleData(object = seurat.pathway, features = rownames(seurat.pathway))
   #pcs.compute = length(pathway.genes_) -1
   #seurat.pathway <- RunPCA(object = seurat.pathway, features =pathway.genes_, do.print = FALSE, npcs = pcs.compute,maxit =maxit) #no print
   pca.res = princomp(t(as.matrix(seurat.pathway[["RNA"]]@data)))
   cor.matrix = cor(t(as.matrix(seurat.pathway[["RNA"]]@data)))

   if(length(grep("rand",which.pathway))>0){ #save all random pathways in a subfolder
     save(seurat.pathway, pca.res, cor.matrix,file = paste("../datasets/TabulaMuris_bmp/random/", which.pathway, "_clusteredOK_NoVarGenes_04082019.rda",sep=""))

   }else{
     save(seurat.pathway, pca.res,cor.matrix, file = paste("../datasets/TabulaMuris_bmp/", which.pathway, "_clusteredOK_NoVarGenes_04082019.rda",sep=""))
   }
   print( paste( "DONE ", which.pathway, "...\n"))

   return(seurat.pathway)
}
#############################################
#  n.cores  = 10
#  cl <- parallel::makeForkCluster(n.cores)
#  doParallel::registerDoParallel(cl)

do.pca.all = function(pathway.list){

    results  = foreach(p = pathway.list) %dopar% do.pca(which.pathway =p)

}


do.pca.all.list = function(pathway.list){

    results  = foreach(p = pathway.list) %dopar% do.pca.from.list(which.pathway =p)

}


#    parallel::stopCluster(cl)

#######
#######
# DATA analysis
extractSeuratList<-function(all.random = "", pc.cutoff = 0.5){
    npc = array(0, length(all.random))
    n.clusters = array(0, length(all.random))
    tot.var = array(0, length(all.random))
    npc.list.nat = list()
    for(i in 1:length(all.random)){
        load(all.random[i])
        scaled = seurat.pathway[['RNA']]@scale.data
        pca.scaled = princomp(t(scaled))
        pct.var =pca.scaled$sdev^2/sum(pca.scaled$sdev^2)
        tot.var[i] =sum(pca.scaled$sdev^2)
        npc[i] =  length(which(cumsum(pct.var)< pc.cutoff))
        n.clusters[i]  =length(unique(seurat.pathway$seurat_clusters))
        npc.list.nat[[i]] = pct.var
        rm(seurat.pathway)
    }

    ngenes =do.call("rbind",lapply(npc.list.nat, length))
    cluster.data = data.frame(clusters = n.clusters, ngenes = ngenes, npc = npc, tot.var = tot.var)
    cluster.data =as_tibble(cluster.data)

    return(list(cluster.data, npc.list.nat))
}

## use as:   > p = plot.nat.vs.rand(cluster.stats,cluster.data.nat)
plot.nat.vs.rand <-function(cluster.stats,cluster.data.nat,which.var = "clusters"){

    if(which.var =="clusters"){
       p2 = ggplot(data = cluster.data.nat, aes(x = ngenes,y=clusters,label = name)) +
            geom_point() + scale_color_manual(values=c('#999999'))+ geom_text_repel()+ theme_classic(base_size =18)

       p2 = p2 + geom_point(data=cluster.stats,aes(ngenes,clusters,colour ="red")) +
          geom_line(data=cluster.stats,aes(ngenes,clusters,colour ="red")) +
              geom_errorbar(data=cluster.stats,aes(ymin=clusters-sd.clusters,ymax=clusters+sd.clusters),width=0.2,position=position_dodge(0.05))

    }else if(which.var =="tot.var"){
        p2 = ggplot(data = cluster.data.nat, aes(x = ngenes,y=tot.var,label = name)) +
             geom_point() + scale_color_manual(values=c('#999999'))+ geom_text_repel()+ theme_classic(base_size =18)

        p2 =  p2 + geom_point(data=cluster.stats,aes(ngenes,tot.var,colour ="red")) +
           geom_line(data=cluster.stats,aes(ngenes,tot.var,colour ="red")) +
               geom_errorbar(data=cluster.stats,aes(ymin=tot.var-sd.tot.var,ymax=tot.var+sd.tot.var),width=0.2,position=position_dodge(0.05))


    }else if(which.var=="npc"){
        p2 = ggplot(data = cluster.data.nat, aes(x = ngenes,y=npc,label = name)) +
             geom_point() + scale_color_manual(values=c('#999999'))+ geom_text_repel()+ theme_classic(base_size =18)

        p2 = p2 + geom_point(data=cluster.stats,aes(ngenes,npc,colour ="red")) +
           geom_line(data=cluster.stats,aes(ngenes,npc,colour ="red")) +
               geom_errorbar(data=cluster.stats,aes(ymin=npc-sd.npc,ymax=npc+sd.npc),width=0.2,position=position_dodge(0.05))


    }

 return(p2)
}
