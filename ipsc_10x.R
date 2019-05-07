


setwd( "/home/agranado/MEGA/Caltech/rnaseq/10x_iPSC/countMat_cellranger")


filterRawCounts<-function(rawdata, min.reads = 2000, min.cells.frac = 0.005, quantile.mito = 0.85, filter.mt.low = 0.05, remove.mt = T,projectid=""){

    rawdata = rawdata[,Matrix::colSums(rawdata)>min.reads]
    rawdata = rawdata[Matrix::rowSums(rawdata>0  )> dim(rawdata)[2] *min.cells.frac,]

    mito.genes <- grep(pattern = "^MT-", x = rownames(x = rawdata), value = TRUE)


    percent.mito <- Matrix::colSums(rawdata[mito.genes, ]) / Matrix::colSums(rawdata)
    rawdata  = rawdata[,percent.mito<quantile(percent.mito,probs = quantile.mito) & percent.mito>filter.mt.low]
    #remove the genes from the dataset?
    if(remove.mt) rawdata = rawdata[!rownames(rawdata) %in% mito.genes,]

    hist(percent.mito)


    data.diff <- CreateSeuratObject(counts = rawdata,project = projectid,min.features = 200,min.cells = 20,
                  meta.data = as.data.frame(percent.mito))


    return(data.diff)
}

### READ AND LOAD SEURAT ############################
if(mode ==1){
          rawdata.diff <- Read10X(data.dir="Diff/")
          data.diff<-filterRawCounts(rawdata.diff, remove.mt = T, min.reads = 15000,quantile.mito = 0.9, projectid = "iPSC_diff")


          rawdata.plur <- Read10X(data.dir="Plur/")
          data.plur = filterRawCounts(rawdata.plur,remove.mt = T,min.reads = 15000,quantile.mito = 0.9, projectid = "iPSC_plur")

          ### DATA PRE PROCESSING ##############################

          data.diff<-NormalizeData(object = data.diff, verbose=F)
          data.diff<-FindVariableFeatures(object=data.diff,selection.method ="vst",nfeatures =2000)

          data.plur<-NormalizeData(object = data.plur, verbose=F)
          data.plur<-FindVariableFeatures(object=data.plur,selection.method ="vst",nfeatures =2000)


          ## DATA INTEGRATION ###################################
          ipsc.anchors <- FindIntegrationAnchors(object.list = list(data.plur, data.diff), dims = 1:50)
          ipsc.combined <- IntegrateData(anchorset =ipsc.anchors, dims = 1:50)



          # Run the standard workflow for visualization and clustering
          #with regression

          #ipsc.combined <- ScaleData(object = ipsc.combined, verbose = FALSE, vars.to.regress ="percent.mito")
          #no regression :
          ipsc.combined <- ScaleData(object = ipsc.combined, verbose = FALSE ,vars.to.regress =c("nCount_RNA","percent.mito"))
          ipsc.combined <- RunPCA(object = ipsc.combined, npcs = 100, verbose = FALSE)

          #
          x11()
          VizDimLoadings(object = ipsc.combined, dims = 1:4, reduction = "pca")
          x11()
          ElbowPlot(ipsc.combined,ndims = 50)
          # t-SNE and Clustering
          ipsc.combined <- RunUMAP(object = ipsc.combined, reduction = "pca", dims = 1:40)
          ipsc.combined <- FindNeighbors(object = ipsc.combined, reduction = "pca", dims = 1:40)
          ipsc.combined <- FindClusters(ipsc.combined, resolution = 0.8)

          #switch to RNA as default assay:
          DefaultAssay(object = ipsc.combined) <- "RNA"

          #Find conserved markers:
          FindConservedMarkers(object = ipsc.combined, ident.1 = 7, grouping.var = "orig.ident",
                                verbose = FALSE)

          p1 <- DimPlot(object = ipsc.combined, reduction = "umap", group.by = "orig.ident")
          p2 <- DimPlot(object = ipsc.combined, reduction = "umap", label = TRUE)
          x11()
          #overlap samples and clusters :
          plot_grid(p1, p2)
          x11()
          #split by sample
          DimPlot(object = ipsc.combined , reduction = "umap", split.by = "orig.ident")

}





# bmp profiles :
bmp.ligands<-c("Bmp2","Bmp3","Bmp4","Bmp5","Bmp6","Bmp7",
                "Bmp8a","Gdf9","Gdf10","Gdf15")
bmp.receptors<-c("Bmpr1a","Bmpr1b","Acvr1","Acvrl1","Acvr1b","Tgfbr1","Acvr1c","Acvr2a","Acvr2b","Bmpr2","Tgfbr2","Alk1")

notch.all<-c( "Dll1",   "Dll3"   ,"Dll4"   ,"Dtx1",   "Jag1"   ,"Jag2",   "Notch1","Notch2" ,"Notch3", "Notch4" ,"Mfng",   "Rfng",   "Lfng")



x11()
mean.mat = avg.matrix(ipsc.combined,toupper(bmp.receptors))
mean.mat[mean.mat>1.5] =1.5
ann.markers = get.clus.ident(ipsc.combined,integrated = T)

pheatmap(mean.mat, annotation_row= ann.markers)

#ligands
x11()
mean.mat = avg.matrix(ipsc.combined,toupper(bmp.ligands))
mean.mat[mean.mat>1.5] =1.5
ann.markers = get.clus.ident(ipsc.combined,integrated = T)

pheatmap(mean.mat, annotation_row= ann.markers)



#NOTCH

x11()
mean.mat = avg.matrix(ipsc.combined,toupper(notch.all))

ann.markers = get.clus.ident(ipsc.combined,integrated = T)
mean.mat[mean.mat>1.5] =1.5
pheatmap(mean.mat, annotation_row= ann.markers)

################
x11()
markers<-read.csv("markersNico.tsv",header= F, sep="\t")
mean.mat = avg.matrix(ipsc.combined,toupper(markers$V1))

ann.markers = get.clus.ident(ipsc.combined,integrated = T)
mean.mat = mean.mat[,colSums(mean.mat)>1]; #filter low expressed genes
mean.mat = mean.mat[,colSums(mean.mat)<35] #filter too highly expressed genes
pheatmap(mean.mat, annotation_row= ann.markers)















####### METHOD 2

################


###############

################ FROM SCRATCH USinG concatenation and SCT normlization (gives 2 clusters)

#USAGE
#filterRawCounts<-function(rawdata, min.reads = 2000, min.cells.frac = 0.005, quantile.mito = 0.85, remove.mt = T){

rawdata.diff <- Read10X(data.dir="Diff/")
data.diff<-filterRawCounts(rawdata.diff, remove.mt = T, min.reads = 15000,quantile.mito = 0.9, projectid = "iPSC_diff")


rawdata.plur <- Read10X(data.dir="Plur/")
data.plur = filterRawCounts(rawdata.plur,remove.mt = T,min.reads = 15000,quantile.mito = 0.9, projectid = "iPSC_plur")





#####
ipsc.combined <- merge(x = data.plur, y = data.diff, add.cell.ids = c("Plur", "Diff"), project = "iPSC")

#new normalization in seurat 3
regress.mito = T
if(regress.mito){
    ipsc.combined <- SCTransform(object = ipsc.combined, vars.to.regress = "percent.mito", verbose = FALSE)
  }else{
    ipsc.combined <- SCTransform(object = ipsc.combined,verbose = FALSE)
}

####

ipsc.combined <- RunPCA(object = ipsc.combined, npcs = 100, verbose = FALSE)
x11()
VizDimLoadings(object = ipsc.combined, dims = 1:4, reduction = "pca")
x11()
ElbowPlot(ipsc.combined)
# t-SNE and Clustering
ipsc.combined <- RunUMAP(object = ipsc.combined, reduction = "pca", dims = 1:40)
ipsc.combined <- FindNeighbors(object = ipsc.combined, reduction = "pca", dims = 1:40)
ipsc.combined <- FindClusters(ipsc.combined, resolution = 0.8)


p1 <- DimPlot(object = ipsc.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = ipsc.combined, reduction = "umap", label = TRUE)
x11()
#overlap samples and clusters :
plot_grid(p1, p2)
x11()
#split by sample
DimPlot(object = ipsc.combined , reduction = "umap", split.by = "orig.ident")

# # # # # #
 # # # # # #
# read markers and make plots
plot.features.markers<-function(ipsc.combined) {
  markers<-read.csv("markersNico.tsv",header= F, sep="\t")
  names(markers)<-c("gene","celltype")
  markers$gene<-as.character(markers$gene)
  celltypes = unique(markers$celltype)
  #Plot all features
  for(i in 1:length(celltypes) ){
    markers %>% filter(celltype==celltypes[i]) -> aa ; aa$gene
    FeaturePlot( ipsc.combined, features = aa$gene)
    ggsave(filename=paste("plots/",celltypes[i],"_test.pdf",sep=""),width  = 12, height = 12)
  }

}

avg.matrix<-function(seurat.obj,genes.list,by= "ident",upper.case = T){
  cellnames = rownames(seurat.obj@meta.data)
  genenames = rownames(seurat.obj)


  if(upper.case) genes.list = toupper(genes.list)

  genes.list = genenames[which(genenames %in% genes.list)]

  data.to.plot = FetchData(seurat.obj, c(genes.list, by))
  data.to.plot$cell = rownames(data.to.plot)
  #we added 2 new fields, get the gene names by excluding them (or get the before)...
  genes.plot = colnames(data.to.plot)[1:(length(colnames(data.to.plot))-2)]

  data.to.plot %>% gather( key =genes.plot, c(genes.plot), value = expression) -> data.to.plot

   data.to.plot %>% dplyr::group_by_at(c(by, "genes.plot")) %>% dplyr::summarize(avg.exp = mean(expression)) %>% spread(genes.plot,avg.exp) -> mean.expr.matrix
  mean.mat =as.matrix( mean.expr.matrix[,-1]); rownames(mean.mat)<-unlist(mean.expr.matrix[,by])

  data.to.plot %>% dplyr::group_by_at(c(by, "genes.plot")) %>% dplyr::summarize(avg.exp = mean(expression),sd.exp = sd(expression))

  return(mean.mat)
}


get.clus.ident<-function(ipsc.combined,integrated = T){

    n.clusters = levels(ipsc.combined[[]]$seurat_clusters)
    sample.orig = array("",length(n.clusters))
    all.clust.composition = c()
    for (i in 1:length(n.clusters)){
        ipsc.combined@meta.data %>% filter(seurat_clusters == n.clusters[i]) %>% select(orig.ident) %>% table() -> clust.composition
        sample.orig[i]  = names(clust.composition)[which.max(clust.composition)]

        all.clust.composition = rbind(all.clust.composition, clust.composition)
    }
    ann.markers = data.frame(sample = sample.orig)
    rownames(ann.markers)<-levels(ipsc.combined[[]]$seurat_clusters)


    if(integrated){
        return(ann.markers)
    }else{

        cell.numbers = as.data.frame(all.clust.composition)
        cluster.type = rep("mix",dim(cell.numbers)[1])

        cluster.type[cell.numbers$iPSC_diff/cell.numbers$iPSC_plur>=2 ] = "Diff"
        cluster.type[ cell.numbers$iPSC_diff/cell.numbers$iPSC_plur<=0.5 ] = "Plur"

        ann.markers = data.frame( sample = cluster.type)
        rownames(ann.markers)<-levels(ipsc.combined[[]]$seurat_clusters)

        return(ann.markers)
    }
}
