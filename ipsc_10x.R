


setwd( "/home/agranado/MEGA/Caltech/rnaseq/10x_iPSC/countMat_cellranger")


filterRawCounts<-function(rawdata, min.reads = 2000, min.cells.frac = 0.005, quantile.mito = 0.85, remove.mt = T){

    rawdata = rawdata[,Matrix::colSums(rawdata)>min.reads]
    rawdata = rawdata[Matrix::rowSums(rawdata>0  )> dim(rawdata)[2] *min.cells.frac,]

    mito.genes <- grep(pattern = "^MT-", x = rownames(x = rawdata), value = TRUE)


    percent.mito <- Matrix::colSums(rawdata[mito.genes, ]) / Matrix::colSums(rawdata)
    rawdata  = rawdata[,percent.mito<quantile(percent.mito,probs = quantile.mito)]
    if(remove.mt) rawdata = rawdata[!rownames(rawdata) %in% mito.genes,]

    hist(percent.mito)

    return(rawdata)
}

### READ AND LOAD SEURAT ############################

rawdata.diff <- Read10X(data.dir="Diff/")
rawdata.diff<-filterRawCounts(rawdata.diff, remove.mt = T)


data.diff <- CreateSeuratObject(counts = rawdata.diff,
                           project = "iPSC_diff",min.features = 200,min.cells = 20)



rawdata.plur <- Read10X(data.dir="Plur/")
rawdata.plur = filterRawCounts(rawdata.plur)
data.plur <- CreateSeuratObject(counts = rawdata.plur,
                           project = "iPSC_plur",min.features = 200,min.cells = 20)

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

ipsc.combined <- ScaleData(object = ipsc.combined, verbose = FALSE, vars.to.regress =c("nCount_RNA","percent.mito"))
#no regression :
#ipsc.combined <- ScaleData(object = ipsc.combined, verbose = FALSE ) #, vars.to.regress =c("nCount_RNA","percent.mito"))
ipsc.combined <- RunPCA(object = ipsc.combined, npcs = 50, verbose = FALSE)

#
VizDimLoadings(object = ipsc.combined, dims = 1:2, reduction = "pca")
# t-SNE and Clustering
ipsc.combined <- RunUMAP(object = ipsc.combined, reduction = "pca", dims = 1:30)
ipsc.combined <- FindNeighbors(object = ipsc.combined, reduction = "pca", dims = 1:30)
ipsc.combined <- FindClusters(ipsc.combined, resolution = 0.8)


p1 <- DimPlot(object = ipsc.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = ipsc.combined, reduction = "umap", label = TRUE)
x11()
#overlap samples and clusters :
plot_grid(p1, p2)
x11()
#split by sample
DimPlot(object = ipsc.combined , reduction = "umap", split.by = "orig.ident")

################


####### METHOD 2

################


################ FROM SCRATCH
n.min.reads =3000
min.genes.per.cell =200
rawdata.diff <- Read10X(data.dir="Diff/")
rawdata.diff<-filterRawCounts(rawdata.diff, min.reads = n.min.reads, remove.mt = F)


data.diff <- CreateSeuratObject(counts = rawdata.diff,
                           project = "iPSC_diff",min.features = min.genes.per.cell,min.cells = 20)



rawdata.plur <- Read10X(data.dir="Plur/")
rawdata.plur = filterRawCounts(rawdata.plur,min.reads = n.min.reads, remove.mt = F)
data.plur <- CreateSeuratObject(counts = rawdata.plur,
                           project = "iPSC_plur",min.features = min.genes.per.cell,min.cells = 20)




#####
ipsc.combined <- merge(x = data.plur, y = data.diff, add.cell.ids = c("Plur", "Diff"), project = "iPSC")

ipsc.combined <- PercentageFeatureSet(object = ipsc.combined, pattern = "^MT-", col.name = "percent.mt")
#new normalization in seurat 3
ipsc.combined <- SCTransform(object = ipsc.combined, vars.to.regress = "percent.mt", verbose = FALSE)

####

ipsc.combined <- RunPCA(object = ipsc.combined, npcs = 50, verbose = FALSE)

VizDimLoadings(object = ipsc.combined, dims = 1:2, reduction = "pca")
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


# read markers and make plots
markers<-read.csv("markersNico.tsv",header= F, sep="\t")
names(markers)<-c("gene","celltype")
markers$gene<-as.character(markers$gene)
celltypes = unique(markers$celltype)

for(i in 1:length(celltypes) ){
  markers %>% filter(celltype==celltypes[i]) -> aa ; aa$gene
  FeaturePlot( ipsc.combined, features = aa$gene)
  ggsave(filename=paste("plots/",celltypes[i],"_test.pdf",sep=""),width  = 12, height = 12)
}
