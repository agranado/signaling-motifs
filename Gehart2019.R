
#library(Matrix)
library(Seurat)
library(data.table)
library(stringr)
rm(list=ls())
#data from intestine single cell rna seq
data.path = "/home/agranado/MEGA/Caltech/rnaseq/datasets/Gehart2019/GSE113561_RAW/rawcounts/"
files = paste(data.path,list.files(data.path),sep="")

f1<-fread(files[1])
  gene.names.raw =  f1[,1] #names in first colum #same order for all csv, so we can take it as then assign it to the big matrix
f1<-f1[,-1]

raw.data.list = list(Matrix(as.matrix(f1),sparse=T))
for (i in 2:length(files)){
  f2<-fread(files[i])
  f2<-f2[,-1]
  raw.data<-Matrix(as.matrix(f2),sparse=T)
  raw.data.list<-append(raw.data.list,raw.data)
#  f1<-cbind(f1,f2)
}

raw.data <- do.call(cbind, raw.data.list)

row.names(raw.data)<-gene.names.raw$V1


#before doing anything, we need to include the meta data annotated by the authors:
meta.data<-read.csv("/home/agranado/MEGA/Caltech/rnaseq/datasets/Gehart2019/GSE113561_RAW/GSE113561_SCIntOrg.CELLIDMetatable.csv")
#there is a duplicate entry here, lets remove it :
meta.data=meta.data[!duplicated(meta.data$CELLID),]

#for some reason there are more cells in the raw counts that in the meta.data, so let's remove them:
common.cells = intersect(colnames(raw.data),meta.data$CELLID)
#remove them from the count matrix
raw.data= raw.data[,which(colnames(raw.data) %in% common.cells)]
#there is one duplicated cell in the meta.data file :
meta.data=meta.data[which(meta.data$CELLID %in% common.cells), ]

#now we can match the entries by name:
row.names(meta.data)<-meta.data$CELLID
meta.data = meta.data[colnames(raw.data),]
#they have their own "exclude column"
raw.data = raw.data[,!meta.data$exclude]
#7507 cells after filtering


#plot the histogram of reads across cells:
hist(log10(Matrix::colSums(raw.data)))
#First filter: at least 2000 UNIQUE transcripts
high.count.cells = Matrix::colSums(raw.data>0)>2000 #this threshold they use in the paper
raw.data  = raw.data[,high.count.cells]


#The data consists of intestine and organoid cells so lets create two objects:
organoid<-grep(pattern="SCOrg.*",  colnames(raw.data),value=F)
raw.organoid<-raw.data[,organoid]
raw.data<-raw.data[,-organoid]



#look for more about ERCC here http://tools.thermofisher.com/content/sfs/manuals/cms_086340.pdf
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
# there are genes (rows) that come from the ERCC control, so we can find them and calibrate the quantifitation
# there should be 92 ERCC transcripts

#percent is the ratio, for each cell, between the sum of ERCC detection divided by the total count
percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
sum.ercc<-Matrix::colSums(raw.data[erccs, ])
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)

raw.data <- raw.data[-ercc.index,] #remove the ERCC sequences


#remove mitocondrial genes
#in this dataset, cells with high mitochondrial content (>50%) have been already removed
# this matrix is basically already being filtered and QC'd

mt.index <- grep(pattern = "^mt-", x = rownames(x = raw.data), value = FALSE)
percent.met <- Matrix::colSums(raw.data[mt.index, ])/Matrix::colSums(raw.data)
raw.data<-raw.data[-mt.index,]

#remove specific genes (listed by the authors)
#These genes are known to cause problem in clustering (artifacts)
remove.genes<-c("Rn45s", "Malat1", "Kcnq1ot1", "A630089N07Rik","Gm17821")
remove.pattern<-paste(remove.genes,collapse="|")
raw.data<-raw.data[-which(row.names(raw.data) %like% remove.pattern),]
#####

#this datset has the chromose number for each each, let's remove that:
#gene.names<-str_match(row.names(raw.data),"(.*)__chr.*")[,2]
#row.names(raw.data)<-gene.names
# remove exclude now? >meta.data[colnames(raw.data),]$exclude
#TIME data only exist for 1735 cells, so that's why their "exclude" filter might be important
#extract time data for each cell
cell.times<-time.data[colnames(raw.data),]$Time

#after removing the chromosome we have duplicated genes, let's remove them:
#raw.data = raw.data[-which(duplicated(row.names(raw.data))),]
# START Seurat:

# # # # ## SEURAT OBJECT
tiss <- CreateSeuratObject(raw.data = raw.data)

tiss <- AddMetaData(object = tiss, meta.data)
tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")
tiss <- AddMetaData(object = tiss, percent.met, col.name = "percent.mt")
tiss <- AddMetaData(object = tiss, cell.times,  col.name = "cell.time" )
# Change default name for sums of counts from nUMI to nReads
colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads' # this is not UMI data so Seurat calculates only the number of reads


x11()
tiss <- NormalizeData(object = tiss, scale.factor = 1e4) #default normalization by Seurat
tiss <- ScaleData(object = tiss)
tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = 4, y.cutoff = 0.5, x.low.cutoff = 0.0125)



tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 100,genes.print = 5)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)

x11()
#```{r, echo=FALSE, fig.height=4, fig.width=8}
PCHeatmap(object = tiss, pc.use = 1:3, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, num.genes = 8)
#```
x11()
PCElbowPlot(object = tiss, num.pc = 100)
n.pcs = 30 #based on PCElbowPlot
res.used <- 0.8 # mid range value for this parameter is related to the number of clusters goes from 0.6 - 1.2

x11();VizPCA(object = tiss, pcs.use = 1:4,font.size=1)
x11();PCAPlot(object = tiss, dim.1 = 1, dim.2 = 2)

#resolution: Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
    resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc=T,plot.SNN=T) #DONE



tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs, perplexity=40,
                    check_duplicates = F)
x11()
TSNEPlot(tiss)
#

tiss@meta.data['cluster'] <- tiss@ident
tiss@meta.data['cell'] <- rownames(tiss@meta.data)

#look for markers that authors found in their analysis:
features.plot = c("Sct", "Agr2", "Spink4", "Tff3",
        "Muc2", "Lyz1", "Defa17", "Dll1", "Neurog3","Neurod1","Reg4","Chga","Tac1","Tph1")
#another set of markers they mention in the paper
features.plot <- c("Cck","Gcg","Nts","Sst","Ghrl","Tac1","Tph1","Gip","Neurog3","Neurod1","Isl1","Reg4","Chga")
features.plot =paste(features.plot,"_",sep="") #to avoid similar named genes
cellMarkers<- paste(features.plot,collapse="|")
plot.markers = row.names(raw.data)[which(row.names(raw.data) %like% cellMarkers)]
x11()
FeaturePlot(object = tiss, features.plot = plot.markers, cols.use = c("grey", "blue"),
    reduction.use = "tsne")

# find markers for all clusters:
tiss.markers <- FindAllMarkers(object = tiss, only.pos = TRUE, min.pct = 0.25,
                                   thresh.use = 0.25)
tiss.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
##
top10 <- tiss.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top10$gene<-str_match(top10$gene,"(.*)__chr.*")[,2] #remove the chr from the gene name

# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
x11();DoHeatmap(object = tiss, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)



### CALL cell types in clusters based on known markers provided by the authors :
cell.markers<-read.csv("/home/agranado/MEGA/Caltech/rnaseq/datasets/Gehart2019/GSE113561_RAW/1-s2.0-S009286741831643X-mmc3.csv")

all.clusters = unique(top10$cluster)
cell.types = colnames(cell.markers)
cell.type.call = matrix(0,length(all.clusters),length(cell.types))
for(i in all.clusters){
  i = as.numeric(i)
  for(this.type in 1:length(cell.types)){
    cell.type.call[i+1,this.type] = length(which(top10$gene[top10$cluster==i] %in% cell.markers[,cell.types[this.type]]))
  }
}
colnames(cell.type.call)<-cell.types
row.names(cell.type.call)<- as.character(all.clusters)

###
#BMP profiles for clusters:
 bmp.receptors<-c("Bmpr1a","Bmpr1b","Acvr1","Acvrl1","Acvr1b","Tgfbr1","Acvr1c","Acvr2a","Acvr2b","Bmpr2","Tgfbr2")
 bmp.ligands<-c("Bmp2","Bmp3","Bmp4","Bmp5","Bmp6","Bmp7",
                "Bmp8a","Gdf3","Gdf9","Gdf10","Gdf11","Gdf15")
 bmp.smads<-c("Smad1" ,"Smad2" ,"Smad3", "Smad4", "Smad5", "Smad6", "Smad7", "Smad9")
 features.plot = c(bmp.receptors,bmp.ligands,bmp.smads)


 features.plot =paste(features.plot,"_",sep="") #to avoid similar named genes
 cellMarkers<- paste(features.plot,collapse="|")
 plot.markers = row.names(raw.data)[which(row.names(raw.data) %like% cellMarkers)]


## NOTCH

  notch.all<-c(
  "Dll1","Dll3","Dll4","Dtx1","Jag1","Jag2","Adam10","Psen1","Psen2","Psenen","Notch1","Notch2","Notch3","Notch4","Mfng","Lfng","Rfng")
  features.plot = notch.all

  features.plot =paste(features.plot,"_",sep="") #to avoid similar named genes
  cellMarkers<- paste(features.plot,collapse="|")
  plot.markers = row.names(raw.data)[which(row.names(raw.data) %like% cellMarkers)]
  x11();FeaturePlot(tiss,plot.markers,cols.use = c("lightgrey","blue"))


###### MONOCLE
cds<-importCDS(tiss,import_all=T)

pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$cluster),
                                        c("0" = 'Late_progenitors',
                                          "1" = 'Early_progenitors',
                                        "2" = 'EC_early',
                                        "3" = 'L_I_N_cells',
                                        "4" = 'EC_late',
                                        "5" = 'K_cells',
                                        "6" = 'EC_late',
                                        "7" = 'Goblet_cells',
                                        "8" = 'NA',
                                        "9" = 'X_cells',
                                        "10" = 'Delta_cells'))

cell_type_color <- c("Late_progenitors" = "#E088B8",
                    "Early_progenitors" = "#46C7EF",
                    "EC_early" = "#EFAD1E",
                    "L_I_N_cells" = "#8CB3DF",
                    "EC_late" = "#53C0AD",
                    "K_cells" = "#4EB859",
                    "Goblet_cells" = "#D097C4",
                    "X_cells" = "#ACC436",
                    "Delta_cells" = "#F5918A",
                    'NA' = '#000080')


DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

cds <- estimateSizeFactors(cds)

cds <- preprocessCDS(cds, num_dim = 50)

cds <- reduceDimension(cds, reduction_method = 'UMAP')

cds <- partitionCells(cds)

cds <- learnGraph(cds,  RGE_method = 'SimplePPT')
x11()
plot_cell_trajectory(cds,
                     color_by = "cell_type2") +
                     scale_color_manual(values = cell_type_color)
#set the root
# a helper function to identify the root principal points:
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)

  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

MPP_node_ids = get_correct_root_state(cds,
                                      cell_phenotype =
                                        'cell_type2', "Early_progenitors")
cds <- orderCells(cds, root_pr_nodes = MPP_node_ids)
x11()
plot_cell_trajectory(cds)


## takes a bit
 pr_graph_test <- principalGraphTest(cds, k=3, cores=8)

###
#plot and svae all trajectories:
for (i in 1:length(plot.markers)){

  pdf(paste(plots.path,"notch/","trajectory_",plot.markers[i],".pdf",sep=""))
  plot_cell_trajectory(cds, markers = c(plot.markers[i]), use_color_gradient = TRUE)
  dev.off()
}

bmp.all = pathway.genes(pathway = "bmp")
bmp.indexes=match(bmp.all,rownames(raw.data))
