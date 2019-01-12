library(tidyverse)
library(stringr)
library(Seurat)
library(here)

## ------------------------------------------------------------------------
load(file=here("00_data_ingest", "11_global_robj", "FACS_all.Robj"))
#object name is tiss_FACS and it is fully annotated by the tabula muris people (after running the main script)

bmps<-grep(pattern = "Bmp", x = rownames(x = tiss_FACS@data), value = TRUE)
#extract the data for all bmp genes
bmp.counts<-tiss_FACS@data[bmps,]
#type 1
bmp.receptors<-c("Bmpr1a","Bmpr1b","Acvr1","Acvrl1","Acvr1b","Tgfbr1","Acvr1c","Acvr2a","Acvr2b","Bmpr2","Tgfbr2")
bmp.ligands<-c("Bmp2","Bmp3","Bmp4","Bmp5","Bmp6","Bmp7","Bmp10","Bmp15",
            "Bmp8a","Gdf2","Gdf1","Gdf3","Gdf5","Gdf6","Gdf7","Gdf9","Gdf10","Gdf11","Gdf15")


plots.folder = "/home/agranado/MEGA/Caltech/rnaseq/bmp-profiles"

#plot the bmp receptors on top of the tSNE created for the Nature paper
#this tSNE used ALL genes and 43 principal components
x11();FeaturePlot(tiss_FACS,bmp.receptors,cols.use = c("lightgrey","blue"),nCol=4)
ggsave((plots.folder,"bmpReceptors_tSNE.pdf",sep="")  , width = 14, height = 7, units = "in")


# # # # https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html
# # # # https://satijalab.org/seurat/conversion_vignette.html
# # # # https://satijalab.org/seurat/pbmc3k_tutorial.html
# # # # https://github.com/farrellja/URD/blob/master/Analyses/QuickStart/URD-QuickStart-AxialMesoderm.md


# Alternatively, we start from scratch and cluster only using BMP
# what I should have done from the beggining..

FACS_files = list.files(here("00_data_ingest","00_facs_raw_data","FACS"), full.names = TRUE)
#kidney data, from : https://hemberg-lab.github.io/scRNA.seq.course/tabula-muris.html

dat=read.delim(FACS_files[grep("Kidney",FACS_files)],sep=",",header=T)
rownames(dat)<-dat[,1]
dat<-dat[,-1]



cellIDs <- colnames(dat)
cell_info <- strsplit(cellIDs, "\\.")
Well <- lapply(cell_info, function(x){x[1]})
Well <- unlist(Well)
Plate <- unlist(lapply(cell_info, function(x){x[2]}))
Mouse <- unlist(lapply(cell_info, function(x){x[3]}))

#read annotations
ann=read.table("./00_data_ingest/00_facs_raw_data/annotations_FACS.csv",sep=",",header=T)
ann<-ann[match(cellIDs,ann[,1]),]
celltype<-ann[,3]

#cell type composition of tissue:
summary(factor(celltype))

#quality filters
cell_anns<-data.frame(mouse = Mouse,well=Well,type=celltype)
rownames(cell_anns)<-colnames(dat)
sce<-SingleCellExperiment(assays=list(counts = as.matrix(dat)),colData = cell_anns)

isSpike(sce, "ERCC") <- grepl("ERCC-", rownames(sce))
#remove genes that are not expressed in any cell
keep_feature<-rowSums(counts(sce)>0) >0
sce<-sce[keep_feature,]

# # # # # # # # # #
 # # # # # # # # #
# # # # # # # # # #

#clustering of bmp profiles using Seurat native functions

bmp.receptors<-c("Bmpr1a","Bmpr1b","Acvr1","Acvrl1","Acvr1b","Tgfbr1","Acvr1c","Acvr2a","Acvr2b","Bmpr2","Tgfbr2")
bmp.ligands<-c("Bmp2","Bmp3","Bmp4","Bmp5","Bmp6","Bmp7","Bmp10","Bmp15",
            "Bmp8a","Gdf2","Gdf1","Gdf3","Gdf5","Gdf6","Gdf7","Gdf9","Gdf10","Gdf11","Gdf15")
bmp.smads<-c("Smad1" ,"Smad2" ,"Smad3", "Smad4", "Smad5", "Smad6", "Smad7", "Smad9")

#load main tiss object from bmp_FACS_Notebook.Rmd
#Extracted from DotPlot function Seurat in Downloads/Seurat/plotting.R
genes.plot = bmp.ligands

cols.use = c("lightgrey", "blue")
  col.min = -2.5
  col.max = 2.5
  dot.min = 0
  dot.scale = 6
  scale.by = 'radius'
  scale.min = NA
  scale.max = NA
  group.by
  plot.legend = FALSE
  do.return = FALSE
  x.lab.rot = FALSE

scale.func <- switch(EXPR = scale.by, size = scale_size,
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))


retrieve.genes<-function(tiss,genes.plot){
  #retrieves a list of genes from the Seurat object and calculates basic statistics
  #returns a tidy data.frame will all the information
      data.to.plot <- data.frame(FetchData(object = tiss, vars.all = genes.plot))
      colnames(x = data.to.plot) <- genes.plot
      data.to.plot$cell <- rownames(x = data.to.plot)
      data.to.plot$id <- object@ident #extract tSNE cluster
      data.to.plot$ontology = object@meta.data$cell_ontology_class #extract ontology class
      data.to.plot$tissue = object@meta.data$tissue #tissue for each cell
      # # # # # # # # #
      data.to.plot %>% gather(
        key = genes.plot,
        value = expression,
        -c(cell, id,ontology,tissue)
      ) -> data.to.plot

      data.to.plot %>%
        group_by(id, tissue, ontology,genes.plot) %>% # NOTE doing it like this groups all variables might work but then you can index easily
        summarize(
          avg.exp = mean(expm1(x = expression)), # e^x -1
          pct.exp = PercentAbove(x = expression, threshold = 0) #any cell with expression >0
        ) -> data.to.plot

      data.to.plot %>%
        ungroup() %>%
        group_by(genes.plot) %>%
        mutate(avg.exp.scale = scale(x = avg.exp)) %>% #scale average expression (not log)
        mutate(avg.exp.scale = MinMax(   #make all values > abs(2.5) = 2.5
          data = avg.exp.scale,
          max = col.max,
          min = col.min
        )) ->  data.to.plot
      data.to.plot$genes.plot <- factor(
        x = data.to.plot$genes.plot,
        levels = rev(x = genes.plot)
      )
    # data.to.plot$genes.plot <- factor(
    #   x = data.to.plot$genes.plot,
    #   levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
    # )
    data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
    data.to.plot$pct.exp <- data.to.plot$pct.exp * 100

    return(data.to.plot)
}

plot.dot.expression <-function(data.to.plot,genes.plot){
p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
  geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  x11();p
}

#takes the tidy data frame, converts the key-values to matrix and performs clustering
cluster.variable<-function(data.to.plot,variable="pct.exp"){
    all.clusters.id = unique(data.to.plot$id) # # # # NOTE HERE the index is the one you used for retrieving the data, so be careful
    dat.matrix = matrix(0,length(all.clusters.id),length( levels(data.to.plot$genes.plot)  ))
    colnames(dat.matrix)<-levels(data.to.plot$genes.plot)
    row.names(dat.matrix)=levels(data.to.plot$id)

    for (c in as.numeric(all.clusters.id)){
      dat.matrix[c,]=as.numeric(spread(data.to.plot[data.to.plot$id==c-1,c("genes.plot",variable)],genes.plot,variable))
    }

    return(dat.matrix)
}

#choose one variable and plot the heatmap
dat.matrix= cluster.variable(data.to.plot,"avg.exp.scale")
x11();heatmap.2(dat.matrix,trace = "none",col=brewer.pal(9,"Blues"))

#plot heatmap using pheatmap
# kmean_k groups rows to make the heatmap smaller, it does plots (i think) the average levels for the group
x11();pheatmap(dat.matrix,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
         cex=1, clustering_distance_rows="euclidean", cex=1,
         clustering_distance_cols="euclidean", clustering_method="complete",kmeans_k=14)





#
