library(tidyverse)
#library(stringr)
#library(Seurat)
#library(here)

library(tidyr)
library(dplyr)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(grid)
library(gridExtra)

sce.analysis<-function(){
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
    ggsave(paste(plots.folder,"bmpReceptors_tSNE.pdf",sep="")  , width = 14, height = 7, units = "in")


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
}
# # # # # # # # # #
 # # # # # # # # #
# # # # # # # # # #

#clustering of bmp profiles using Seurat native functions

#load main tiss object from bmp_FACS_Notebook.Rmd
#Extracted from DotPlot function Seurat in Downloads/Seurat/plotting.R

cols.use = c("lightgrey", "blue")
  col.min = -2.5
  col.max = 2.5
  dot.min = 0
  dot.scale = 6
  scale.by = 'radius'
  scale.min = NA
  scale.max = NA
  #group.by
  plot.legend = FALSE
  do.return = FALSE
  x.lab.rot = FALSE

scale.func <- switch(EXPR = scale.by, size = scale_size,
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

#local seurat functions
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}
MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

#TUtorial for filtering and heatmap https://bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops2/McDavid/MAITAnalysis.html
#this function needs the big Seurat object
#the data.frame returned can be saved and loaded so this function wont be necessary every time
fetch.data<-function(tiss,genes.plot){
  data.to.plot <- data.frame(FetchData(object = tiss, vars.all = genes.plot))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- tiss@ident #extract tSNE cluster
  data.to.plot$ontology = tiss@meta.data$cell_ontology_class #extract ontology class
  data.to.plot$tissue = tiss@meta.data$tissue #tissue for each cell
  return(data.to.plot)
} # bug Jan 16 : until here works fine bugID namesBmp4

#this funciton will calculate the average across a particular variable. In normal DotPlot will do by cluster
retrieve.genes<-function(data.to.plot,genes.plot,which.var){
  #retrieves a list of genes from the Seurat object and calculates basic statistics
  #returns a tidy data.frame will all the information

      # # # # # # # # #
      #Lets make it tidy
      #here the genes are distributed as columns, cells as rows.
      #WE want to create a new column called "genes.plot"
      #The values of the DF are already the RNA counts so we just create the name "expression" and GATHER will take those values
      #Extra columns correspond to CELLS and therefore we leave them as they are, they need to map to the cells in the original rows
      # data.to.plot %>% gather(
      #   key = genes.plot,
      #   value = expression,
      #   -c(cell, id,ontology,tissue) #leave extra columns out
      # ) -> data.to.plot

      data.to.plot %>% gather(
        key=genes.plot, c(genes.plot) , #gather only the columns for genes
        value = expression) -> data.to.plot
      #NOTE if genes.plot has less genes than those in the data.to.plot object, the non selected
      #genes will NOT be gathered and therefore we need to remove them:


      # if there are genes that were not grouped bc the list is shorter
      # select only relevant fields
      data.to.plot %>% select(cell, id, cell_type,ontology,tissue,genes.plot,expression) -> data.to.plot

      data.to.plot %>% #group_by_at uses strings as arguments for indexing the data frame
        group_by_at(c(which.var,"genes.plot")) %>% # NOTE doing it like this groups all variables might work but then you can index easily
        dplyr::summarize(
          avg.exp = mean(expm1(x = expression)), # e^x -1
          pct.exp = PercentAbove(x = expression, threshold = 0), #any cell with expression >0
          avg.log.exp = mean(expression)
        ) -> data.to.plot

      # THIS doesnt have to change
      # Only affects ave.exp.scale
      data.to.plot %>%
        ungroup() %>%
        group_by(genes.plot) %>%
        mutate(avg.exp.scale = scale(x = avg.exp)) %>% #scale average expression (not log) #mutate adds a new variable to df
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


#takes the tidy data frame, converts the key-values to matrix and performs clustering
cluster.variable<-function(data.to.plot,variable="pct.exp"){
    names(data.to.plot)[1]<-"id"
    #some cells don't have ontology class assigned and they show NA, so we make it string
    data.to.plot$id[is.na(data.to.plot$id)] = "NA"

    all.clusters.id = unique(data.to.plot$id) # # # # NOTE HERE the index is the one you used for retrieving the data, so be careful
    dat.matrix = matrix(0,length(all.clusters.id),length( levels(data.to.plot$genes.plot)  ))
    colnames(dat.matrix)<-levels(data.to.plot$genes.plot)
    row.names(dat.matrix)=unique(data.to.plot$id)

    for (c in 1:length(all.clusters.id)){
      dat.matrix[c,]=as.numeric(spread(data.to.plot[data.to.plot$id==all.clusters.id[c],c("genes.plot",variable)],genes.plot,variable))
    }

    return(dat.matrix)
}

check.n.extract<-function(df.file){
  if(file.exists(df.file)){
    load(df.file)
  }else{
    #create the file by fetching from the Seurat object
    print("Fetching Data Frame from Seurat object .../n")
    data.to.plot=fetch.data(tiss,genes.plot)
    save(data.to.plot, file = df.file) #and save
  }
}

##############
#MAIN function:

load.data<-function(pathway="bmp",  which.var = "ontology" ,quant.var = "pct.exp"){
  #choose one variable and plot the heatmap, load previously saved data.to.plot

  if(pathway == "bmp"){
    f.index =1

  }else if(pathway=="notch"){
    f.index = 2

  }else
    stop("pathway does not exist ")



  df.files =c("data.to.plot","data.to.notch") #extracted from tiss based on a list of genes
  df.file = paste(df.files[f.index],".rdata",sep="")
  load(df.file)
  data.to.plot = eval(parse(text=df.files[f.index])) #the way I saved the variables have different names
  #let's create a new variable combining cell type and tissue (as Sarah did)
  data.to.plot %>% unite(col="cell_type",c(tissue,ontology),sep=",",remove=F) -> data.to.plot

  #manual annotation from tabula muris paper
  #percent of cells with positive expression values, counts fraction of cells > 0 counts
  genes.plot = pathway.genes(pathway)
  print("Creating tidy data frame.../n")
  data.to.plot  = retrieve.genes(data.to.plot,genes.plot,which.var) #bug ID namesBmp4: until here fine
  dat.matrix = cluster.variable(data.to.plot,quant.var)
  return(list(data.to.plot,dat.matrix))
}
#get the annotated genes (manually curated lists)
pathway.genes<-function(pathway ="bmp"){
  bmp.receptors<-c("Bmpr1a","Bmpr1b","Acvr1","Acvrl1","Acvr1b","Tgfbr1","Acvr1c","Acvr2a","Acvr2b","Bmpr2","Tgfbr2")
  bmp.ligands<-c("Bmp2","Bmp3","Bmp4","Bmp5","Bmp6","Bmp7",
              "Bmp8a","Gdf3","Gdf9","Gdf10","Gdf11","Gdf15")
  bmp.smads<-c("Smad1" ,"Smad2" ,"Smad3", "Smad4", "Smad5", "Smad6", "Smad7", "Smad9")

  notch.all<-c(
  "Dll1",
  "Dll3",
  "Dll4",
  "Dtx1",
  "Jag1",
  "Jag2",
  "Adam10",
  "Psen1",
  "Psen2",
  "Psenen",
  "Notch1",
  "Notch2",
  "Notch3",
  "Notch4")


    if(pathway =="bmp"){
      genes.plot = c(bmp.receptors,bmp.ligands,bmp.smads)
    }else if(pathway=="notch"){
      genes.plot = notch.all
    }

    return (genes.plot)
}
#optional:
#threshold is the sum across tissues for pct.exp  cells, it is a problem with bmp:
barplot.sort.filter<-function(dat.matrix,threshold = 100){
    df.matrix =as.data.frame(colSums(dat.matrix)) #total pct expressed
    colnames(df.matrix)[1]="sum.pct" #name variable

    df.matrix$gene.name = row.names(df.matrix) #create gene column
    df.matrix =df.matrix[with(df.matrix,order(sum.pct)), ] #sort based on the statistic

    #factor the column such that ggplot does not plot it in alphabetical other (default behavior)
    #the order in the factor/levels will be the order that ggplot will use in the x-axis
    df.matrix$gene.name = factor(df.matrix$gene.name,levels = df.matrix$gene.name)

    x11()
    # Basic barplot
    p<-ggplot(data=df.matrix, aes(x=gene.name, y=sum.pct)) +
      geom_bar(stat="identity")
    p
    # Horizontal bar plot
    p + coord_flip()

    #filter and returned
    pass.genes = as.vector(df.matrix[df.matrix$sum.pct>100,]$gene.name)


    new.matrix = dat.matrix[,pass.genes]
    return(new.matrix)
}

# simple bar plot
# for a list of genes and one attribute it will plot a sorted barplot
simple.barplot.sort<-function(df){
  colnames(df)<-c("x")
  df$gene.name = row.names(df)
  df=df[with(df,order(x)), ]
  df$gene.name = factor(df$gene.name, levels = df$gene.name)
  p<-ggplot(data =df, aes(x=gene.name,y = x))  +
  geom_bar(stat = "identity")
  p
  p + coord_flip()

}


##############
##############
#plot functions

plot.heatmap.2 <-function(dat.matrix){
    x11();heatmap.2(dat.matrix,trace = "none",col=brewer.pal(9,"YlGnBu"))
}

plot.dot.expression <-function(data.to.plot,genes.plot){
  names(data.to.plot)[1]<-"id"
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
  geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),text=element_text(size =20))
  x11();p
}

#plot heatmap using pheatmap
# kmean_k groups rows to make the heatmap smaller, it does plots (i think) the average levels for the group
plot.pheatmap<-function(dat.matrix,mode="single",kmeans_k=10,cluster.cols=T){
    scale.which = "none"

    p1=pheatmap(dat.matrix,
             show_rownames=T, cluster_cols=T, cluster_rows=T, scale=scale.which,
             cex=1.3, clustering_distance_rows="euclidean", cex=1,
             clustering_distance_cols="euclidean", clustering_method="complete",kmeans_k=kmeans_k,
              cluster.cols=cluster.cols)

    p2=pheatmap(dat.matrix,
            show_rownames=T, cluster_cols=T, cluster_rows=T, scale=scale.which,
            cex=1.3, clustering_distance_rows="euclidean", cex=1,
            clustering_distance_cols="euclidean", clustering_method="complete",
            cluster.cols=cluster.cols)

    x11()
    if(mode=="single"){
      p2
    }else if(mode=="kmeans"){
      p1
    }else if(mode=="both"){
      plot_list = list(p1[[4]],p2[[4]])

      g<-do.call(grid.arrange,plot_list)
    }
}
#bug ID nameBmp4 SOLVED Jan 16th
### let's filter the matrix
#


#### CLUSTER HEATMAP BY TYPE / TISSUE
heatmap.annotations<-function(dat.matrix,color.file1 ="celltype.colors.csv" ,color.file2="tissuecolor.csv"){

  class.combinations = as.data.frame(do.call(rbind,strsplit(row.names(dat.matrix),",")))
  row.names(class.combinations)<-row.names(dat.matrix)
  names(class.combinations)<-c("Tissue","Type")

  #colors
  annotation.colors = list()
  celltype.colors<-as.vector(read.csv(color.file1,header=F))
  celltype.colors = celltype.colors$V1
  names(celltype.colors)<-levels(class.combinations$Type)

  annotation.colors$Type  =celltype.colors

  tissue.colors<-as.vector(read.csv(color.file2,header=F))
  tissue.colors = tissue.colors$V1
  names(tissue.colors)<-levels(class.combinations$Tissue)

  annotation.colors$Tissue  =tissue.colors

  return(list(class.combinations,annotation.colors))
}


fancy.heatmap <-function(dat.matrix, nclusters = 15,cex.row = 2,cex.col = 3){

  #get annotations: (with default colors)
  ann.list = heatmap.annotations(dat.matrix)
  class.combinations = ann.list[[1]]
  annotation.colors = ann.list[[2]]
  #Format: cell type (tissue)
  row.names(class.combinations)<-paste(class.combinations$Type," (",class.combinations$Tissue,")",sep="")
  row.names(dat.matrix)<-paste(class.combinations$Type," (",class.combinations$Tissue,")",sep="")
  scale.which = "none"
  x11();
  p2=pheatmap(dat.matrix,
             show_rownames=T, cluster_cols=T, cluster_rows=T, scale=scale.which,
             cex=1, clustering_distance_rows="euclidean", cex=2,
             clustering_distance_cols="euclidean", clustering_method="complete",
             annotation_row = class.combinations,cutree_rows = nclusters,
             annotation_colors = annotation.colors,
             fontsize_row = cex.row, fontsize_col = cex.col,
              annotation_legend=F)
}

make.legend <-function(annotation.colors,idx=1,text.size=18){
  all.tissues = names(annotation.colors[[idx]]);
  color.table = data.frame(tissue =all.tissues,ymin = 0:(length(all.tissues)-1),ymax =0:(length(all.tissues)-1))
  color.table$ymax = color.table$ymax +1
  color.table$color = as.vector(annotation.colors[[idx]])
  color.table$xmin = rep(0,length(all.tissues))
  color.table$xmax = rep(1,length(all.tissues))
  x11();ggplot(color.table, aes (xmin = xmin,xmax =xmax, ymin = ymin, ymax = ymax))+
  geom_rect(colour="white",alpha=0.9,aes(fill=tissue)) +
  scale_fill_manual(values = color.table$color,guide=FALSE) +
  geom_text(aes(x=xmax+0.1, y=(ymin+ymax)/2, label=tissue),hjust =0, size = text.size) +
  coord_cartesian(xlim = c(0, 5), # This focuses the x-axis on the range of interest
                      clip = 'off') +   # This keeps the labels from disappearing
      theme(plot.margin = unit(c(1,3,1,1), "lines")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) + xlab("") + ylab("") +
          ggtitle("Legend")
}
# # # # #
 # # # # 


pca.bi.plot<-function(pathway = "notch", which.var = "ontology",quant.var = "pct.exp"){
  #1 load data
  test.list = load.data(pathway=pathway,which.var = which.var, quant.var = quant.var)
  data.to.plot = test.list[[1]]
  dat.matrix = test.list[[2]]
  #2 load annotations
  ann.list = heatmap.annotations(dat.matrix) #with default color file parameters for NOTCH only

  groups = FALSE

  if(groups){
      g <- ggbiplot(prcomp(dat.matrix,center=T),choices=1:2, obs.scale = 1, var.scale = 1,
              groups = df.notch$Tissue, ellipse = F,
              circle = F,varname.size =5,varname.adjust = 1.8)
              g <- g + scale_color_discrete(name = '')
              g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top');
      x11();g
            }else{
      g <- ggbiplot(prcomp(dat.matrix,center=T),choices=1:2, obs.scale = 1, var.scale = 1,
              ellipse = F,
              circle = F,varname.size =5,varname.adjust = 1.8)
              g <- g + scale_color_discrete(name = '')
              g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top');
      x11();g

            }

}




distinct.colors = c("#5cc69a",
"#c05ac5",
"#5dbb4d",
"#7362cf",
"#b5b233",
"#7183ca",
"#dd8c30",
"#48afd5",
"#d24c3d",
"#379179",
"#d34787",
"#42843d",
"#b46ea9",
"#99b25f",
"#c2646f",
"#73732b",
"#d29d61",
"#9e5f2e")




#plot heatmap using pheatmap
# kmean_k groups rows to make the heatmap smaller, it does plots (i think) the average levels for the group
# dat.matrix = new.matrix
# scale.which = "none"
# p1=pheatmap(dat.matrix,
#          show_rownames=T, cluster_cols=T, cluster_rows=T, scale=scale.which,
#          cex=1, clustering_distance_rows="euclidean", cex=1,
#          clustering_distance_cols="euclidean", clustering_method="complete",kmeans_k=14)
#
# p2=pheatmap(dat.matrix,
#         show_rownames=T, cluster_cols=T, cluster_rows=T, scale=scale.which,
#         cex=1, clustering_distance_rows="euclidean", cex=1,
#         clustering_distance_cols="euclidean", clustering_method="complete")
#
# plot_list = list(p1[[4]],p2[[4]])
# x11();
# g<-do.call(grid.arrange,plot_list)
