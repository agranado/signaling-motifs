# General script for all single cell analysis sessions

library(pheatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(stylo )
library(viridis)
library(Seurat)
library(ggplot2)
library(Spectrum)
library(ClusterR)
library(tibble)
library(cluster)


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


#a = brewer.pal(9,'Blues')
#colors =colorRampPalette(a)

blues_pal<-colorRampPalette(brewer.pal(n = 9, name = 'BuPu'))

# First step: calculate the average expression profile for each global cluster
# For a given scRNA object AND a list of genes, it retrieves the average expression profiles for each cluster
# FetchData is using SCT normalization from Seurat
# All expression values used in the paper are SCT normalized. See the Seurat SCT vignette for additional information.
avg.matrix<-function(seurat.obj,genes.list,by= "seurat_clusters",upper.case = T,
                      min.expr = 0.02, min.genes.pr = 0.2, filter_clusters = F){
  cellnames = rownames(seurat.obj@meta.data)
  genenames = rownames(seurat.obj)


  if(upper.case) genes.list = toupper(genes.list)

  genes.list = genenames[which(toupper(genenames) %in% toupper(genes.list))]
  # Active Assay is SCT (you can in principle do it using the RNA slot but that is only Log normalized, though to be fair the values are almost identical)
  data.to.plot = FetchData(seurat.obj, c(genes.list, by))
  data.to.plot$cell = rownames(data.to.plot)
  #we added 2 new fields, get the gene names by excluding them (or get the before)...
  genes.plot = colnames(data.to.plot)[1:(length(colnames(data.to.plot))-2)]

  data.to.plot %>% gather( key =genes.plot, c(genes.plot), value = expression) -> data.to.plot

  data.to.plot %>% dplyr::group_by_at(c(by, "genes.plot")) %>% dplyr::summarize(avg.exp = mean(expression)) %>% spread(genes.plot,avg.exp) -> mean.expr.matrix
  mean.mat =as.matrix( mean.expr.matrix[,-1]); rownames(mean.mat)<-unlist(mean.expr.matrix[,by])

  #data.to.plot %>% dplyr::group_by_at(c(by, "genes.plot")) %>% dplyr::summarize(avg.exp = mean(expression),sd.exp = sd(expression))

  return(mean.mat)
}

# Normalization and Scaling

# scaling function min.max per column
min.maxNorm<-function(x){
    maxs = apply(x,2,max)
    mins = apply(x,2,min)
    for(i in 1:dim(x)[2]){
        x[,i] = (x[,i] - mins[i]) / (maxs[i] - mins[i])
    }
    return(x)

}


# qualitative palette


makeQualitativePal <- function(n, rand_order = T, skip = 0, tail_colors = F){

  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) # only 70 qualitative colors 
  #pie(rep(1,n), col=sample(col_vector, n))
  if(n<=length(col_vector)){
      if(rand_order ==T){
        return(sample(col_vector, n))
      }else{
        # to add diversity we can get the last n colors of the array. Useful when plotting two pathways
        if(tail_colors){
            x_col = tail(col_vector,n)
        }else{
            x_col  =col_vector[(1+skip):(n+skip)]
        }
        return(x_col)
      }
  }else{
    return(c(col_vector,col_vector[1:(n-length(col_vector))]))
  }
}


# Color palette for Tabula Muris Adult
# This palette includes colors for both 10x and FACS
# It is based on the original tissue colors from the TM paper
# Additional colors for tissues in 10x not present in FACS were added
my_tissue_colors = list(Tissue = c(

    'Bladder'  = "#aec7e8",
    'Brain_Myeloid'=  "#ff7f0e",
    'Brain_Non-Myeloid' = "#ffbb78",
    'Diaphragm' = "#2ca02c",
    'Fat' =  "#98df8a",
    'BAT' =  "#98df8a",
    'SCAT' =  "#98df8a",
    'GAT' =  "#98df8a",
    'MAT' =  "#98df8a",
    "Heart_and_Aorta"  =            "#d62728",
    'Heart' = "#d62728",
    "Kidney" =            "#ff9896",
    "Large_Intestine"  =  "#9467bd",
    "Limb_Muscle"  =      "#c5b0d5",
    "Liver"        =      "#8c564b",
    "Lung"         =       "#c49c94",
    "Mammary_Gland"     = "#e377c2",
    "Marrow"            = "#f7b6d2",
    "Pancreas"          = "#7f7f7f",
    "Skin"              = "#c7c7c7",
    "Spleen"            = "#bcbd22",
    "Thymus"            = "#dbdb8d",
    "Tongue"            = "#17becf",
    "Trachea"           = "#9edae5",
    'Skeletal' ="#D9D9D9" ,
    'Muscle' = "#c5b0d5"
        ))


# Load the list of pathway genes
#
# pathways = read.csv('/home/agranado/MEGA/Caltech/rnaseq/groupMeeting2020/knn_oldmouse/pathway_list2_oct2020.csv', header = T, colClasses = 'character')
#
# i = 4
# all_pathways = pathways$pathway %>% unique()
# pathways %>% filter(pathway == all_pathways[i]) %>% select(gene) -> this_pathway
#
# #

#pathways_file_short = '/home/agranado/MEGA/Caltech/rnaseq/groupMeeting2020/knn_oldmouse/pathway_list2_oct2020_shortname.csv'

#all_pathways = read.csv(pathways_file_short, header = T, colClasses = 'character')
#pathway_genes = all_pathways$gene
# # Fix the order for ploting Fig 1
# bmp_mat = read.csv('Bmp_adult_raw_SCRAN.csv',header =T )
# bmp_mat %>% select(all_of(bmp.order), tissue, cell_type, dataset) %>% write.csv('Bmp_adult_raw_SCRAN.csv', quote = F, row.names = F)
