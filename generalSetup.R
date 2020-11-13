# General script for all single cell analysis sessions

library(pheatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(stylo )
library(viridis)
library(Seurat)


a = brewer.pal(9,'Blues')
colors =colorRampPalette(a)


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


# Load the list of pathway genes

pathways = read.csv('/home/agranado/MEGA/Caltech/rnaseq/groupMeeting2020/knn_oldmouse/pathway_list2_oct2020.csv', header = T, colClasses = 'character')

i = 4
all_pathways = pathways$pathway %>% unique() 
pathways %>% filter(pathway == all_pathways[i]) %>% select(gene) -> this_pathway
