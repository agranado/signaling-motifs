# first round of figures


# first round of figures



res = foreach(i =0:2, .packages = c("dplyr","mixtools","ggplot2","stringr") ) %dopar%{

    tabula %>% select(cell,tissue,cell_ontology_class,seurat_clusters) %>%
        runFullPipeline(input_matrix = notch_magic_seurat,gene.list = notch.genes,control_type=i) -> motifs_by_cluster_control1

}



## # # # #

# Pipeline 1
# Execute in parallel with foreach. Runs for 1 list of genes (pathway)
# Returns a list of dataframes grouped by seurat_clusters (cell types)
# Run the GMM and labels for a given set of genes
# Apply 2 types of controls:
#   1. Randomize the counts for each gene across cells in the matrix
#   2. Randomize the cell type (cluster ID) for all cells
# Apply k-means on the profiles such that if a cell type has two dominant motifs but they are very similar they become a single motif
# Sum over the frequencies with which each motif appears and return a data.frame
gene_list = bmp.receptors
expr_matrix = sce.seurat[['RNA']]@data[gene_list,]

# outputs 3 objects: original matrix + 2 controls 
res = foreach(i =0:2, .packages = c("dplyr","mixtools","ggplot2","stringr") ) %dopar%{

  tabula %>% select(cell,tissue,cell_ontology_class,seurat_clusters) %>%
    pipelineNov2019(gene.list = gene_list, expr_matrix = expr_matrix ,control = i) %>%
      kMeansOnMotifs(class="seurat_clusters") %>% group_by(seurat_clusters,motif_class) %>%
          summarise(freq = sum(freq)) %>% arrange(seurat_clusters, desc(freq)) -> motifs_by_cluster_kmeans


}

plotCumMotifNumbers<-function( res, pathway_name= ""){

  vals = seq(0,1,0.01)
  res_mat = do.call(rbind,lapply(res, atLeastN_motifsWithPercent,vals = vals))
  lwd = 2
  type = "l"
  plot(res_mat[1,],lwd= lwd ,type =type, col = "darkslategray4", main ="Number of clusters with at least 1 motif with % dominance", xlab="Dominance (%cells with motif)",ylab="Number of global clusters")
  lines(res_mat[2,], lwd= lwd ,type =type, col = "chocolate")
  lines(res_mat[3,], lwd= lwd ,type =type, col = "gray")
}

# should be the same as res[[1]]
plotHeatmapMotifs<-function(motifs_by_cluster_kmeans,min.pct = 0.1){

  motifs_by_cluster_kmeans %>% filter(freq>min.pct) %>% spread(key = seurat_clusters,value= freq,fill =0) -> clust_vs_motifs;
  # common for all classes
  clust_vs_motifs %>% select(-motif_class) %>% as.data.frame() -> clust_vs_motifs_mat ;

  row.names(clust_vs_motifs_mat)  = clust_vs_motifs$motif_class ;
  p= pheatmap(clust_vs_motifs_mat,fontsize = 12);
  x11(); plot(p[[4]]) ;
}

## # # # #
# Pipeline 2
tabula %>% select(cell,tissue,cell_ontology_class,seurat_clusters) %>%
  pipelineNov2019(gene.list = gene_list, expr_matrix = expr_matrix ,control = 0) %>%
    kMeansOnMotifs(class="seurat_clusters") -> motifs_by_cluster_kmeans
