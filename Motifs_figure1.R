
library(dplyr)
library(tidyr)

# Figure1: +
# tiss.norma is a Seurat object with the Cell Atlas
# we extract the UMAP coordinates and gene expression values for two pathway components
# By using a divergent color scale we can plot two components at the same time with quantitative information
plot2genesUMAP<-function(gene1,gene2,alpha2=1, alpha3 = 0.5, col1 = 'blue',col2 ='red'){
    aa<-FetchData(tiss.norm, c(gene1,gene2,"UMAP_1","UMAP_2"))

    names(aa)<-c("gene1","gene2","UMAP_1","UMAP_2")

    aa %>% ggplot(aes(x =UMAP_1,y =UMAP_2)) + geom_point(size = 0.1,color='lightgray') +
        theme_minimal() + theme(text=element_text(size =22)) + theme(legend.position = "none") +
        ylab('Umap     1') + xlab('Umap 2') ->p

    size_expr = 0.1
    p = p + geom_point(data = aa %>% filter(gene1>0), aes(x = UMAP_1, y=UMAP_2,color = -gene1), size = size_expr) +
        geom_point(data = aa %>% filter(gene2>0), aes(x = UMAP_1, y=UMAP_2, color = gene2), size = size_expr, alpha = alpha2) +
        geom_point(data = aa %>% filter(gene1>0 & gene2>0),aes(x = UMAP_1,y = UMAP_2,color = -gene1),size = size_expr, alpha = alpha3) +
        scale_color_gradient2(midpoint = 0, low = col1, high=col2, mid = "lightgray")

    return(p)
}

# Make:
# x11() ; plot2genesUMAP("Acvrl1","Bmpr1a", alpha2=1, alpha3 = 0.2,col1 = "blue",col2 = "red")

my_tissue_colors = list(Tissue = c(
    'Aorta' = "#1f77b4",
    'Bladder'  = "#aec7e8",
    'Brain_Myeloid'=  "#ff7f0e",
    'Brain_Non-Myeloid' = "#ffbb78",
    'Diaphragm' = "#2ca02c",
    'Fat' =  "#98df8a",
   "Heart"  =            "#d62728",
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
   "Trachea"           = "#9edae5"))

#Generate a heatmap of the expression matrix, coloured by tissue
# filtered by min expression (rowSums)
# with cell type annotations
generateCellTypeAnnotation<-function(meta_data = data.frame() ,remove = 45, remove_repeats = T ){

  meta_data  %>% group_by(cell_ontology_class,tissue,seurat_clusters) %>% count() %>% group_by(seurat_clusters) %>%   mutate(freq = n / sum(n)) %>% arrange(desc(freq,seurat_clusters)) %>% top_n(1,freq) %>% arrange(seurat_clusters) -> cell_type_top_label

  if(remove_repeats)
    cell_type_top_label[-remove,] -> cell_type_top_label # cluster 45 show 2 top cell types, choose one of them drop the other one.


  # Create 2 complex cell type labels:
  cell_type_top_label %>% mutate(cell_type = paste(seurat_clusters," ", tissue,": ",cell_ontology_class,sep="")) -> cell_type_top_label
  cell_type_top_label %>% mutate(cell_type2 = paste(seurat_clusters," ",cell_ontology_class,sep="")) -> cell_type_top_label

  # For pheatmap to work properly:
  cell_type_ann =  cell_type_top_label %>% select(tissue, seurat_clusters, cell_type2, cell_ontology_class) %>% as.data.frame()
  row.names(cell_type_ann)<- 0:(dim(cell_type_ann)[1]-1) # This works because the clusters are arranged by seurat_clusters! See above. Otherwise, be careful.
  # Finished general part: The above code is general for all pathways

  return(cell_type_ann)
}

# Sep 29th 2020.
# More general version of the annotation function:
# For Tabula Senis, or when re-clustering data
# This function works when the meta data contains the column tissue_cluster which is the result of independently cluster each tissue and then just labeling the cells
generateCellTypeAnnotation2<-function(meta_data = data.frame() ,remove = 45, remove_repeats = T ){

  meta_data  %>% group_by(cell_ontology_class,tissue,tissue_cluster) %>%
                count() %>% group_by(tissue_cluster) %>%
                mutate(freq = n / sum(n)) %>% arrange(desc(freq,tissue_cluster)) %>%
                top_n(1,freq) %>% arrange(tissue_cluster) -> cell_type_top_label

  if(remove_repeats)
    cell_type_top_label[-remove,] -> cell_type_top_label # cluster 45 show 2 top cell types, choose one of them drop the other one.


  # Create 2 complex cell type labels:
  # tissue cluster numbers
  cell_type_top_label$tiss_cluster_n = lapply(str_split(cell_type_top_label$tissue_cluster, "_"), function(x){return(x[length(x)])})
  cell_type_top_label$tiss_cluster_n = 1:dim(cell_type_top_label)[1]

  cell_type_top_label %>% mutate(cell_type = paste(tissue_cluster,": ",cell_ontology_class,sep="")) -> cell_type_top_label
  cell_type_top_label %>% mutate(cell_type2 = paste(tiss_cluster_n," ",cell_ontology_class,sep="")) -> cell_type_top_label

  # For pheatmap to work properly:
  cell_type_ann =  cell_type_top_label %>% select(tissue, tissue_cluster, cell_type2, cell_ontology_class) %>% as.data.frame()
  row.names(cell_type_ann)<- cell_type_ann$tissue_cluster # This works because the clusters are arranged by seurat_clusters! See above. Otherwise, be careful.
  # Finished general part: The above code is general for all pathways

  return(cell_type_ann)
}




quantile_breaks <- function(xs, n = 50) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
}

# This function generates a heatmap with the name of cell types as rows and the tissue as color annotaiton
# For Tabula Muris, it uses the my_tissue_color palette
exprMatrixAsAnnotatedHeatmap<-function(expr.mat = bmp.mat, annotation_df = cell_type_top_label, min.expr = 2,
                            ann_colors = my_tissue_colors, cluster_col = T, quantile_norm = F, k = 10,
                            n_breaks = 10, cluster_row = T, return_plot = F , silent_plot =F, fontsize = 7, col_genes = T){


      # The matrix now was cell types as row.names
      # We use the simple cell type since we are adding colours by tissue
      row.names(expr.mat )<-annotation_df$cell_type2
      expr.mat.fil = expr.mat [ rowSums(expr.mat ) >=min.expr,]




      # After filtering, we loose the row.names!! Not desired, but we can fix by renaming the rows again
      cell_type_ann_fil = annotation_df %>% filter(cell_type2 %in% row.names(expr.mat.fil))
      cell_type_ann_pheat <- data.frame(Tissue = cell_type_ann_fil$tissue)
      row.names(cell_type_ann_pheat) <- cell_type_ann_fil$cell_type2


      # Ready for pheatmap
      # Here we filter the bmp.matrix for cell types with at least 2 average log counts (sum over the 7 profiles)
      #  This filter gets rid of cell types with very low (almost noise) expression that DO affect the cosine distance!

      if(length(ann_colors)==0){
        my_tissue_colors = NA
      }


      if(quantile_norm){
        blues_pal<-colorRampPalette(brewer.pal(n = 9, name = 'Blues'))

        mat_breaks <- quantile_breaks(expr.mat.fil, n = n_breaks)
        p2  = pheatmap(
            mat               = expr.mat.fil,
            color             = blues_pal(length(mat_breaks) - 1),
            cluster_cols      = cluster_col,
            breaks            = mat_breaks,
            border_color      = NA,
            clustering_distance_rows = dist.cosine(expr.mat.fil),
            show_rownames     = T,
            drop_levels       = TRUE,
            fontsize          =fontsize,
            annotation_row = cell_type_ann_pheat,
            annotation_colors = my_tissue_colors,
            cutree_rows = k
        )


      }else{
        # Make heatmaps using genes (features) as columns. Default
        if(col_genes == T){
          p2 =   pheatmap(expr.mat.fil, clustering_distance_rows = dist.cosine(expr.mat.fil),
          annotation_row = cell_type_ann_pheat, annotation_colors = my_tissue_colors, fontsize = fontsize, cluster_cols = cluster_col,
            cluster_rows = cluster_row, cutree_rows = k, silent = silent_plot)
        }else{
          # Transpose matrix
          # This is better when we have many clusters (>100)
          p2 =   pheatmap(t(expr.mat.fil), clustering_distance_cols = dist.cosine(expr.mat.fil),
          annotation_col = cell_type_ann_pheat, annotation_colors = my_tissue_colors, fontsize = fontsize, cluster_rows = cluster_col,
            cluster_cols = cluster_row, cutree_cols = k, silent = silent_plot)
        }

      }

      if(return_plot){
        return(p2)
      }


}


### Create variables:

meta_data = tiss.norm@meta.data
heatmap_ann = generateCellTypeAnnotation(meta_data = meta_data)
