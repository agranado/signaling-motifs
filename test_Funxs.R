compareARI<-function(p, dge_res_list = res.list_DEsingleAWS , clust_method1_ = 'complete',clust_method2_ = 'complete', all_clusters_ = c(), min_exp_ = 1.5, dist_method_ = 'fold' , dge_= 'desingle'){
    #clust_method1 for DGE
    #clust_mehtod2 for Cosine matrix
    # test different p_values for DE metric
    pv = c(10^-1,10^-3,10^-6, 10^-9)
    all_data = list()
    for(p_v in 1:4){
        mi_list = confusionMutualInformation(res.list_ = dge_res_list, which.genes = all_pathways[[1]][[p]], test.clusters =all_clusters_, p_val_threshold = pv[p_v],
                                            min_exp = min_exp_, dist_method = dist_method_, max_k = 40, clust_method1 = clust_method1_, clust_method2 = clust_method2_, dge = dge_)
        all_data[[p_v]] = mi_list
    }
    #
    # # What to plot? c(MI, MI2, ARI)
    # idx = 3
    #
    # x11(); par(mfrow = c(1,2))
    # cols = c('black', 'red', 'green','cyan')
    # plot(all_data[[1]][[idx]], type = 'o', col = cols[1], ylim=c(0,1), main = all_pathways[[2]][[p]])
    #
    # for(i in 2:4){
    #     lines(all_data[[i]][[idx]], col = cols[i])
    #
    # }
    #
    #
    # all_data_mat = c()
    # for(i in 1:length(all_data)){
    #     all_data_mat = rbind(all_data_mat, all_data[[i]][[idx]])
    # }
    # apply(all_data_mat, 2,mean) %>% plot(type ='o')


    return(all_data)
}


silhoutteDGE<-function(input_matrix,k.max, clust_method = 'kmeans',hclust_method ='complete'){

  if(clust_method =='kmeans'){
  wss <- sapply(2:k.max,
                function(k){res = kmeans(input_matrix, k, nstart=50,iter.max = 30 );ss =silhouette(res$cluster,dist(input_matrix));return(mean(ss[,3]))})
  }else if(clust_method=='pam'){
  wss <- sapply(2:k.max, function(k){res = pam(input_matrix, k );ss =silhouette(res$clustering,dist(input_matrix));return(mean(ss[,3]))})

  }else if(clust_method =="spectral"){
  wss <- sapply(2:k.max,
                  function(k){res = spectralKmeans(input_matrix, k );ss =silhouette(res$cluster,dist(input_matrix));return(mean(ss[,3]))})
              }else if(clust_method =="cosine_kmeans"){
  wss <- sapply(2:k.max,
                    function(k){res = kmeans(dist.cosine(input_matrix) %>% as.matrix(), k );ss =silhouette(res$cluster,dist(input_matrix));return(mean(ss[,3]))})

              }else if(clust_method=="cutree"){
                # This function calculates a hierarchical clustering using cosine distance
                # Then makes partitions of the resulting tree for different values of k
                # The silhouette score is calculated using Euclidean distance
                # p = pheatmap(input_matrix, clustering_distance_rows = dist.cosine(input_matrix), silent = T)
                p = pheatmap(input_matrix,clustering_method = hclust_method , ,silent = T)
                wss <- sapply(2:k.max,
                    function(k){res = treeClust(p, k );ss =silhouette(res,input_matrix %>% dist() );return(mean(ss[,3]))})

              }else if(clust_method =="cutree_cosine"){
                  p = pheatmap(input_matrix, silent = T, clustering_distance_rows = dist.cosine(input_matrix))
                  wss <- sapply(2:k.max,
                    function(k){res = treeClust(p, k );ss =silhouette(res,dist(input_matrix));return(mean(ss[,3]))})

              }else if(clust_method=="cutree_cor"){
                  p = pheatmap(input_matrix, silent = T, clustering_distance_rows = as.dist(1-cor(t(input_matrix))))
                  wss <- sapply(2:k.max,
                  function(k){res = treeClust(p, k );ss =silhouette(res,dist(input_matrix));return(mean(ss[,3]))})

              }

  return(wss)
}

treeClust<-function(p,k){
  #p2 = pheatmap(input_matrix, clustering_distance_rows = dist.cosine(input_mat), silent = T)
  aa = as.hclust(p$tree_row)
  return(cutree(aa, k))
}



# Mutual information with control for two clusterings DGE and Cosine
# cluste_method2 is for cosine hierarchical
# For parallel:
confusionMutualInformation2<-function(p =1,res.list_ = list(), test.clusters =c(),  p_val_threshold = 10^-6, min_exp = 2,   k = 8, dist_method = 'fold',
                          clust_method1 = 'complete', clust_method2 = 'complete', max_k = 20, y_min_exp = 2, dge ='desingle', n_rand =100){
  library(infotheo)

  which.genes = all_pathways[[1]][[p]]
  # We first calculate the Differntial Expression Matrix:

  # We compute both matrices
  # A differential expression matrix (i,j = N of differentially expressed genes between cell type i and j, scaled by fold change)
  # The average expression matrix


  # This matrix is in order 0:N
  bmp.mat.raw = avg.matrix(tiss.norm,which.genes, by='seurat_clusters')

  # Filter for cell types with not enough expression
  x_fil <- rowSums(bmp.mat.raw)>=min_exp
  # Filter genes that don't have enough expression (e.g. for other pathways)
  y_fil <- colSums(bmp.mat.raw) >= y_min_exp
  # Filter the main average matrix
  bmp.mat.raw.fil = bmp.mat.raw[x_fil,y_fil]
  # Re-define the gene list
  # only the genes that passed the min expresssion threshold
  which.genes = names(which(y_fil))
  # The new function for DEsingle actually computed pairwise differential expression using 0,1,2,...N order for the clusters, so we DONT use the all_clusters list
  # but rather just the numbers 0:N-1 as per Seurat convention
  # The Matrix rows will be named correctly with the cluster IDs
  if(dge =='desingle'){
    bmp.dge.mat = distMatDEsingle(res.list_,  p_val = p_val_threshold, which.genes = which.genes, cluster_list = 0:(length(test.clusters)-1), which.method = dist_method, log_transform = F)
  }else if(dge=='wilcoxon'){
    # res.list_ object must correspond to Wilcoxon output
    # all_clusters order is important here!
    # p_val default = 10^-6
    bmp.dge.mat = distMatrix(res.list_, all_clusters, which.genes = which.genes, p_val_threshold, subset_genes = T)
  }
  # Calculate the usual average matrix for all profiles

  # Clustering for both DGE and Cosine matrices
  # Filter in DGE matrix both rows and columns
  # Make pheatmap using clust_method1
  bmp.dge.mat[names(which(x_fil)),names(which(x_fil))] %>% pheatmap(border_color = NA, clustering_method =clust_method1, col = viridis(100), silent = T) -> p_dge

  # Plot average matrix in DGE order

  #bmp.mat.raw.fil = bmp.mat.raw[p_dge$tree_row$labels[p_dge$tree_row$order], ]
  bmp.mat.raw.fil %>% pheatmap(clustering_distance_rows = dist.cosine(bmp.mat.raw.fil), clustering_method = clust_method2, silent = T) -> p_cosine

  # Mutual information for different values of k
  # Using confusion matrix as an approximation of P(x,y)
  mi = rep(0, max_k)
  mi2 = rep(0, max_k)
  ari = rep(0, max_k) # Adjusted Random Index
  ami = rep(0,max_k) # Adjusted Mutual Information package::aricode
  nmi = rep(0,max_k) # package::aricode
  cs = rep(0,max_k)  # package::clusteval
  mi_rand = rep(0, max_k)

  # Matrix for multiple repeats of randomiztion
  ami_rand = matrix(0,n_rand, max_k)
  ari_rand = matrix(0,n_rand,max_k)

  for( k in 3:max_k){

      # Partition the same tree with a different value of k in each iteration
      dge_classes = p_dge$tree_row %>% cutree(k)
      cosine_classes = p_cosine$tree_row %>% cutree(k)

      # Now we have defined the partitions into k classes for each clustering


      # calculate a confusion matrix
      # How well the classes from hierarchical clustering match those predicted by DGE matrix
      # A perfect matrix would have diagonal structure meaning perfect correspondence between the two clustering methods

      conf_matrix = matrix(0,length(unique(dge_classes)),length(unique(cosine_classes)))

      for(i in 1:k){
          map_classes<-cosine_classes[names(which(dge_classes==i))]
          for(j in 1:k){
              conf_matrix[i,j ] = sum(map_classes==j)#/length(map_classes)
          }
      }
      # Mutual information between arrays of features, it does not use the confusion matrix
      mi[k] = mutinformation(cosine_classes[names(dge_classes)], dge_classes)
      # Alternative methods that directly apply Shannon's formula:
      # Normalize by the number of cell types to create a joint distribution P(x,y)
      p_x_y = conf_matrix / sum(sum(conf_matrix))
      mi2[k] = mi.plugin(p_x_y, unit = 'log2') # Calculate MI directly from confusion matrix

      ari[k] = adjustedRandIndex(cosine_classes[names(dge_classes)], dge_classes) # Adjusted Rand Index
      ami[k] = AMI(cosine_classes[names(dge_classes)], dge_classes)
      nmi[k] = NMI(cosine_classes[names(dge_classes)], dge_classes)
      cs[k] = cluster_similarity(cosine_classes[names(dge_classes)], dge_classes, 'jaccard') # package::clusteval

      # Controls with randomized clustering
      # Randomize labels but keep the sizes of the classes
      rand_clustering = sample(dge_classes)
      names(rand_clustering) = names(dge_classes)

      mi_rand[k] = mutinformation(rand_clustering, dge_classes)


      for( r in 1:n_rand){
        #Randomize each iteration
        # This function will return a matrix for the control
        rand_clustering = sample(dge_classes)
        names(rand_clustering) = names(dge_classes)

        ami_rand[r, k] = AMI(rand_clustering, dge_classes)
        ari_rand[r,k] = adjustedRandIndex(rand_clustering,dge_classes)

      }
  }

  # Normalize by max MI possible for each number of classes:
  mi = mi / log(1:max_k)
  mi2 = mi2 / log2(1:max_k)
  mi_rand = mi_rand /log(1:max_k)

  # Silhouette score calculation on distance matrix
  silh =silhoutteDGE(bmp.dge.mat[names(which(x_fil)),names(which(x_fil))], k.max =max_k, clust_method ='cutree', hclust_method = clust_method1)


  # Returns a vector of mutual information for each value of k
  # Controls start at index 8
  return(list(mi, mi2, ari,    ami, nmi,cs,     silh, mi_rand, ami_rand,   ari_rand))

}


# Sep 24th
# Goal: generate more cell type diversity from the cell atlas by clustering each tissue independently.
# We can then add a new column to the main tiss.norm object with new labels for each tissue_cluster
# Ideally we would like ~10 clusters per tissue for a total of ~180 clusters.
clusterByTissue<-function(tiss.norm  = c() ){
  tissues <- tiss.norm@meta.data$tissue %>% unique()

  clustered_data = c()
  for(t in tissues){
    tiss.norm@meta.data %>% filter(tissue ==t) ->meta.tissue
    #get all cells in this tissue
    subset_cells = meta.tissue$cell
    subset_matrix = tiss.norm[['RNA']]@counts[,subset_cells]
    # now we can create a Seurat object with the subset of cells.
    tiss_tissue= CreateSeuratObject(counts =subset_matrix, meta.data = meta.tissue )
    # Seurat pipeline standard
    tiss_tissue %>% SCTransform(vars.to.regress = 'percent.ercc') %>%
                  RunPCA(verbose = F,npcs = 50) %>% RunUMAP(dims = 1:30) %>%
                  FindNeighbors(n.start =10) %>%
                  FindClusters(resolution = 1.5) -> tiss_tissue

    tiss_tissue@meta.data %>% select(cell, seurat_clusters, tissue) %>%
      mutate(tissue_cluster = paste(tissue,'_', seurat_clusters, sep = "")) -> clusters_tissue

    clustered_data = rbind(clustered_data, clusters_tissue)

  }

  return(clustered_data)
}


# Oct 7th 2020
# Export count matrices to csv for analysis in python
# For each pathway, we need annotated count matrix along with 2 house-keeping genes
# Differential gene expression matrix
all_pathway_file_names = c('Bmp_l','Notch', 'Bmp_Tgfb','Bmp','Wnt','Eph_r','Eph_l', 'Fgfr', 'Srsf','Lpa', 'Wnt_l')
house_keeping =  c('Malat1', 'Gapdh', 'Sdha', 'Actb')


exportExpressionData<-function(i=1, pathway_list =list() , res.list_ = list(), p_val_threshold = 10^-6, min.expr =2, k_clust = 10,
                               house_keeping =  c('Malat1', 'Gapdh', 'Sdha', 'Actb') , all_clusters = c(),
                               file.path = '/home/agranado/MEGA/Caltech/rnaseq/groupMeeting2020/knn_oldmouse/expression_export/',
                               pathway_names = all_pathway_file_names){
  # k can be further optimized in python based on cosine clustering or DGE
  # Run confusionMatrix2 function to get annotated count matrix and DGE
  #
  conf_notch_ = confusionMatrix2(res.list_ =res.list_DEsingleAWS,
                               which.genes = pathway_list[[i]],
                               test.clusters = all_clusters, k = k_clust,
                               dist_method ='euclidean',
                               clust_method1 = 'ward.D2', min_exp = min.expr)

  # Average expression for house-keeping genes
  house_matrix = avg.matrix(tiss.norm, house_keeping)

  # Left join with the pathway data before exporting
  house_df<- house_matrix %>% as.data.frame() %>% mutate(seurat_clusters =row.names(house_matrix))
  export_df <- conf_notch_[[2]] %>% left_join(house_df, by ='seurat_clusters')

  # write expression matrix
  write.csv(export_df, file =paste(file.path,pathway_names[i], '_ann_expr_matrix.csv', sep =''),quote = F, row.names = F)
  # write DGE
  write.csv(conf_notch_[[3]], file = paste(file.path, pathway_names[i], '_DGE.csv',sep = ""), quote = F)
}

# To convert gene names to standard format (first letter cap, then lower case)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Oct 7th
# Export senis data into .csv for metric learning
# All pathways
# run the function for each dataset. Careful with the cell type annotation. All other parameters should be the same.
exportSenis<-function(tiss.norm, pathway_list,
                      house_keeping =  c('Malat1', 'Gapdh', 'Sdha', 'Actb'),
                      file.path = '/home/agranado/MEGA/Caltech/rnaseq/groupMeeting2020/knn_oldmouse/expression_export/Senis/',
                      pathway_names = all_pathway_file_names,
                      cell_type_ann = data.frame(), data_id = '' ){


  for(p in 1:length(pathway_list)){
    expr.mat= avg.matrix(tiss.norm , pathway_list[[p]], by ='seurat_clusters')
    df=as.data.frame(expr.mat)
    df$seurat_clusters = row.names(expr.mat)
    # left join with metadata
    # specific for each Seurat object
    export_df = left_join(df, cell_type_ann %>% select(seurat_clusters, tissue, cell_type), by ='seurat_clusters')


    # merge with housekeeping genes
    # Average expression for house-keeping genes
    house_matrix = avg.matrix(tiss.norm, house_keeping)

    # Left join with the pathway data before exporting
    house_df<- house_matrix %>% as.data.frame() %>% mutate(seurat_clusters =row.names(house_matrix))
    export_df <- export_df %>% left_join(house_df, by ='seurat_clusters')

    write.csv(export_df, file = paste(file.path, pathway_names[p], '_ann_expr_matrix_Senis_',data_id ,'.csv', sep=''), quote = F, row.names = F)
  }

}

# Oct 8th 2020
# A more general function to export all pathways as .csv starting from a Seurat object
# Pathway genes are read from a .txt file
# We are using the general format of developmental datasets (same as Niv)
pathways_file = '/home/agranado/MEGA/Caltech/rnaseq/groupMeeting2020/knn_oldmouse/pathway_list2_oct2020.csv'
exportSeuratPathways<-function(tiss.norm = c() , pathway_file = pathways_file, all_pathway_file_names = c() , cell_type_ann = data.frame()  , upper.case = F){

    # Read the list of genes for the pathways of interest + splicing family(negative control)
    all_pathways = read.table(pathway_file, sep=',', header = T, colClasses = 'character')


    expr.mat = avg.matrix(tiss.norm , all_pathways$gene, by='seurat_clusters', upper.case = upper.case)
    if(upper.case)
      colnames(expr.mat) <- firstup(tolower(colnames(expr.mat)))


    #check for non present genes and fill with 0 matrix
    not_present = all_pathways$gene[ !(all_pathways$gene %in% firstup(tolower(colnames(expr.mat)))) ]
    not_present_zero = matrix(0, dim(expr.mat)[1], length(not_present))

    colnames(not_present_zero) <-not_present
    row.names(not_present_zero) <- row.names(expr.mat)

    #bind both matrices
    expr.mat = cbind(expr.mat , not_present_zero)

    df=as.data.frame(expr.mat)
    df$seurat_clusters = row.names(expr.mat)
    # left join with metadata
    # specific for each Seurat object
    export_df = left_join(df, cell_type_ann %>% select(seurat_clusters, tissue, cell_type), by ='seurat_clusters')


    # merge with housekeeping genes
    # Average expression for house-keeping genes
    house_matrix = avg.matrix(tiss.norm, house_keeping, upper.case = T)
    if(upper.case)
      colnames(house_matrix) <- firstup(tolower(colnames(house_matrix)))

    # Left join with the pathway data before exporting
    house_df<- house_matrix %>% as.data.frame() %>% mutate(seurat_clusters =row.names(house_matrix))
    export_df <- export_df %>% left_join(house_df, by ='seurat_clusters')


    return(export_df)

}
