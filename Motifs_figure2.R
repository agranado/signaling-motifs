library(viridis)
# This functions were optimized in AWS instance
# Small modifications from original versions include:
#   distMatrix: Same function
#   nGenesDifferent: Distance metric as N of differentially expressed genes is now scaled by the log-fold change, which produces better results
#   confusionMatrix: new function that calculates a confusion matrix for a pathway, using two clustering methods
#   plotDGE: small function that plots DGE matrix as heatmap

distMatrix<-function(res.list,all_clusters,which.genes,p_val,subset_genes = T){
    dist_mat = matrix(0,length(all_clusters),length(all_clusters));
    for(i in 1:(length(all_clusters)-1)){
        for(j in (i+1):length(all_clusters)){
            dist_mat[i,j] = nGenesDifferent(res.list,i,j,p_val,which.genes = which.genes, subset_genes = subset_genes)
        }
    }
    x = dist_mat + t(dist_mat)
    colnames(x)<-all_clusters
    rownames(x)<-all_clusters
    return(x)
}


nGenesDifferent<-function(res.list = list(),i = 1,j =2, p_val_low=10^-6,which.genes = bmp.receptors, subset_genes = T){
    df = res.list[[i]][[j]];
    # If is not empty
    if(is.data.frame(df)){
        # Filter only genes of interest
        if(subset_genes ==T){
            df %>% mutate(gene = row.names(df)) %>% filter(gene %in% which.genes ) ->df
        }

        # return N DGE scaled by log fold change
        return(df %>% mutate(gene = row.names(df)) %>% filter(p_val_adj <= p_val_low) %>% summarize(d = sum(abs(avg_logFC)))  %>% pull(d) )

    }else{return(0)}

}




# Works with Wilcoxon (Seurat) Differential expression
confusionMatrix<-function(res.list_ = list(), which.genes  = seven.receptors, clust_method1 = 'complete',  min_exp = 2,   k = 8){
    # We compute both matrices
    # A differential expression matrix (i,j = N of differentially expressed genes between cell type i and j, scaled by fold change)
    # The average expression matrix
    bmp.dge.mat = distMatrix(res.list_, all_clusters, which.genes = which.genes, 10^-6, subset_genes = T)
    bmp.mat.raw = avg.matrix(tiss.norm,which.genes, by='seurat_clusters')

    # Filter for cell types with enough expression
    x_fil <- rowSums(bmp.mat.raw)>=min_exp
    # Filter in DGE matrix
    bmp.dge.mat[names(which(x_fil)),names(which(x_fil))] %>% pheatmap(cutree_rows  = k, border_color = NA, clustering_method = clust_method1, col = viridis(100)) -> p_dge
    dge_classes = p_dge$tree_row %>% cutree(k)

    # Plot average matrix in DGE order
    bmp.mat.raw.fil = bmp.mat.raw[p_dge$tree_row$labels[p_dge$tree_row$order], ]
    bmp.mat.raw.fil %>% pheatmap(clustering_distance_rows = dist.cosine(bmp.mat.raw.fil),cutree_rows = k) -> p_cosine
    cosine_classes = p_cosine$tree_row %>% cutree(k)

    # calculate a confusion matrix
    # How well the classes from hierarchical clustering match those predicted by DGE matrix
    # A perfect matrix would have diagonal structure meaning perfect correspondence between the two clustering methods

    conf_matrix = matrix(0,length(unique(dge_classes)),length(unique(cosine_classes)))

    for(i in 1:k){
        map_classes<-cosine_classes[names(which(dge_classes==i))]
        for(j in 1:k){
            conf_matrix[i,j ] = sum(map_classes==j)/length(map_classes)
        }
    }
    pheatmap(conf_matrix, clustering_method = 'ward.D2', col = inferno(100))

}


# Sep 9th 2020
# Works with DEsingle output
# New version that incorporates DEsingle output as distance matrix
# This should be the preferred method from now on.
confusionMatrix2<-function(res.list_ = list(), which.genes  = seven.receptors, test.clusters =c(),  p_val_threshold = 10^-6,
                          min_exp = 2,   k = 8, dist_method = 'fold', clust_method1 = 'complete', clust_method2 = 'complete'){
    # We compute both matrices
    # A differential expression matrix (i,j = N of differentially expressed genes between cell type i and j, scaled by fold change)
    # The average expression matrix

    # The new function for DEsingle actually computed pairwise differential expression using 0,1,2,...N order for the clusters, so we DONT use the all_clusters list
    # but rather just the numbers 0:N-1 as per Seurat convention
    # The Matrix rows will be named correctly with the cluster IDs
    bmp.dge.mat = distMatDEsingle(res.list_,  p_val = p_val_threshold, which.genes = which.genes, cluster_list = 0:(length(test.clusters)-1), which.method = dist_method, log_transform = F)
    # Calculate the usual average matrix for all profiles
    # This matrix is in order 0:N
    bmp.mat.raw = avg.matrix(tiss.norm,which.genes, by='seurat_clusters')

    # Filter for cell types with enough expression
    x_fil <- rowSums(bmp.mat.raw)>=min_exp
    # Filter in DGE matrix
    bmp.dge.mat[names(which(x_fil)),names(which(x_fil))] %>% pheatmap(cutree_rows  = k, border_color = NA, clustering_method = clust_method1, col = viridis(100)) -> p_dge
    dge_classes = p_dge$tree_row %>% cutree(k)

    # Plot average matrix in DGE order
    bmp.mat.raw.fil = bmp.mat.raw[p_dge$tree_row$labels[p_dge$tree_row$order], ]
    bmp.mat.raw.fil %>% pheatmap(clustering_distance_rows = dist.cosine(bmp.mat.raw.fil),cutree_rows = k, clustering_method = clust_method2) -> p_cosine
    cosine_classes = p_cosine$tree_row %>% cutree(k)

    # calculate a confusion matrix
    # How well the classes from hierarchical clustering match those predicted by DGE matrix
    # A perfect matrix would have diagonal structure meaning perfect correspondence between the two clustering methods

    conf_matrix = matrix(0,length(unique(dge_classes)),length(unique(cosine_classes)))

    for(i in 1:k){
        map_classes<-cosine_classes[names(which(dge_classes==i))]
        for(j in 1:k){
            conf_matrix[i,j ] = sum(map_classes==j)/length(map_classes)
        }
    }


    row.names(conf_matrix) <-1:k
    colnames(conf_matrix) <- 1:k

    pheatmap(conf_matrix, clustering_method = 'ward.D2', col = inferno(100))

    export_data = as.data.frame(bmp.mat.raw[x_fil,])
    export_data$seurat_clusters = row.names(bmp.mat.raw[x_fil,])
    export_data$dge = dge_classes
    export_data$cosine = cosine_classes


    return(list(conf_matrix, export_data,bmp.dge.mat[names(which(x_fil)),names(which(x_fil))]))
}

# Calculate the confusion matrix for different values of k
# We can then compute the Mutual information between the two clustering methods: DGE matrix vs Hierarchical clustering on cosine distance
# clust_method1 is for the DGE
# cluste_method2 is for cosine hierarchical
confusionMutualInformation<-function(res.list_ = list(), which.genes  = seven.receptors, test.clusters =c(),  p_val_threshold = 10^-6, min_exp = 2,   k = 8, dist_method = 'fold',
                          clust_method1 = 'complete', clust_method2 = 'complete', max_k = 20, y_min_exp = 2, dge ='desingle'){
  library(infotheo)

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
  ami_rand = rep(0,max_k)
  ari_rand = rep(0,max_k)

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


      ami_rand[k] = AMI(rand_clustering, dge_classes)
      ari_rand[k] = adjustedRandIndex(rand_clustering,dge_classes)
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


#
# Based on the optimal clustering (number of classes)
# This is for reference Tabula Muris!
# Create a data frame with the average expression matrix AND the label from DGE and cosine distance as two independent clsutering
# This data frame can be then used to train a KNN (or any supervised clustering method) to classify new cell types e.g. from old mouse
labelCelltypes<-function(res.list_ = list(), pathway = seven.receptors, min_exp= 2, k = 8,
                          cell_type_ann = cell_type_ann_tm){

  # 0.1 A differential expression matrix (i,j = N of differentially expressed genes between cell type i and j, scaled by fold change)
  # The average expression matrix
  bmp.dge.mat = distMatrix(res.list_, all_clusters, which.genes = pathway, 10^-6, subset_genes = T)
  bmp.mat.raw = avg.matrix(tiss.norm,pathway, by='seurat_clusters')

  # 0.2 Filter for cell types with enough expression
  x_fil <- rowSums(bmp.mat.raw)>=min_exp

  # 1. compute labels for cosine distance
  bmp.mat.raw.fil = bmp.mat.raw[x_fil,]
  bmp.mat.raw.fil %>% pheatmap(clustering_distance_rows = dist.cosine(bmp.mat.raw.fil),cutree_rows = k, silent = T) -> p_cosine
  cosine_classes = p_cosine$tree_row %>% cutree(k)

  # 2. compute labels for DGE
  # Subset by x_fil
  bmp.dge.mat[names(which(x_fil)),names(which(x_fil))] %>% pheatmap(cutree_rows  = k, border_color = NA, clustering_method = 'complete', col = viridis(100),silent = T) -> p_dge
  dge_classes = p_dge$tree_row %>% cutree(k)

  # 3. Map to original matrix and return data frame

  ann_data = cbind(bmp.mat.raw.fil %>% as.data.frame(), data.frame(seurat_clusters = names(cosine_classes), cosine = cosine_classes, dge = dge_classes) )

  # 4. Add cell type information

  # from pre-existing data frame annotation:
  # This data frame can be generated with generateCellTypeAnnotation
  # Original cell type has a the cluster N in front:
  cell_type_noN = lapply(cell_type_ann$cell_type2, function(x){ a<-str_split(x," ", simplify = T);a[-1]   })
  lapply(cell_type_noN, paste,collapse = " ") %>% unlist() -> cell_type_noN
  # new cell type field with no number
  cell_type_ann$cell_type = cell_type_noN
  # We are going to left join the data matrix with the cell types
  ann_data %>% left_join(cell_type_ann %>% select(seurat_clusters, cell_type),by='seurat_clusters') -> ann_data


  # 5. Plot the annotated heatmap
  ann_data_pheatmap = ann_data %>% select(seurat_clusters, cosine)
  row.names(ann_data_pheatmap) <- ann_data_pheatmap$seurat_clusters
  ann_data_pheatmap$cosine = ann_data_pheatmap$cosine

  bmp.mat.raw.fil %>% pheatmap(clustering_distance_rows = dist.cosine(bmp.mat.raw.fil),cutree_rows = k,
                      annotation_row = ann_data_pheatmap %>% select(cosine) )

  # Return
  # Expression matrix as data frame with meta data included in the last columns
  # cluster, cosine_class, dge_class (k as defined by hyperparameter)
  return(ann_data)
}

# Generate avg.matrix as an annotated data frame with cell type information as additional column
annAvgMatrix<-function(which.genes = seven.receptors , seurat_obj = tiss.norm, cell_type_ann = cell_type_ann_tm){
  # Make simple cell type col in data frame (no cluster N)
  cell_type_noN = lapply(cell_type_ann$cell_type2, function(x){ a<-str_split(x," ", simplify = T);a[-1]   })
  lapply(cell_type_noN, paste,collapse = " ") %>% unlist() -> cell_type_noN
  # new cell type field with no number
  cell_type_ann$cell_type = cell_type_noN

  # calculate avg Matrix from seurat_ob
  x = avg.matrix(seurat_obj, which.genes,by="seurat_clusters")

  as.data.frame(x) %>% mutate(seurat_clusters = row.names(x)) %>%
  	left_join(cell_type_ann %>%
  	select(seurat_clusters, cell_type, tissue), by="seurat_clusters") -> x.ann

  return(x.ann)


}

plotDGE<-function(res.list_ = list(), which.genes  = seven.receptors,   min_exp = 2,   k = 8){

  # We compute both matrices
  # A differential expression matrix (i,j = N of differentially expressed genes between cell type i and j, scaled by fold change)
  # The average expression matrix
  bmp.dge.mat = distMatrix(res.list_, all_clusters, which.genes = which.genes, 10^-6, subset_genes = T)
  bmp.mat.raw = avg.matrix(tiss.norm,which.genes, by='seurat_clusters')

  # Filter for cell types with enough expression
  x_fil <- rowSums(bmp.mat.raw)>=min_exp
  # Filter in DGE matrix
  bmp.dge.mat[names(which(x_fil)),names(which(x_fil))] %>% pheatmap(cutree_rows  = k, border_color = NA, clustering_method = 'complete', col = viridis(100))
}

#DGE matrix

tryFindMarkers<-function(i,j,which.pathway){
        result = tryCatch( {
            FindMarkers(tiss.norm,ident.1 = all_clusters[i],ident.2 = all_clusters[j],features = which.pathway,logfc.threshold = 0.2,verbose = F)
        }, warning = function(w){
            print('warning catched')
            return(0)
        }, error = function(e){
            print('error catched')
            return(0)

        }, finally ={

        } )
    return(result)
}

# with min.pct parameter. Should substitute the above version
tryFindMarkers2 <-function(i,j,which.pathway, min.pct.expr = 0 ){
    result = tryCatch( {
        FindMarkers(tiss.norm,ident.1 = all_clusters[i],ident.2 = all_clusters[j],min.pct = min.pct.expr, features = which.pathway,logfc.threshold = 0.2,verbose = F)
    }, warning = function(w){
        print('warning catched')
        return(0)
    }, error = function(e){
        print('error catched')
        return(0)

    }, finally ={

    } )
    return(result)
}


# Update August 10th 2020
# DEsingle
# Calculate differential expression using DEsingle package

# Note:  if using other datasets: check for the cellID colname
tryFindMarkers_DEsingle<-function(brain_lineage  = c() , cluster1 = 0, cluster2 = 1, BPPARAMs = c() ){
  brain_lineage@meta.data %>% filter(seurat_clusters %in% c(cluster1,cluster2)) %>% select(ID, seurat_clusters) -> which_cells
  res_desingle = DEsingle(brain_lineage[['RNA']]@counts[, which_cells$ID ], group = droplevels(which_cells$seurat_clusters), parallel = T, BPPARAM= BPPARAMs )

}

# Calculate distance between cell types based on their DEsingle DGE profiles

# Aug 13th 2020 OK
# test.clusters includes the list of seurat_clusters as they were tested pairwise by DEsingle
# This version assumes that the DEsingle data.frame output includes all genes in the pathway list
cosineDGE<-function(res.list.desingle =list(), p_value_low = 10^-2 ,genes = pathways[[1]], i, j , test.clusters, method = 'cosine' ){
  # For a pair of clusters we define a distance metric that relies on the DGE test
  # It is basically a weighted cosine distance
  # The weight array theta is binary, with theta =1 when p_value is significant and 0 otherwise
  # Therefore, we calculate the cosine distance taking into account only those genes that pass the p_value threshold
  ii = min(which(test.clusters==i),which(test.clusters==j))
  jj = max(which(test.clusters==i),which(test.clusters==j))
  res.list.desingle[[ii]][[jj]] -> x
  # if directly from DEsingle output (all genes included)
  x %>% mutate(gene = row.names(x)) %>% filter(gene %in% genes) %>% select(gene, norm_total_mean_1, norm_total_mean_2, norm_foldChange,pvalue.adj.FDR) -> x
  row.names(x)<-x$gene
  # When filtering using %in% the order is changed to that of the genes array, but this shouldn't be a problem
  # sort in the original order
  #x<-x[genes,]

  # if a gene is missing from the data frame, that row will become NA if we subset using the full gene list
  # FIX THIS
  genes = x$gene # from now on keep only those genes that actually appear in the data frame

  # make mask vector: theta = 1 iff p_value < min_threshold
  theta = rep(0, length(genes))
  names(theta) <- genes
  # make mask based on the p_values
  x %>% filter(pvalue.adj.FDR<=p_value_low) -> DEgenes
  theta[DEgenes$gene] = 1
  # norm_total_mean in DEsingle is not log normalized.
  dist_array = rbind(log(x$norm_total_mean_1+1) * theta, log(x$norm_total_mean_2+1) * theta)

  # if one or more genes did not appear at all in the output of DEsingle, they woudl show NA values
  # we can replace them by zeros such that the distance metric ignores them (both cosine and euclidean)
  # The genes are in order so, the replacement should be trivial
  dist_array[is.na(dist_array)] = 0

  if(method=='cosine'){
    d_i_j =  dist.cosine( dist_array)
    if(is.na(d_i_j)){d_i_j = 0}
  }else if(method=='euclidean'){
    d_i_j =  dist( dist_array)
  }else if(method == 'fold'){
    # Sep 9th 2020
    # fold change is the ratio of un-normalized data
    # by taking the log of the ratio we make it symmetric for up and downregulation
    # Absolute value makes the metric more like a distance
    #d_i_j = sum(abs(log(x$norm_foldChange) * theta))
    d_i_j = sum(abs(dist_array[1,] - dist_array[2,]))

  }


  if(sum(theta)==0){
      return(0)
    }else{
      return(d_i_j)
    }
}

# Aug 13th 2020 build a distnace matrix from pairwise DGE comparisons using DEsingle.
# OK
distMatDEsingle<-function(res.list = list() , p_val = 10^-2, which.genes = pathways[[1]],
                          cluster_list = c(), which.method = 'cosine', log_transform = F){


  dist_mat = matrix(0,length(cluster_list),length(cluster_list));
  for(i in 1:(length(cluster_list)-1)){
      for(j in (i+1):length(cluster_list)){
          dist_mat[i,j] = cosineDGE(res.list,p_value_low = p_val, genes = which.genes, cluster_list[i],cluster_list[j], cluster_list, method = which.method)
      }
  }


  x = dist_mat + t(dist_mat)

  if(log_transform){
    x = log10(x+1)
  }


  colnames(x)<-cluster_list
  rownames(x)<-cluster_list
  return(x)
}


# Estimate number of clusters from the output of silhPathwayControl
# Using the silhouette score, we can estimate the approximate number of clusters for real pathways
# and for random pathways
numberClustersPathways<-function(silh_cutree, pathway_list){


    p = 8 # which pathway

    n_clust_rand = data.frame()

    for(p in 1:length(silh_cutree)){
      x = silh_cuttree[[p]][[1]]
      fil_window = 7 #length of smoothing filter

      f11 = rep(1/fil_window,fil_window)

      n_clusters = c()
      for(i in 1:dim(x)[1]){

          fil_x = stats::filter(x[i,], f11, sides = 2)  #apply filter to silh series

          n_clusters[i] = which.max(fil_x) + 1 # which k has maximum silh score
      }

      # compute confidence interval for the mean
      a = mean(n_clusters)
      s = sd(n_clusters)

      n = length(n_clusters)
      # Assumes normal distribution (is not!)
      error <- qnorm(0.975)*s/sqrt(n)
      left <- round(a-error, 2)
      right <-round(a+error, 2)


      # plot histogram with CI in the title
      # x11()
      # hist(n_clusters %>% log() , main =paste(all_pathways_names[p], '95% CI (',left,',',right ,')'))

      # make data frame
      n_clust_rand = rbind(n_clust_rand, data.frame(n_genes = length(pathway_list[[p]]) , mean = a, ci_low =left,ci_high =right))

    }


    return(n_clust_rand)
}
