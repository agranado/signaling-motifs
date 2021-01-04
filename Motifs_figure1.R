
library(dplyr)
library(tidyr)

# Nov 2020 update
# Fig 1

# 1. General idea as shematics (illustrator)
#
# 2. Introduction of the dataset: heatmaps of expression
#       The motif idea is introduced as the clustering of the expression heatmap
#       This has to look great to the eye. The summary statistics will be presented in Fig 2
#       Number of expressed components for each pathway
# 3. Diversity scores for motifs
#       Global dendrogram colored by pathway motif
#       Boxplots of cell-type diversity for each pathway motif (ranked and colored mapped)
#





# Spectrum clustering

# Nov 9th 2020
# Functions to compute pathway-pathway correlations
# They all run from a list of .csv files with the average count matrices (normalized by house-keeping genes)
# The pipeline consists of:
# 1. Spectral clustering to determine the optimal number of clusters, the distance matrix and the cluster labels.
# 2. Confusion matrix between pairs of pathways
# 3. Adjusted Mutual Information for all pairs of pathways
#
# Control 1: A marker list can be used to define the transcriptome-pathway from which we can compute the correlation to the different patwhays
# Control 2: Srsf as an anti-motif that is not correlated to the transcriptome or to any individual pathway

# Load a pathway from the .csv exports from the jupyter notebook
# For each pathway: wide format with 3 annotation columns at the end
# Works for developmental datasets and the standard format from the jupyter notebooks export
load_path <- function(which_pathway = 'Bmp', which_datasets = c(), file_id  = '_raw_development.csv', meta_cols = 4,
    filter_pct = F, min.pct = 0.05){

    aa = read.csv(paste( which_pathway,file_id , sep=''))
		if(length(which_datasets)!=0){
				aa %>% filter(dataset %in% which_datasets) -> aa
		}
		#aa %>% select(-dataset) -> aa
    # Remove meta.data columns (not necessary right now)
    aa[,1:(dim(aa)[2]-meta_cols)] -> aa

    if(filter_pct){

        # filter_pct will read a file with the percetange of expressing cells per cluster per gene
        # The idea is to denoise the data by filtering minimum % of expression values
        # From a pct vs mean plot, we can see that around 10% of expression corresponds to 0.5 log expression which we can use as the low limit
        p = read.csv(paste( which_pathway,'_pct' ,file_id , sep=''))
        if(length(which_datasets)!=0){
            p %>% filter(dataset %in% which_datasets) -> p
        }
        # Same format as the countn matrix files
        # Remove meta.data columns (not necessary right now)
        p[,1:(dim(p)[2]-meta_cols)] -> p

        # Once we have both matrices, we can filter aa with the min.pas values
        nopass.min <-p<min.pct
        # for each column, we set all values that don't pass the threshold to zero
        # everything else stays the same
        for(i in 1:dim(aa)[2]){
            aa[ nopass.min[,i] ,i] = 0
        }

   }


    return(aa)
}


# Performs spectrum clustering an all pathways from the list
# It also plots the heatmap of the expression profiles using the Spectrum distance matrix for the heatmap dendrogram
# Returns a data-frame with columsn as pathways and values as the cluster-label for that pathway
# Rows are the cell types
# Spec method parameter is 2 by default, though some pathways behave better with method 1 (default parameter in Spectrum function)
# spec_methods array indicates, for each patwhay the preferred clustering method
clusterAllPathways <- function( test_pathways, file_id = '_raw_development.csv', save_plot = T,
                                plot_id = 'heatmap_Spectrum',   spec_methods_ = c(2,2,2,1,2,1,1,2,2,2,2), filter_data_string = 'Chan', cols_meta = 4,
                                use_manual_spectral = F, dist_methods = rep(0, length(test_pathways)),
                                quantile_norm = F, k_spectral = rep(10,length(test_pathways)),
                                return_list = F,saturate_val = 0.5, scale_data_type = 0,
                                plot_cosine = F, manual_2 = F, remove_outliers_list = rep(0, length(test_pathways)), outliers_list = c() , filter_min_vals = c(),
                                frac_out = rep(0, length(test_pathways))  , cluster_methods = rep('ward.D2', length(test_pathways)),
                                kernels = rep('stsc', length(test_pathways))){
  # return_list = F will return only the data.frame with the final labels for each cell type. DEFAULT
  # return_list = T will return the ann.data.frame + expression_mat + dist_mat for each pathway
  all_labels = list()
  # Read the first file to extract basic information
  aa = read.csv(paste( test_pathways[1], file_id , sep=''))

  # The dataFrame will include multiple datasets, we can filter by name to include only specific expression matrices
  # To include all, we can use word Tabula, which will include all adult datasets
  datasets<- aa$dataset %>% unique() %>% as.character()
  chan <-datasets %>% as.character() %>% grepl(pattern = filter_data_string)
  filter_datasets <- datasets[chan]


  all_dist_matrices = list() # For each pathway save the final distance matrix to return as a list
  all_path_mat  = list() # For each pathway return the list of expression matrices

  if(save_plot)
    pdf(paste('plots/',plot_id, '.pdf',sep=''), width = 10, height =4)

  for( i in 1:length(test_pathways)){

      # 1. Data pre-processing

      # 1.1 Load data (before, this was done inside spectrumClustering, but it did not make too much sense. So I moved it outside)
      # Each pathway gets an optimal number of clusters
      # WE need to compute this k from somewhere else (Spectrum maybe?)
      path_mat <- load_path(test_pathways[i], filter_datasets, file_id , cols_meta) %>% as.matrix()
      # Spectrum needs row names. These names come from Seurat. They mean nothing and always have an order. Since we have not filtered the data at all yet.
      row.names(path_mat) <- 1:dim(path_mat)[1]

      # 1.2 Filter genes that are not expressed on any cell type
      # Under some clustering, some genes might have 0 expression so let's filter for that first
      path_mat = path_mat[,which(colSums(path_mat)>0)]
      path_mat_outliers = path_mat # the original path_mat will be overwritten after applyig spectral clustering, so we save this one with all data, including outliers.

      # 1.3 Filter non-expresing cell types as outliers
      # if remove_min is >0 then we use remove_min as the minimum value of the sum of expression across all pathway genes in a cell type
      # if remove_min == 0, we don't filter for minimum expression
      remove_min = ifelse(filter_min_vals[i] >0,  filter_min_vals[i], 0)

      # 1.4 Are we removing outliers
      # check if there are outliers to remove
      # Alternatively, in the future there will be a function to automatically detect outliers
      remove_outliers = remove_outliers_list[i]
      # If outliers were provided manually we just take those. Otherwise, we compute them using the Mahalanobis distance metric
      if(remove_outliers)
        outliers = outliers_list[[i]]
        # Here, outliers could be empty, if that's the case. We are doing automated outliers detection inside spectrumClustering



      # We feed all the parameters above as input for Spectral clustering

      # 2. Spectrum
      # this function call the loading function
      # And clusters the ACM from there. So this is where data is actually loaded:
      # We need Spectrum for the distance matrix
      # REMOVING OUTLIERS TAKES PLACE HERE

      # frac_out is used to compute the Mahalanobis distance. It is the fraction of the data that will be used to compute the centroid from which outliers will be defined
      # it takes values from 0.75 (more outliers) to 0.95 (less outliers)
      list_res = spectrumClustering(path_mat ,spec_methods_[i], scale_data_type, remove_outliers, outliers, remove_min, frac_out[i],
                                    spec_kernel = kernels[i])
      labels = list_res[[2]]
      path_mat = list_res[[1]] # This matrix has not outliers
      dist_mat = as.dist(list_res[[3]]) # Matrix directly from Spectrum (converted to distance) No outliers (if remove_outliers  == T)

      # For manual spectral clustering we get the spectrum results object
      spec_res = list_res[[4]]
      outliers = list_res[[5]] # in case outliers were detected inside the function

      num_spec = length(unique(labels$label))

      # 3. Manual calculation of Graph Laplacian and clustering
      # Dec 1st:  This update applies Laplacian directly on the matrix that Spectrum returns
      # the only parameter is the number of eigen-vectors to consider:
      # we can use the old parameter k_spectral as the number of eigen-vectors to use
      # Which is highly related to the number of clusters we expect
      # manual2 note:
      if(manual_2){ # Now we compute all manual stuff here. This will be replaced by dist_methods[i]
        res_manual = manualSpectral2(A = spec_res$similarity_matrix, labels = labels, k_opt = k_spectral[i])
        # Compare the two clusterings using the count matrix
        # Distance comes from similarity matrix computed manually Q_t
        labels = res_manual[[1]] # updated data.frame with both spectrum and manual labels

        # manual2 note: Now this is dist(Y), where Y is the eigen matrix from the manual method
        Q_t = res_manual[[3]] # eigen-matrix
        num_clust = k_spectral[i] # specified by the number of eigen-vectors to consider
        dist_mat = dist(Q_t) # cluster on eigen-matrix
        clustering_method_heatmap = cluster_methods[i]


     }else{ # Default: Spectrum
        clustering_method_heatmap ='complete'
        # keep default behaviour
        all_labels[[i]] = labels$label
        num_clust = length(unique(labels$label))
     }



     # Spectrum somehow screwes the random seed for the random number generator
     # We can not use random number to set the seed because the random.seed is locked!
     # This sets a new random seed using clock data (?) which is independent of the curret random seed
     set.seed(NULL)
     set.seed(round(runif(1,0,1)*65000))

    # For plotting the heatmap and make comparisons
    # Depending on the spectralClsutering function, this matrix could be normalized (min.max or scaled)
    if (plot_cosine )
      dist_mat = dist.cosine(path_mat)



    # 4. MANUAL SPECTRAL LABELS
    # Nov 19th
    # Before making the actual heatmap, we use pheatmap for h-clustering using the manual distance matrix  + ward.D2
    p = pheatmap(path_mat, silent=T, clustering_distance_rows = dist_mat, clustering_method =clustering_method_heatmap )
    # now cut the tree to get a third label set and save it on data.frame labels
    labels$ward = treeClust(p, num_clust ) %>% as.character()
    # NOTE: OVER WRITE the FINAL LABEL with the ward.D2 clustering from the distance matrix
    all_labels[[i]] = labels$ward # assign final labels here

    # 5.  VISUALIZATION

     #dist_mat = dist.cosine(path_mat)
     # Make heatmap normalized by quantiles with a blue color scale for Fig 1
      if(save_plot){  # Save pdfs with clustering for each pathway
        blues_pal<-colorRampPalette(brewer.pal(n = 9, name = 'BuPu'))
        # For visualization, we apply saturation such that 50% expression of Gapdh is the highest value
        # This prevent visual artifacts where highly expressed genes bias the color palette range
        path_mat[path_mat>saturate_val] = saturate_val # 50% of Gapdh since this is normalized data


        # COLORS
        # Make colors for all labels
        # Make palette of distinct colors
        cols_manual  = makeQualitativePal(length(labels$manual %>% unique()))
        names(cols_manual) <- labels$manual %>% unique() %>% sort()
        cols_spectrum = makeQualitativePal(length(labels$label %>% unique() ))
        names(cols_spectrum) <- labels$label %>% unique()
        # Since we are going to map this to the global clusters, let's choose colors in deterministic order
        tail_colors = ifelse(i==2, T, F) # For Notch
        skip_colors = ifelse(i==3, 47,0) # For Wnt

        cols_ward = makeQualitativePal(length(labels$ward %>% unique() ), rand_order = F, tail_colors = tail_colors, skip = skip_colors)
        names(cols_ward) <- labels$ward %>% unique() %>%  sort()

        col_annotation = list(manual = cols_manual, label = cols_spectrum, ward = cols_ward)


        # Rotate heatmap for final figure in the paper (only plot ward partitioning of Q_dist)
        pheatmap(t(path_mat), annotation_col = labels %>% dplyr::select(ward), # choose here which labels to plot
                 clustering_distance_cols = dist_mat,
                 cutree_cols = num_clust, show_colnames = T, fontsize =15,
                 clustering_method = clustering_method_heatmap,
                  annotation_colors = col_annotation,border_color      = NA,
                  drop_levels       = TRUE, cluster_rows      = F,
                  color = blues_pal(20))
        # Plot the distance matrix for inspection
        plot_d = as.matrix(dist_mat)
        plot_d = -(plot_d-max(max(plot_d)))
        pheatmap(plot_d, annotation_row = labels %>% dplyr::select(ward),

                 cutree_rows = num_clust, show_colnames = T, fontsize =5,
                 clustering_method = clustering_method_heatmap,
                 clustering_distance_rows = dist_mat, clustering_distance_cols = dist_mat,
                  annotation_colors = col_annotation,border_color      = NA,
                  drop_levels       = TRUE,
                  color = blues_pal(20))



    }

    # Dec 1 2020
    # I we filtered the matrix for outlieres, for completeness, we add them to the list of labels with a labels of 0
    if(remove_outliers){
      labels_outliers = data.frame(ward = rep('S', length(outliers)))
      row.names(labels_outliers) = outliers
      labels_ = labels %>% dplyr::select(ward) # From the original labels data.frame we select only the label we are going to use as motif label
      labels_ = rbind(labels_ , labels_outliers)
      labels_ = labels_[as.character(1:dim(labels_)[1]), ] # we re-order using character indexing
      # Re assign the all_labels list. NOTE: this list has no row.names, it is assumed to be in oder already.
      all_labels[[i]] = labels_ # we sort the labels again in the right order: outliers where they are supposed to be
      # WE are going to concatenate these labels across pathways. We need the order to be the same.

      # make a heatmap with the outliers
      if(save_plot)
        pheatmap(t(path_mat_outliers[outliers, ]),color = blues_pal(20), cluster_rows = F, fontsize = 15)

    }


     k_manual = length(unique(labels$manual))
     print(paste('Pathway ', test_pathways[i], 'k_spec = ', toString(num_spec), ' k_manual = ', toString(k_manual) ))

     all_dist_matrices[[i]] = dist_mat # For each pathway save the final distance matrix to return as a list
     all_path_mat[[i]] = path_mat

  }

  if(save_plot)
    dev.off()

  motif_labels = do.call(cbind,all_labels)
  colnames(motif_labels) <- test_pathways

  return(motif_labels)
}

# Main function:
# Perform spectral clustering using Spectrum and the average count matrix
# Returns a list with the count matrix, the cluster labels and the distance matrix from Spectrum
# filter_datsets is the list of dataset names that we will consider for the analysis
# frac_out is the fraction of datapoints that we use to compute a centroid from which outliers will be compared (only used when removing outliers)
spectrumClustering <- function(path_mat = c(), spec_method = 1 , scale_type  = 0, remove_outliers = F, outliers = c() ,  remove_min = 0, frac_out = 0.80 , spec_kernel = 'stsc' ){

  # Dec 3rd: starts from path_mat object
  # This function will apply outlier detection and Spectral clustering on the filtered matrix
  # Returns the results from spectral clustering + the filtered matrix and cluster labels for each data point
  # This is the data matrix we feed into Spectrum
  if(scale_type ==0){
    x = t(path_mat)
  }else if(scale_type ==1){
    scaled_path_mat = scale(path_mat)
    x = t(scaled_path_mat)

  }else if(scale_type ==2 ){
    scaled_path_mat =min.maxNorm(path_mat)
    x = t(scaled_path_mat)
  }


  # Dec 1 2020 filter outliers at this point. We don't want to include them in the spectral clustering
  # The graph for spectral clustering is based on nearest neighbors, so outliers will still have neares neighbours, but that does not mean they are similar.
  # Dec 3: Outlier detection is not sensitive to data normalization or scaling. The same outliers are detected on scaled(x) or min.maxNorm(x)
  if(remove_outliers){
    if(length(outliers)==0) # if this array is empty, it means that it was not manually specified
      outliers = detectOutliers(path_mat, frac_data = frac_out) # then we do automated detection of outliers

    # For some pathways, cell types don't express enough of any individual gene.
    # We can label them as outliers from the beginning and not include them into the analysis
    if(remove_min>0){
      sum_a =apply(path_mat, 1, sum )
      #names(sum_a) <-row.names(path_mat) # we don't need to add names, because until here the data is still ordered numerically by row.names
      min_outliers<-which(sum_a<remove_min)
      outliers = c(outliers, min_outliers) %>% unique()
    }


    x = t(path_mat[-outliers, ]) # colnames are respected here
    path_mat = t(x) # we return the expression matrix after we removed the outliers. WE can then concatenate outliers outside this function
    # But the idea here is that this is the main function for spectral clustering and there if will only include the relevant profiles
  }

  # Return the scaled matrix? mostly for visualization. Since all calculations are done on the Q_T or dist_mat
  #if(scale_type >0) path_mat = scaled_path_mat

  # kernel = 'stsc' or 'density'
	res_spectrum  = Spectrum(x, maxk = 30, method = spec_method, showres = F, tunekernel = T, kerneltype =spec_kernel) # or kenerl = 'density'
	# Make data.frame with all labels
	labels = data.frame(label = res_spectrum$assignments %>% as.character())
	row.names(labels) <- colnames(x) # we assign the names of those data points that were actually processed

	# Take distance matrix from Spectrum results
  # Similarity matrix from Spectrum respects the row names form the original matrix
	dist_mat = max(max(res_spectrum$similarity_matrix )) - res_spectrum$similarity_matrix


  # Return looks fine
  # If data was filetered, the distance matrix has lower dimension than the original data
	return(list(path_mat, labels, dist_mat, res_spectrum, outliers))
}

test_pathways = c('Bmp', 'Notch', 'Wnt', 'Fgfr','Lpa', 'Srsf','Wnt_l','Bmp_l',
                  'Bmp_Tgfb','Eph_r','Eph_l')

# Dec 2nd 2020:  automated detection of outliers in data using mahalanobis distance and covariance
# WE compute the distance of all data point to the mean of the data but using only a fraction of the data,
# the most central fraction of data points. In this way we identify more outliers at the same level of significance.
detectOutliers <- function(path_mat = c(), sig_level = 0.0001, frac_data = 0.80){
    alpha <- sig_level
    cutoff <- (qchisq(p = 1 - alpha, df = ncol(path_mat)))
    output75 <- cov.mcd(path_mat, quantile.used = nrow(path_mat)* frac_data )
    mhmcd75 <- mahalanobis(path_mat, output75$center, output75$cov)
    names_outlier_MCD75 <- which(mhmcd75 > cutoff)
    return(names_outlier_MCD75)
}


# Manual clustering #############################################################################3
##################################################################################################
# # # # # # # # # # # # # # # # # # # # # ## # # # # # # #



# Starts from the similarity matrix output by Spectrum and performs Spectral clustering using the normalized graph Laplacian
# this is the output from the spectral clustering function
# list[[1]] = data.frame with annotation and labels
# list[[2]] = Spectrum object
# list[[3]] = data matrix
manualSpectral <- function(A=c(), labels =data.frame(),k_neigh = 10, n_graph_iter=5){
 	# A = spectral_hvg[[2]]$similarity_matrix
	# spectral_hvg[[1]] %>% select(spectrum) -> labels
	# 2. K-nearest similarity matrix: k = 10 as in the Spectrum paper
  #k_neigh  = 7
	for(i in 1:dim(A)[1]){
		A[i,] -> neighbor_distances
		neighbor_distances %>% sort(decreasing = T) -> sorted_neighbors
		neighbor_distances[neighbor_distances<sorted_neighbors[k_neigh]] = 0
		A[i,] = neighbor_distances /sum(neighbor_distances)
		}

	# 3. Graph difussion: 5 iterations
	Q_t = A # Initial condition
	I = diag(nrow =dim(A)[1], ncol = dim(A)[2]) # identity matrix
	for(p in 2:n_graph_iter){
		Q_t_minus = Q_t
		Q_t = A %*% Q_t_minus %*% t(A) + I
		}

	# 4. Make Degree matrix D: row sums of Q_t
	D = diag(rowSums(Q_t), nrow =dim(A)[1], ncol = dim(A)[2])

	# 5. Normalized Graph Laplacian (taken from the web)
	# [http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html](http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html)
	"%^%" <- function(M, power)
	with(eigen(M), vectors %*% (values^power * solve(vectors)))
	# Normalized Laplacian
	L <- (D %^% (-1/2)) %*% Q_t %*% (D %^% (-1/2)) # normalized Laplacian

	# 6. Eigenvectors and values from Graph Laplacian
	evL = eigen(L, symmetric = T)
	evL$values %>% diff() -> eigen_diff

	# 7. Identify optimal k of clusters based on eigen-gap
	# Consider only the k top eigenvectors

	k_opt = which.min(eigen_diff) + 1 # diff_ removes the first element
	X = evL$vectors[,1:k_opt]
	Y = matrix(0, dim(X)[1], dim(X)[2])
	# Normalize eigen vectors
	for(i in 1:dim(X)[1]){
		for(j in 1:dim(X)[2]){
			Y[i,j] = X[i,j] / sqrt( sum( X[i,]^2) )
		}
	}

	# 8. Gaussian mixture on eigenvector space
	gmm = GMM(Y, k_opt ,seed_mode = "random_subset", km_iter = 10,
	em_iter = 10, verbose = F)
	pr = predict_GMM(Y, gmm$centroids, gmm$covariance_matrices, gmm$weights)

	# 9. Assign labels
	labels$manual = pr$cluster_labels %>% as.character()
	dist_mat = as.dist( max(max(Q_t)) - Q_t )

	return(list(labels, dist_mat, Q_t) )
}

manualSpectral2 <-function(A = c(), labels = data.frame, min_dist_adj = 0 , norm_laplacian= T, k_gmm = 10, k_opt = 12, filter_Qt  =F ){


    Q_t = A
    # 4. Make Degree matrix D: row sums of Q_t
    D = diag(rowSums(Q_t), nrow =dim(A)[1], ncol = dim(A)[2])



      if(filter_Qt){
        Adj_mat = Q_t # from similarity matrix
        diag(Adj_mat)<-0 # Adjacency is zero at diagonal
        Adj_mat[Adj_mat<min_dist_adj] = 0 # binarize the matrix, filter for data.points with low similarity
        Adj_mat[Adj_mat>0] = 1 # Everything else becomes 1 to indicate a connection between those points
      }else{
        Adj_mat = Q_t # alternative, just take the similarity matrix directly with no filtering
        #diag(Adj_mat) = 0
      }

      # Degree matrix
      Deg_mat = diag(rowSums(Adj_mat),nrow = dim(Adj_mat)[1], ncol = dim(Adj_mat)[2]) # Degree matrix comes now naturally

      L = Deg_mat - Adj_mat # Unnormalized graph laplacian (most general form)

      if(norm_laplacian){

        # Here we can compute the normalized Laplacian
        # matrix power operator: computes M^power (M must be diagonalizable)
        "%^%" <- function(M, power)
            with(eigen(M), vectors %*% (values^power * solve(vectors)))

        L <- (Deg_mat %^% (-1/2)) %*% L %*% (Deg_mat %^% (-1/2))  # normalized Laplacian

      }


      evL <- eigen(L, symmetric=TRUE) # Eigen-decomposition


      evL$values %>% diff() -> eigen_diff #optimal N of clusters

      # 7. Identify optimal k of clusters based on eigen-gap
      # Consider only the k top eigenvectors

      if(k_opt ==0){
        k_opt = which.min(eigen_diff) + 1 # diff_ removes the first element
      }
      #X = evL$vectors[,1:10] # Select only top 10 eigen vectors
      # with the parameter k_opt we specify how many eigenvector we want to consider
      # the coice of eigenvectors will almost dictate the number of clusters
      X   <- evL$vectors[,(ncol(evL$vectors)-k_opt+1):ncol(evL$vectors)]


      # Eigen-matrix normalization
      Y = matrix(0, dim(X)[1], dim(X)[2])
      # Normalize eigen vectors
      for(i in 1:dim(X)[1]){
          for(j in 1:dim(X)[2]){
              Y[i,j] = X[i,j] / sqrt( sum( X[i,]^2) )
          }
      }


      # 8. Gaussian mixture on eigenvector space
      gmm = GMM(Y, k_gmm ,seed_mode = "random_subset", km_iter = 10,
                em_iter = 10, verbose = F)
      pr = predict_GMM(Y, gmm$centroids, gmm$covariance_matrices, gmm$weights)

      # 9. Assign labels
      # labels = data.frame(manual =  pr$cluster_labels %>% as.character())
      labels$manuel = pr$cluster_labels %>% as.character()
      row.names(labels) = row.names(A)
      row.names(Y) = row.names(A)




      return(list(labels, Y, dist(Y)))

}


treeClust<-function(p,k){
  #p2 = pheatmap(input_matrix, clustering_distance_rows = dist.cosine(input_mat), silent = T)
  aa = as.hclust(p$tree_row)
  return(cutree(aa, k))
}

# Wrapper de manualSpectral
# With new code for actually including the laplacian graph and eigen-vector decomposition
runSinglePathway<-function(pathway_file,min_sat = 4, k_neigh =10 , n_graph_iter = 5,
                           min_dist_adj = 0.05, filter_Qt = F , scale_data = 0,
                           spectrum_density = 'density', norm_laplacian = T, k_opt = 12, manual_processing =F,
                           k_gmm = 12, clust_method = 'complete', remove_outliers = F, outliers = c() ) {


    aa = read.csv(pathway_file, header = T)

    if(remove_outliers){
       # remove outliers that don't fit into a single cluster
       aa = aa[-outliers,]

    }

    labels = data.frame()
    bmp.mat<-aa[,1:(dim(aa)[2]-3)]
    row.names(bmp.mat) <-row.names(bmp.mat)

    if(scale_data==0){
      x = bmp.mat   %>% t()
    }else if(scale_data ==1){
      x = bmp.mat %>% scale()  %>% t()
    }else if(scale_data ==2){
      x = bmp.mat %>% min.maxNorm()  %>% t()
    }

    colnames(x)<-1:dim(x)[2]

    res.spec = Spectrum(x, maxk = 30, method = 2, tunekernel = T, showres = F, kerneltype = spectrum_density)

    A = res.spec$similarity_matrix

    # Spectrum returns a similarity matrix that can be further processed using their method description
    # However, this matrix could already be processed by spectrum
    if(manual_processing){
      # Keep only the top k-neighbours for each data point, everything else gets 0
      # Each row is nor normalized such that the total sum of similarities  to their neighbours is the same for each data point
      for(i in 1:dim(A)[1]){
          A[i,] -> neighbor_distances
          neighbor_distances %>% sort(decreasing = T) -> sorted_neighbors
          neighbor_distances[neighbor_distances<sorted_neighbors[k_neigh]] = 0
          A[i,] = neighbor_distances /sum(neighbor_distances)
      }

      # 3. Graph difussion: 5 iterations
      Q_t = A # Initial condition
      I = diag(nrow =dim(A)[1], ncol = dim(A)[2]) # identity matrix
      for(p in 2:n_graph_iter){
          Q_t_minus = Q_t
          Q_t = A %*% Q_t_minus %*% t(A) + I
      }

    }else{ # As per Spectrum documentation: the output similarity matrix is already the final similarity matrix Q_t
      # Nov 30 update: there is no need to process the matrix again!
      # This can be used as Q_t and follow from here with the spectral clustering algorithm
      Q_t = A

    } # Else apply Laplacian directly on the matrix they return

      # 4. Make Degree matrix D: row sums of Q_t
      D = diag(rowSums(Q_t), nrow =dim(A)[1], ncol = dim(A)[2])

      row.names(Q_t) = row.names(bmp.mat)
      colnames(Q_t) = row.names(bmp.mat)


      # Visual palette and color saturation
      blues_pal<-colorRampPalette(brewer.pal(n = 9, name = 'BuPu'))


    if(filter_Qt){
      Adj_mat = Q_t # from similarity matrix
      diag(Adj_mat)<-0 # Adjacency is zero at diagonal
      Adj_mat[Adj_mat<min_dist_adj] = 0 # binarize the matrix, filter for data.points with low similarity
      Adj_mat[Adj_mat>0] = 1 # Everything else becomes 1 to indicate a connection between those points
    }else{
      Adj_mat = Q_t # alternative, just take the similarity matrix directly with no filtering
      #diag(Adj_mat) = 0
    }

    # Degree matrix
    Deg_mat = diag(rowSums(Adj_mat),nrow = dim(Adj_mat)[1], ncol = dim(Adj_mat)[2]) # Degree matrix comes now naturally

    L = Deg_mat - Adj_mat # Unnormalized graph laplacian (most general form)

    if(norm_laplacian){

      # Here we can compute the normalized Laplacian
      # matrix power operator: computes M^power (M must be diagonalizable)
      "%^%" <- function(M, power)
          with(eigen(M), vectors %*% (values^power * solve(vectors)))

      L <- (Deg_mat %^% (-1/2)) %*% L %*% (Deg_mat %^% (-1/2))  # normalized Laplacian

    }


    evL <- eigen(L, symmetric=TRUE) # Eigen-decomposition


    evL$values %>% diff() -> eigen_diff #optimal N of clusters

    # 7. Identify optimal k of clusters based on eigen-gap
    # Consider only the k top eigenvectors

    if(k_opt ==0){
      k_opt = which.min(eigen_diff) + 1 # diff_ removes the first element
    }
    #X = evL$vectors[,1:k_opt] # Select only top 10 eigen vectors
    # when selecting the top k_opt (highest eigen-values) clustering is not good.
    # We want the vectors with the smalles eigen-values:
    X   <- evL$vectors[,(ncol(evL$vectors)-k_opt+1):ncol(evL$vectors)]



    Y = matrix(0, dim(X)[1], dim(X)[2])
    # Normalize eigen vectors
    for(i in 1:dim(X)[1]){
        for(j in 1:dim(X)[2]){
            Y[i,j] = X[i,j] / sqrt( sum( X[i,]^2) )
        }
    }


    # 8. Gaussian mixture on eigenvector space
    gmm = GMM(Y, k_gmm ,seed_mode = "random_subset", km_iter = 10,
              em_iter = 10, verbose = F)
    pr = predict_GMM(Y, gmm$centroids, gmm$covariance_matrices, gmm$weights)

    # 9. Assign labels
    labels = data.frame(manual =  pr$cluster_labels %>% as.character())
    row.names(labels) = row.names(bmp.mat)
    row.names(Y) = row.names(bmp.mat)
    pheatmap(Y, annotation_row = labels)

    #min_sat = 4
    x  = bmp.mat
    x[x>min_sat]=min_sat
    pheatmap(t(x), annotation_col = labels, clustering_distance_cols = dist(Y),
            col= blues_pal(20), cutree_cols = k_opt,
            cluster_rows  = F, drop_levels = T, border_color = NA, clustering_method =clust_method )
}

# Main spectral clustering of Tabula Muris 3m (Nov 24th 2020)
# Function to cluster markers
# Other related plots to Fig 1 (I need to update the location of these functions )


tf_list_file = 'Tabula_Senis_markers1kTF.csv'
marker_list_file = 'Tabula_Senis_markersHVG.csv'
# This function reads .csv files with markers for Tabula Muris
# Performs spectral clustering using spectrum (with a specified k)
# dist_methods 'cosine', 'euclidean', 'spectral'
# method = 1, 2, 3 for Spectrum
# filter_min: for a gene, sum over all cells.
clusterMarkersTabula3m <- function(markers_file = 'Tabula_Senis_markersHVG.csv',
                                   method = 1,
                                   meta_cols = 5, make_plot = T,
																	 merge_markers = F, pathway_ann = data.frame() ,
																	 scale_data = T, filter_min = 0.7, spec_method = 1,
																	 tree_clust_method ='ward.D2',
																	 manual_k = 0, dist_method = 'cosine', filter_id = c('3m')){
    # By default this function filter the 3m dataset
    # but we can cluster markers with all three (Senis + 3m) datasets
    if(!merge_markers){
        markers_df = read.csv(file =markers_file, header = T)
        markers_df %>% filter(dataset %in% filter_id) -> markers_df
        # From now on, everything is for Tabula Muris 3m mouse
    }else{
        # we take both datasets and make a single data.frame
        # 3000 HVG + 1k TF
        markers_df = read.csv(file =tf_list_file, header = T)
        markers_df %>% filter(dataset %in% filter_id) -> markers_df

        markers_df2 = read.csv(file = marker_list_file, header = T)
        markers_df2 %>% filter(dataset %in% filter_id) -> markers_df2
        # take the genes from the second data.frame and cbind
        markers_mat = markers_df2[,1:(dim(markers_df2)[2]-meta_cols)]
        markers_df = cbind(markers_mat, markers_df)
    }

    # Fix the capital Tissue field
    if(!'Tissue' %in% colnames(markers_df)){
      markers_df %>% rename(Tissue = tissue ) -> markers_df
    }


    # 0. Count matrix
    markers_mat = markers_df[,1:(dim(markers_df)[2]-meta_cols)]
    # Note: a number of highly variable genes show actually no expression except a random cell
    # which inflates the variance and contributes to clustering in a weird day
    # even after normalization by Gapdh, z-score = 15 which is outrageous
    which(colSums(markers_mat) >filter_min) -> filter_genes
    markers_mat = markers_mat[,filter_genes]

    markers_mat <-t(markers_mat)
    colnames(markers_mat) <- 1:dim(markers_mat)[2]

    #meta = markers_df %>% dplyr::select(Tissue, cell_ontology_class, cell_class, dataset)
    meta = markers_df %>% dplyr::select(Tissue, cell_ontology_class, cell_class, dataset, age)

		#pathway annotation
		if(length(pathway_ann)>0)
	    meta = cbind(meta, pathway_ann)

    # 1. Scale dataset
		if(scale_data){
	    data_mat = scale(t(markers_mat))
		}else{
			data_mat = t(markers_mat) # no scaling
		}

    row.names(data_mat) <- row.names(meta)

    # 2. Run Spectrum
    markers_spectral = Spectrum(markers_mat, method = spec_method, tunekernel = T, maxk = 40, showres = F)

    # 3. Add labels to meta data_frame
    meta$Tissue = meta$Tissue %>% as.character() # Otherwise factors screw everything
    # Add Spectral labels
    meta$spectrum = markers_spectral$assignments %>% as.character()

    # Prepare colors and make heatmap
    n_clust_markers =meta$spectrum %>% unique() %>% length()

    # Use similarity matrix from spectrum to plot the heatmap
    dist_mat = max(max(markers_spectral$similarity_matrix )) - markers_spectral$similarity_matrix

    meta$cell_class = as.character(meta$cell_class)
    # We have too many Fat Tissues, we can rename them
    meta$cell_class[meta$cell_class =='Stromal ']='Stromal' # correct this type I made when saving the .csv
    meta %>% mutate(tissue2 = ifelse(Tissue %in% c('GAT','MAT', 'BAT','SCAT'), 'Fat', Tissue)) -> meta2
    meta2$Tissue = meta2$tissue2
    row.names(meta2) <-meta2$seurat_clusters

    # cell class colors
    cell_class_cols = brewer.pal(length( meta2$cell_class %>% unique() ), 'Set3')
    names(cell_class_cols) <- meta2$cell_class %>% unique()

    # Add to the color list we had before
    my_tissue_colors$cell_class = cell_class_cols

		# PATHWAY MOTIF ANNOTATION Nov 20 2020
		# If we have pathway annotation, we select those for the row colors
    if(length(pathway_ann)>0){
			ann_heatmap = meta2 %>% dplyr::select( Bmp,Notch,Wnt,Tissue, cell_class)

			cols_bmp  = makeQualitativePal(length(pathway_ann$Bmp %>% unique()), rand_order = F)
      names(cols_bmp) <- pathway_ann$Bmp %>% unique() %>% sort()
			# we start the palette from the next available color
      cols_notch = makeQualitativePal(length(pathway_ann$Notch %>% unique() ),
									 rand_order = F, skip = 0, tail_colors = T)
      # notch colors are in reverse order so we have to move the S label to the beggining of the array
      # mind that this is a factor
      col_names = pathway_ann$Notch %>% unique() %>% sort()
      names(cols_notch) <- factor(c(col_names[length(col_names)] %>% as.character(), col_names[-length(col_names)] %>% as.character() ))

			cols_wnt = makeQualitativePal(length(pathway_ann$Wnt %>% unique() ),
									 rand_order = F, skip = 47, tail_colors = F)
      names(cols_wnt) <- pathway_ann$Wnt %>% unique() %>% sort()

      # Fix colors for outliers Dec 1st 2020
      # the label for outliers is S, so we can check here if there are outliers
      # We don't need to specify an outliers parameter, as long as we check the color names
      if(!is.na(cols_bmp['S'])) cols_bmp['S'] = '#FFFFFF'
      if(!is.na(cols_notch['S'])) cols_notch['S'] = '#FFFFFF'
      if(!is.na(cols_wnt['S'])) cols_wnt['S'] = '#FFFFFF'

			my_tissue_colors$Bmp = cols_bmp
			my_tissue_colors$Notch = cols_notch
			my_tissue_colors$Wnt = cols_wnt

    }else{
			ann_heatmap = meta2 %>% dplyr::select( Tissue, cell_class, age,dataset)
		}



		if(dist_method == 'cosine'){
			heatmap_dist_mat = dist.cosine(data_mat)
		}else if(dist_method =='spectral'){
			heatmap_dist_mat = as.dist(dist_mat )
		}else if(dist_method == 'euclidean'){
			heatmap_dist_mat = dist(data_mat)
		}


    # 4. Plot heatmap: scaled data, spectral dist matrix, ward.D2
    p = pheatmap(data_mat, annotation_row = ann_heatmap,
                 clustering_distance_rows = heatmap_dist_mat,
                 clustering_method = tree_clust_method, show_rownames = F,
                 show_colnames = F, cutree_rows = ifelse(manual_k>0, manual_k,n_clust_markers),
                 annotation_colors = my_tissue_colors,treeheight_col =0, silent = !make_plot)


		# Make aplot of distance boxplots in cell type space for each motif

    # Return annotated data.frame
    return(list(meta2 , markers_spectral, data_mat, p, my_tissue_colors))

}


# MOTIF diversity
# Once we found the pathway classes, we compute the diversity (in transcriptome) of the
# cell types using that particular motif
# We do that using euclidean distance on the 3k marker genes from Seurat
# Update Dec 1 2020: WE need to fix this function to process names of motifs that include letters and numbers
diversityDistance <- function(spectral_hvg, k = 15, scale_data = T){
	all_dists = data.frame()

	for(p in c('Bmp', 'Notch','Wnt')){
	    d_motifs = c()

	    this_path = p
			if(scale_data){
		    cell_type_mat = spectral_hvg[[3]] %>% scale()
			}else{
				cell_type_mat = spectral_hvg[[3]]
			}
      this_pathway_motifs = unique(spectral_hvg[[1]][,this_path])

	    for(i in 1:length(this_pathway_motifs)){
	        motif1_cell_types = cell_type_mat[spectral_hvg[[1]][,this_path]==this_pathway_motifs[i],]

          # Update Dec 3: Here I need a matrix in the pathway space to compute some kind of silhouette score or within-motif distance
	        dist_motifs1 = dist(motif1_cell_types)

	        d = as.matrix(dist_motifs1)
	        all_dists = rbind(all_dists, data.frame(pathway =this_path , motif = this_pathway_motifs[i] , dist_pairs =  d[upper.tri(d)], m = median( d[upper.tri(d)])))

	    }
	}
  # all_dists data frame now contains a tidy list of all pairs of datapoints grouped by motif
  # we can then calculate the distribution of pairwise distance within a motifs in the transcriptome space
  # IN principle, we can also calculate here the distance in the pathway space, as to create a plot of distance in pathway vs distance in transcriptome for the grouped data points


	# Calculate distance within cell classes
	markers_clust =cutree(spectral_hvg[[4]]$tree_row, k)
	markers_clust_df = data.frame(seurat_clusters = as.character(names(markers_clust)), cell_class = markers_clust  )
	control_dists = data.frame()
	for(i in 1:length(unique(markers_clust_df$cell_class))){
	    motif1_cell_types = cell_type_mat[markers_clust_df$cell_class==i,]
			if(sum(markers_clust_df$cell_class==i)>1){
				dist_motifs1 = dist(motif1_cell_types)
		    d = as.matrix(dist_motifs1)
				distr = d[upper.tri(d)]
			}else{
				d = 0
			}
	    control_dists = rbind(control_dists, data.frame(pathway= 'Cell_types',motif = i, dist_pairs =  distr , m = median( distr)))
	}

	#all_dists = rbind(all_dists, control_dists)
	# Compute the quantiles for the expectation from cell types
	control_dists$dist_pairs %>% summary() -> quantiles_cell_type

	return(list(all_dists, quantiles_cell_type))

}

# Make box-plots Nov  24th
# One for each motif colored by the same color-id from spectral clustering
# Draw lines on the quantile range for the expected distribution of related cell types
diveristyBoxPlot <-function(all_dists,this_path = 'Bmp', manual_colors = c() ){
    all_dists %>% filter(pathway == this_path) %>% ggplot(aes(x = reorder(motif, m,FUN = median), y =dist_pairs, fill = motif)) +
    geom_boxplot(outlier.shape = NA) + geom_hline(yintercept=quantiles_cell_type[2], linetype="dashed", color = "red") + geom_hline(yintercept=quantiles_cell_type[5], linetype="dashed", color = "red") + theme_classic()  +  theme(text = element_text(size = 20)) -> p1
    p1  = p1 + ylab('Distance between cell types (transcriptome)') + xlab('Pathway class') + ggtitle(this_path)

    if(length(manual_colors)>0){
        p1 = p1 + scale_fill_manual(values = manual_colors)
    }

    return(p1)
}



# Annotation functions
# Meta data for heatmaps
# Visualization colors


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
  cell_type_ann =  cell_type_top_label %>% dplyr::select(tissue, seurat_clusters, cell_type2, cell_ontology_class) %>% as.data.frame()
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
  cell_type_ann =  cell_type_top_label %>% dplyr::select(tissue, tissue_cluster, cell_type2, cell_ontology_class) %>% as.data.frame()
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
#
# meta_data = tiss.norm@meta.data
# heatmap_ann = generateCellTypeAnnotation(meta_data = meta_data)


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
