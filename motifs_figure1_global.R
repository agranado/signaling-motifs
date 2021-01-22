# Global clustering of transcriptomes and pathway states
# Jan 18 2021
# We have 1.3k leiden clusters from multiple datasets
# From there we create a global UMAP using Seurat by treating the leiden clusters as single cells
# We seek to cluster then the pathway states to find out recurrent expression across cell types


# From the global data frame: extract pathway counts and perform Spectral clustering
# it outputs the clustering by Spectral
clusterPathway<- function(
    devel_adult = data.frame() ,
		which_pathway = 'Notch',
		max_satu = 0.8,
		min_expr = 0.1,
		k_opt = 18, # number of clusters Spectral
		spec_method = 2,
		spec_kernel = 'stsc',
		unique_cell_types = F,
		min_max = T){


	# 2. Basic clustering with cosine distance
	# Returns count matrix and meta data
	# we could remove this to speed up things..
	p = plotAllProfiles(devel_adult = devel_adult, which_pathway = which_pathway, max_satu = max_satu,
											filter_unique = unique_cell_types, min_expr = min_expr,
											make_heatmap = F)

	# return objects
	count_mat =p[[2]] # count matrix directly from devel_adult no filter or transform
	fil_rows = p[[3]] # plotAllProfiles applies the filters
	meta_data = p[[4]] # plotAllProfiles generates meta data from Seurat

	# 3. MinMax normalization before Spectral clustering
	if(min_max)
		count_mat <- min.maxNorm(count_mat)

	# 4. Run spectrum distance matrix
	spec_res = Spectrum(t(count_mat[fil_rows,]), maxk = 40, showres = F,
												method = spec_method,
												kerneltype = spec_kernel)

	# 5. Manual spectral clustering
	manual_res = manualSpectral2(A = spec_res$similarity_matrix, labels = meta_data[fil_rows,], k_opt =k_opt)
	# Plot heatmap with spectral matrix
	spectral_dist = manual_res[[3]]

	count_mat_fil = count_mat
	count_mat_fil[count_mat_fil>max_satu] = max_satu

	# 6. Cluster eigen-vector and create dendrogram
	p_clust <- pheatmap(t(count_mat_fil[fil_rows,]),
											clustering_distance_cols = spectral_dist,
											clustering_method = 'ward.D2',
											silent = T)

	# 7. Label cells with pathway state
	meta_data$pathway_clust = rep(NA, dim(meta_data)[1])

	motif_labels<-cutree(p_clust$tree_col,k_opt)
	meta_data[names(motif_labels),]$pathway_clust = motif_labels
	meta_data$pathway_clust <- as.character(meta_data$pathway_clust)

	# Colors for pathway_clust labels
	colors$pathway_clust = makeQualitativePal(meta_data$pathway_clust %>% unique() %>% length() )
	names(colors$pathway_clust) <- meta_data$pathway_clust %>% unique() %>% as.numeric() %>% sort() %>% as.character()

	# 8. Final heatmap and save dendrogram
	pp = pheatmap(t(count_mat_fil[fil_rows,]),
									clustering_distance_cols = spectral_dist,
							    annotation_col = meta_data[fil_rows,] %>% dplyr::select( Tissue,pathway_clust,dataset),
									clustering_method = 'ward.D2',
									annotation_colors = colors,
									show_colnames = F, cluster_rows = F, fontsize = 12,
									cutree_cols = k_opt, silent = T , col = blues_pal(100))

	return(list(heatmap = pp, colors = colors, counts =  count_mat, fil_rows =  fil_rows, dist_mat  = spectral_dist))
}

# function for recursive clustering
# This version does not need all the pre-processing steps, since it is called after clusterPathway
# should be quicker Jan 22nd 2021
reClusterPathway<- function(count_mat = matrix() ,
                            meta_data = data.frame(),
                            k = 20,
                            spec_kernel = 'density',
                            spec_method =2
                           ){
    # Everything is fine as long as we keep the global cluster labels
    # We need:
    #  1. count matrix (already normalized)
    #  2. meta data with global cluster
    #  3. This will output the metadata with the cluster label

      # 1. Run spectrum distance matrix
    spec_res = Spectrum(t(count_mat), maxk = k, showres = F,
                        method = spec_method,
                        kerneltype = spec_kernel)

    # 5. Manual spectral clustering
    manual_res = manualSpectral2(A = spec_res$similarity_matrix, labels = meta_data, k_opt =k_opt)
    # Plot heatmap with spectral matrix
    spectral_dist = manual_res[[3]]

    # 6. Cluster eigen-vector and create dendrogram
    p_clust <- pheatmap(t(count_mat),
                        clustering_distance_cols = spectral_dist,
                        clustering_method = 'ward.D2',
                        silent = T)

    # 7. Label cells with pathway state
    meta_data$pathway_clust = rep(NA, dim(meta_data)[1])

    motif_labels<-cutree(p_clust$tree_col,k)
    meta_data[names(motif_labels),]$pathway_clust = motif_labels
    meta_data$pathway_clust <- as.character(meta_data$pathway_clust)

    return(meta_data)
}

# This function performs a basic clustering with cosine distance on the raw dataset
# It returns a list of filtered data points that can be further processed
# All Filtering should happen here
# It also return the count_mat that will be used by Spectrum so we can in principle do normalization here
plotAllProfiles<-function(devel_adult = data.frame(), which_pathway ='Bmp', max_satu = 0.5,
													dist_method ='cosine', min_expr = 0.2,
													k_cut = 20, filter_unique =T, make_heatmap  =F, filter_type = 1){

		my_colors_all<-colors
		# warning: devel_adult is global
		devel_adult %>% dplyr::select(Tissue, dataset, age, cell_ontology_class) -> devel_adult_meta
		# Rows are named in numerical order
		row.names(devel_adult_meta)<-row.names(devel_adult)
		row.names(devel_adult)<-row.names(devel_adult)
		# Gene selection
		#which_pathway = 'Bmp_Tgfb'
    #warning:  all_pathways is global
		pathway_genes <- all_pathways %>% dplyr::filter(pathway ==which_pathway) %>% pull(gene)
		x = devel_adult[,pathway_genes] %>% as.matrix()
		# Scaling and saturation
		x_fil = x
		row.names(x_fil) <- row.names(devel_adult)
		#max_satu = 0.5
		x_fil[x_fil>max_satu] = max_satu # For plotting saturation

    # FILTERING #####
    # Basic filtering by sum of all genes in the pathway
    # Filter non-expressing cell types
    if(filter_type ==1){
  	   fil_rows = rowSums(x)>min_expr # For removing cell types
		     fil_rows<-row.names(x)[fil_rows]
		# Normalize by cell-type (use max value for a receptor as reference)
    }else if(filter_type ==2){
    # Filter
    }
    # Jan 2021: Normalize by
    # x_fil is the matrix used by the heatmap
    # Here we normalize by max gene in pathway or by sum of pathway genes
    # we apply the normalization to all profiles. But we do filter by min.expresison when running the heatmaps
    # we also return the fil row list of fitered cell types
    x_fil = apply(x, 1, function(x_r){x_r/sum(x_r)}  )
    colnames(x_fil) <- row.names(x) # apply returns transposed matrix
    x_fil = t(x_fil) # should work
		#for(i in 1:dim(x)[1])
		#	x[i,] = x[i,]/max(x[i,]
		# Find unique cell types
		# Create unique id
		if(filter_unique){
				devel_adult_meta %>%
				mutate(unique_id = paste(Tissue, dataset,age, cell_ontology_class, sep="_")) -> devel_adult_meta
				# Get the first index of a cell type if it appears multiple times
				devel_adult_meta$unique_id -> vec
				tapply(seq_along(vec), vec, identity)[unique(vec)] -> aa
				lapply(aa, function(x){sample(x,1)}) %>% unlist() %>% as.character() -> unique_idx
				# Filter unique cell types that pass the min expression threshold
				fil_rows = intersect(fil_rows, unique_idx)
				#fil_rows = unique_idx
		}

		if(make_heatmap){

      if(dist_method=='cosine'){
        distance_metric = dist.cosine(x[fil_rows,])
      }else if(dist_method =='pearson'){
        distance_metric = as.dist(1-cor(t(x[fil_rows,]), method ='pearson'))
      }else if(dist_method =='spearman'){
        distance_metric = as.dist(1-cor(t(x[fil_rows,]), method ='spearman'))
      }

			p = pheatmap(t(x_fil[fil_rows,]),
			clustering_distance_cols = distance_metric,
			annotation_col = devel_adult_meta[fil_rows,] %>% dplyr::select( Tissue,dataset),
			annotation_colors = my_colors_all,
			show_colnames = F, cluster_rows = F, fontsize = 12,
			cutree_cols = k_cut, silent = T , col = blues_pal(100))
		}else{p = NULL}
		# Return notes 15012021:
		# x is the count matrix directly extracted from the normalized data
		# p in the pheatmap object that contains the dendrogram from which we can apply clustering
		# in practice this function only is used for visualization of a simple clustering on the data
		return(list(p, x_fil, fil_rows, devel_adult_meta))
}


# This function gets the UMAP coords from scratch
# We use this function to export an annotated data.frame for 3D plotting on the shiny app
# We need to make sure that the annotations correspond to gene expression and that the mouse event catches exactly the
# cel type that we are selecting.

plotMotif3D <- function(p_clust = p , which_motif = 1,
										scatter.data = data.frame(), k_opt = 20, export_csv = F ){
	# 3D UMAP
	# 1. Export from Seurat
	umap_coords<-Embeddings(master_seurat,reduction = 'umap') %>% as.data.frame()
	# cell_id column
	umap_coords %>% mutate(cell_id = row.names(umap_coords)) -> umap_coords

	# 2. Make scatter data with all UMAP info (here we can put pathway profiles later)
	left_join(umap_coords, meta_master, by = 'cell_id') -> scatter.data

	# Create clustering labels from Spectral dendrogram
	# after cutting, each leaf still has its original label
	motif_id<-cutree(p_clust$tree_col, k_opt)
	# Color by motif ID
	# Does not change the original scatter.data data.frame
	scatter.data$motif_id = 0
	# Note 01152021: every operation is done by name so we don't need to subset the original data.frame
	# just assigning data to the rows whose name appears in the tree
	# this means we can keep using the whole scatter.data
	scatter.data[names(which(motif_id==which_motif)),]$motif_id = 1

	if(!export_csv){
		# This function makes colors for dataset, Tissue, age (viridis ramp)
		colors<-makeColorsAll(meta_master, list() )
		# Interactive (does not work in Ubuntu for some reason)
		fig <- plot_ly(scatter.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
		marker = list(size = 3), color = ~motif_id,
		colors = c('grey','red'),
		text=~paste("Tissue:",Tissue,"<br>Age:",age,"<br>dataset:",dataset,"<br>Cell type:", cell_ontology_class),
		hoverinfo = 'text')
		fig %>% add_markers()
	}else{
		scatter.labeled = scatter.data
		scatter.labeled$motif_label = -1
		# we create new column: motif_label that will include all clusters
		# motif_id: we use this to highlight with color a single cluster
		scatter.labeled[names(motif_id), ]$motif_label = motif_id

		return(scatter.labeled)
	}
}


# Spectral clustering
# Custom method that takes the distance matrix output from Spectrum
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
