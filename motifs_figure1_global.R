# Global clustering of transcriptomes and pathway states
# Jan 18 2021
# We have 1.3k leiden clusters from multiple datasets
# From there we create a global UMAP using Seurat by treating the leiden clusters as single cells
# We seek to cluster then the pathway states to find out recurrent expression across cell types


# General pipeline
makeMainDataFrame <-function(this_pathway){
  # We fetch the gene data from the seurat_obj
  # And merge with meta data
  pathway_matrix<-FetchData(master_seurat, this_pathway) # Log norm
  devel_adult <- cbind(pathway_matrix, master_seurat@meta.data %>% dplyr::select(Tissue, age,dataset,cell_ontology_class))
  row.names(devel_adult)<-1:dim(devel_adult)[1]

  return(devel_adult)
}

# Clustering functions:

# From the global data frame: extract pathway counts and perform Spectral clustering
# it outputs the clustering by Spectral
clusterPathway<- function(
    devel_adult = data.frame() ,
		which_pathway = c('Bmpr1a', 'Bmpr2') , # list of genes
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
	spec_res = Spectrum(t(count_mat[fil_rows,]), maxk = k_opt, showres = F,
												method = spec_method,
												kerneltype = spec_kernel)

	# 5. Manual spectral clustering
	manual_res = manualSpectral2(A = spec_res$similarity_matrix,
                              labels = meta_data[fil_rows,],
                              k_opt =k_opt)
	# Plot heatmap with spectral matrix
	spectral_dist = manual_res[[3]] # eigen matrix

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
plotAllProfiles<-function(devel_adult = data.frame(),
                          which_pathway =c('Bmpr1a', 'Bmpr2') , max_satu = 0.5,
													dist_method ='cosine', min_expr = 0.2,
													k_cut = 20, filter_unique =T, make_heatmap  =F,
                          filter_type = 1){

		my_colors_all<-colors
		# warning: devel_adult is global
		devel_adult %>% dplyr::select(Tissue, dataset, age, cell_ontology_class) -> devel_adult_meta
		# Rows are named in numerical order
		row.names(devel_adult_meta)<-row.names(devel_adult)
		row.names(devel_adult)<-row.names(devel_adult)
		# Gene selection
		#pathway_genes <- all_pathways %>% dplyr::filter(pathway ==which_pathway) %>% pull(gene)
    pathway_genes = which_pathway # we change the function to get the list of genes Jan 25 2021
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

# Diversity score based on the UMAP coordinates
makeUmapStats <- function(scatter_export2 = data.frame(),
                          silh_scores = data.frame(),
                          dist_method = 'umap', user_dist_matrix = matrix() ){
	# make a distance matrix in the UMAP space
	umap_mat <- scatter_export2 %>% dplyr::select(UMAP_1, UMAP_2, UMAP_3)
	row.names(umap_mat) <- scatter_export2$global_cluster

	# calculate silhouette score
	if(length(data.frame() ) ==0)
		rand_silh_scores  = cluster_silhouette(scatter_export2)

	#compute the distance in the umap space
  if(dist_method =='umap'){
	   dist_umap <- dist(umap_mat) %>% as.matrix()
  }else if(dist_method =='user') {
      #compute distance on transcriptome markers (from Seurat n = 600)
      # here we could still use euclidean, cosine, etc
      dist_umap <- user_dist_matrix
  }


	m_umap_dist = c()
	sd_umap_dist = c()
	# -1 is the label for unclustered profiles
	n_motifs <- max(scatter_export2$motif_label %>% unique() )
	for(i in 1:n_motifs){

		# which global_clusters in this motif
		which_motif <- scatter_export2 %>%
				dplyr::filter(motif_label == i) %>%
				pull(global_cluster)

		# subset distance matrix
		d_motif<-dist_umap[which_motif, which_motif]
		# summary statistics
		m_umap_dist[i] <- mean(d_motif[upper.tri(d_motif)])
		sd_umap_dist[i] <- sd(d_motif[upper.tri(d_motif)])


	}

	umap_stats <- data.frame(motif_label =1:n_motifs, umap_dist = m_umap_dist, umap_dist_sd = sd_umap_dist)
	# join with silhouette scores
	silh_scores %>% left_join(umap_stats, by="motif_label") -> umap_stats

	return(umap_stats)

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
      # gmm = GMM(Y, k_gmm ,seed_mode = "random_subset", km_iter = 10,
      #           em_iter = 10, verbose = F)
      # pr = predict_GMM(Y, gmm$centroids, gmm$covariance_matrices, gmm$weights)
      #
      # # 9. Assign labels
      # # labels = data.frame(manual =  pr$cluster_labels %>% as.character())
      # labels$manuel = pr$cluster_labels %>% as.character()
      # row.names(labels) = row.names(A)
      row.names(Y) = row.names(A)




      #return(list(labels, Y, dist(Y)))
      return(list( data.frame() , Y,dist(Y) ) )

}


# Functions for recursive algorithm
# Make a hatmap of the selected profiles from the main data frame
makeHeatmap<-function(round2 = data.frame(),
											which_labels_pass= c() ,
											this_pathway = c()  ,
											title = ""){

	x<-round2 %>% dplyr::filter(pathway_clust %in% which_labels_pass) %>%
						dplyr::select(c(this_pathway))

	# Make heatmap
	pheatmap(x,
	         clustering_distance_rows = dist.cosine(x %>% as.matrix),
	         cutree_rows = length(which_labels_pass),
	         annotation_row = round2 %>% dplyr::select(pathway_clust),
	         color = magma(100),
	         show_rownames = F,
					 main = title,
					 cluster_cols = F,
					 border_color = NA)
}

###################
# Recursive pipeline
###################
# main wrapper of the recursive algorithm
# inputs are the 1st round from Spectrum + the meta data
recursiveSpectral <- function(p_list = list(), k_opt = 30,n_motifs = 30, silh_cutoff = 0.5, master_clustered = data.frame(), n_iter = 9 ){

      # 6. Compile results from recursive clustering
    recursive_res <- clusterRecursive(p_list, k_opt, min_s_default =silh_cutoff,master_clustered = master_clustered, max_iter = n_iter)
    all_clusters_recursive = recursive_res$recursive_labels
    remaining_profiles = recursive_res$remaining

    # 7. New data frame with unique recursive labels
    recursive_df_list<-makeMasterRecursive(all_clusters_recursive, master_clustered)
    master_recursive = recursive_df_list$master
    master_recursive_tidy = recursive_df_list$tidy
    label_map = recursive_df_list$label_map #data frame with global_cluster, motif_label, pathway_clust


    final_motifs_res = makeFinalMotifs(master_recursive_tidy, n_motifs = n_motifs, label_map)
    scatter_export2 = final_motifs_res$scatter_export
    recursive_master_labels = final_motifs_res$recursive_label

    # 8. Final data.frame
    master_clustered %>%
        rename(motif_label_original = motif_label) %>%
        left_join(recursive_master_labels  %>% dplyr::select(global_cluster,merged_label), by ='global_cluster') %>%
        rename(motif_label = merged_label) -> df


    df<-df[!is.na(df$motif_label), ]
    # df now contains ALL data for those datapoint that were classified during the recursive algorithm
    # UMAP, profiles, meta_data and cluster label

    return(list(master = df, scatter_export = scatter_export2))
}


# runs as a script for now
# calls plotMotif3D which has the meta data
clusterRecursive <- function(p_list, k_opt, min_s_default =0.5, master_clustered = data.frame(),max_iter =  9){
  # 1. Run recursive clustering based on the 1st round of spectral
  # this can start right away after the pipeline
  # prepare parameters for the long run
  round_ids<-c('a','b','c','d','e','f', 'g', 'h','i', 'j')
  round_k <- c(0, 35,30,30, 30,20,20, 20, 20, 20)# test params
  round_k <- c(0, 35,30,25, 25,20,20, 15, 12, 20)# original from Notion
  #round_k <- c(0, 35,30,30, 30,20,20, 20, 20, 20)# test


  min_s <- c(0, 0.5,0.5,0.5, 0.5,0.5,0.5, 0.5, 0.5, 0.5) # test
  min_s <- c(0, 0.4,0.4,0.4, 0.4,0.4,0.4, 0.5, 0.5, 0.5) # original from notion


  #min_s_default = 0.5
  # from the first round (we use 0.5 to make it more stringent)

  # This section runs right aways after the original pipeline

  # 1.1 make master
  #master_clustered <- makeMasterClustered(p_list, k_opt = k_opt)

  # 1.2 compute the silhouette scores for the original clustering
  silh_res_sum <- cluster_silhouette(master_clustered)

  silh_res_sum %>% dplyr::filter(ms<min_s_default) -> not_pass
  not_pass_clustered <- master_clustered %>%
  									dplyr::filter(motif_label %in% not_pass$motif_label)

  # Init master list
  # arrays and list to save outptus:
  # save the first batch of good clusters in the master list
  all_clusters_recursive <- list()
  all_clusters_recursive[[1]] <- master_clustered %>%
  																		dplyr::filter(!motif_label %in% not_pass$motif_label) %>%
  																		dplyr::select(global_cluster, motif_label) %>%
  																		dplyr::mutate(pathway_clust = 0)

  # Init workingn data frame
  remaining_profiles = not_pass_clustered
  pdf("recursive_clustering_heatmaps.pdf")
  # mainFor
  for(i in 2:max_iter){

  	# 1. Re-cluster profiles that did not pass the threshold in the previous iteration
  	# for i =1 the previous iteration is the main pipeline.
  	round2_res = silh_recluster(remaining_profiles, k_r = round_k[i])
  	# save results and new data fram with motif label
    # silhouette score is the average for that label
  	silh_res2 = round2_res$silh
    # Data.frame of individual profiles with their label for the current iteration
  	round2 = round2_res$labels
  	# We take out the profiles that passed the score threshold
    silh_res2 %>% dplyr::filter(sm> min_s[i]) %>% pull(label) %>% as.character( ) -> which_labels_pass
  	# create motif_label column for compatibility with the 1st round of clustering
  	round2 %>%
  		dplyr::filter( pathway_clust %in% which_labels_pass) %>%
  		select(global_cluster, pathway_clust) %>%
  		mutate(motif_label = paste(pathway_clust,round_ids[i], sep =""))-> profiles_out
  	# profiles_out is an annotated data frame with global_cluster and motif_label
  	# which is ALL we care about (for v1.0)

  	# Make heatmap to visualize the quality for clusters during this iteration
  	makeHeatmap(round2, which_labels_pass, this_pathway,
  							title = paste('Clustering round', round_ids[i]) )

  	# Save the profiles that passed the threshold
  	all_clusters_recursive[[i]]<-profiles_out
  	# save the profiles that didn't pass the threshold and iterate again
  	remaining_profiles = round2 %>% dplyr::filter(!pathway_clust %in% which_labels_pass)
  	print( paste('Run', i, 'finished. Profiles out: ', toString(dim(profiles_out)[1])))

  } # For loop

  dev.off()

  return( list(recursive_labels = all_clusters_recursive, remaining= remaining_profiles) )
} # clusterRecursive()

# after recursive clustering we build a data frame with the new labels
makeMasterRecursive <- function(all_clusters_recursive = list(), master_clustered = data.frame() ){
    # rbind the list of data.frames
    # they have unique global_cluster id
    do.call(rbind, all_clusters_recursive) -> aa
    # subset the master data.frame with those profiles that have a motif_label
    master_clustered %>% dplyr::filter(global_cluster %in% aa$global_cluster) -> master_recursive
    # join the motif_label
    master_recursive %>%
      left_join(aa %>% dplyr::select(global_cluster,motif_label), by ='global_cluster') -> master_recursive

    # let's make it tidy
    # which columns to exclude from gather
    meta_cols<-colnames(master_recursive)[which(!colnames(master_recursive) %in% this_pathway)]
    # gather
    master_recursive %>% gather(key ='gene', value ='expr', -meta_cols) -> master_recursive_tidy

    return(list( master = master_recursive, tidy = master_recursive_tidy, label_map = aa  ))
}
# average gene expression profiles for each recursive cluster
# some recursive clusters could be very similar
# Use input from above functions
makeFinalMotifs<-function(master_recursive_tidy = data.frame() , n_motifs =30, label_map = data.frame() ){
      # so the next step will be to apply hierarchical clustering on this matrix
      master_recursive_tidy %>%
      				group_by(motif_label.y, gene) %>%
      				summarise(m_expr = mean(expr)) %>%
      				spread(key =gene, value = m_expr) %>%
      				as.data.frame -> motif_matrix

      print("Main matrix generated")

      # 3. Heatmap of all profiles grouped by motif class
      # motif_matrix contains the gene expression profiles
      # and the motif_label in the format "10f"
      #n_motifs = 30  # how many final motifs do we want, this is important!
      #n_motifs = 35

      x<-motif_matrix[-1]
      row.names(x) <- motif_matrix$motif_label.y

      pdf("final_motifs.pdf", height = 18, width = 6 )
      p_motifs <- pheatmap(x[,this_pathway],
               clustering_distance_rows = dist.cosine(as.matrix(x)),
               color = magma(100),
               border_color = NA,
      				 cluster_cols = F,
      				 cutree_rows = n_motifs)

      dev.off()
      print("Main matrix saved as heatmap")


      # 4. Make the final scatter data.frame that will be used in the app
      # Here we select the final number of motifs
      final_motifs<-cutree(p_motifs$tree_row, k = n_motifs)
      # recursive is the label from each iteration
      # after clustering those by cosine distance we create the merged_label
      # which is going to be the final motif id
      final_motifs_df <- data.frame(motif_label= names(final_motifs), merged_label = final_motifs)

      # aa contains all the profiles that were clustered during the recursive algorithm
      # master_clustered has all profiles from the original round of Spectral
      # basically master_clustered has all profiles that passed the min.expr threshold
      # This new data.frame has final label as merged_label AND global cluster
      label_map %>% left_join(final_motifs_df,by='motif_label') -> recursive_master_labels

      # left_join with the scatter_export that goes into the app
      # For visualization
      # we need to over-write the column motif_label in the scatter_export df
      scatter_export %>% rename(motif_label_original = motif_label) %>%
      							left_join(recursive_master_labels %>%
      							dplyr::select(global_cluster,merged_label), by ='global_cluster') %>%
      							rename(motif_label = merged_label) -> scatter_export2



      # all profiles with no cluster get NA as value
      # for compatibility we assign -1 (used by the app)
      scatter_export2[is.na(scatter_export2$motif_label), ]$motif_label <- -1

      write.csv(scatter_export2, file = "app/global_transcriptome_motifLabeled.csv",
                  quote =F, row.names = T)


      print("Main UMAP coordinates labeled with the new clusters")

      return(list(scatter_export = scatter_export2, recursive_label = recursive_master_labels, motif_matrix = motif_matrix, final_motifs = final_motifs_df))
}
# scatter export has a the UMAP coordinates and now, the new  recursive labels

interactiveHeatmap <- function(master_recursive = data.frame(), this_pathway ){

    # let's make it tidy
    # which columns to exclude from gather
    meta_cols<-colnames(master_recursive)[which(!colnames(master_recursive) %in% this_pathway)]
    # gather
    master_recursive %>% gather(key ='gene', value ='expr', -meta_cols) -> master_recursive_tidy

    # so the next step will be to apply hierarchical clustering on this matrix
    master_recursive_tidy %>%
        group_by(motif_label, gene) %>%
        summarise(m_expr = mean(expr)) %>%
        spread(key =gene, value = m_expr) %>%
        as.data.frame -> motif_matrix

    x<-motif_matrix[-1]
    row.names(x) <-motif_matrix$motif_label
    return(x)
}

# df now contains ALL data for those datapoint that were classified during the recursive algorithm
# UMAP, profiles, meta_data and cluster label

print("Silhouette scores calculated")
# 6. Compare the silhouette scores for Spectral (default) vs recursive
# silhoutte scores for the spectral and recursive methods
s1<-cluster_silhouette(master_clustered)
s2<-cluster_silhouette(df)

plot(s1$ms, ylim = c(0,1),type = "o")
lines(s2$ms, type = "o", col = "red")
# we use the recursive clustering silhouette scores
# scatter_expor2 has the recursive labels for each profile
umap_stats <- makeUmapStats(scatter_export2 , s2)
