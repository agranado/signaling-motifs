# Global clustering of transcriptomes and pathway states
# Jan 18 2021
# We have 1.3k leiden clusters from multiple datasets
# From there we create a global UMAP using Seurat by treating the leiden clusters as single cells
# We seek to cluster then the pathway states to find out recurrent expression across cell types


# General pipeline
makeMainDataFrame <-function(this_pathway, umap_coords = F){
  # We fetch the gene data from the seurat_obj
  # And merge with meta data
  pathway_matrix<-FetchData(master_seurat, this_pathway) # Log norm
  devel_adult <- cbind(pathway_matrix, master_seurat@meta.data %>% dplyr::select(Tissue, age,dataset,cell_ontology_class, Cell_class))
  row.names(devel_adult)<-1:dim(devel_adult)[1]


  devel_adult$global_cluster = 1:dim(devel_adult)[1] # here we haven't filtered anything. Comes directly from Seurat obj
  return(devel_adult)
}

normalizedDevel <- function(this_pathway, sat_val =0.99, fill_zero_rows = F ){
    devel_adult <- makeMainDataFrame(this_pathway) #pink variables go to Shiny

    devel_adult %>% mutate(cell_id = paste(global_cluster, dataset,sep="_")) -> devel_adult

    x =devel_adult[,this_pathway]
    max_sat_gene = apply(x, 2, quantile, sat_val) # starts from x


		# are there any 0 sat values ?
		if(sum(max_sat_gene==0)>0){
			max_val_gene = apply(x, 2, max) # starts from x
			max_sat_gene[max_sat_gene==0] <- max_val_gene[max_sat_gene==0]
			}

    for(s in 1:dim(x)[2])
        x[which(x[,s]>max_sat_gene[s]),s]<- max_sat_gene[s]

    x<- min.maxNorm(x)
		if(fill_zero_rows)
			x[x==0] = 10^-7

		devel_adult[,this_pathway] <- x

		row.names(devel_adult) <- devel_adult$global_cluster
    return(devel_adult)
}
# Make clustered data frame after 1st round of the Spectral pipeline

makeMasterClustered <- function(p_list, k_opt = 40 ){
		norm_counts <- p_list$counts %>% as.data.frame()
		# save row.names as global cluster id
		rownames_to_column(norm_counts, 'global_cluster') -> norm_counts

		#Run the plotMotif3D to extract the 3D UMAP coordinates
		scatter_export = plotMotif3D(p_clust = p_list$heatmap ,
																	which_motif = 1,
																	scatter.data,
																	k_opt = k_opt,
																	export_csv = T)

		# Let's create a single master data frame with motif labels
    # here we add the rownames from the count matrix
		scatter_export$global_cluster <- as.character(norm_counts$global_cluster)

		# this will be already filtered
		left_join(norm_counts, scatter_export, by ='global_cluster') %>%
					dplyr::filter(motif_label>0) -> master_clustered

		return(master_clustered)

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
		min_max = T,
    sum_norm = F,
    sat_quantile = F,
    sat_val = 0.99,
    min_expr_single_cell = 0){

  #filter more one gene: filter cell that only express 1 gene
  #gene_min_expr: min value to consider a gene as ON, 0 by default means no filtering

  # sum_norm: normalize each profile by the sum of the pathway genes.
  # min_max: min max scaling. don't use together with sum_norm

	# 2. Basic clustering with cosine distance
	# Returns count matrix and meta data
	# we could remove this to speed up things..
	p = plotAllProfiles(devel_adult = devel_adult, which_pathway = which_pathway, max_satu = max_satu,
											filter_unique = unique_cell_types, min_expr = min_expr,
											make_heatmap = F, norm_by_sum = sum_norm,
                      saturate_quantile = sat_quantile, quant_threshold = sat_val,
                      gene_min_expr = 0 )

	# return objects
	count_mat =p[[2]] # count matrix directly from devel_adult no filter or transform
	fil_rows = p[[3]] # plotAllProfiles applies the filters
	meta_data = p[[4]] # plotAllProfiles generates meta data from Seurat

	# 3. MinMax normalization before Spectral clustering
	if(min_max)
		count_mat <- min.maxNorm(count_mat)

  # 3.1 After min.max normalization, we can filter cell types that express more than one gene
  # This will prevent random pathways from creating clusters that only express one gene  (high silhouette high diversity)
  # fil_rows is a vector of row names to be included in the downstream analysis
  # if we want to further filter the matrix, we need to do an intersect with the second names vector
  fil_rows = intersect(fil_rows , row.names(count_mat)[which(rowSums(count_mat > min_expr_single_cell ) >1)] )  # at least two genes expressed

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
  # For the heatmap we saturate the signal:
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
      # k is the number of eigen-vectors
    spec_res = Spectrum(t(count_mat), maxk = k, showres = F,
                        method = spec_method,
                        kerneltype = spec_kernel)

    # 5. Manual spectral clustering (on k)
    manual_res = manualSpectral2(A = spec_res$similarity_matrix, labels = meta_data, k_opt =k)
    # Plot heatmap with spectral matrix
    spectral_dist = manual_res[[3]]

    # 6. Cluster eigen-vector and create dendrogram
    # ward.D2 seems to be a good option to cluster the eigen-matrix
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
                          filter_type = 1, norm_by_sum = F,
                          saturate_quantile = F, quant_threshold = 0.99,
                          gene_min_expr = 0.2){

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
    # Second layer of the filter makes sure that cells are expressing at least n genes
    # with a min value indicated in the parameters
    if(filter_type ==1){
      # For removing cell types with no expression OR cell types that only express one gene
  	   fil_rows = rowSums(x)>min_expr & (rowSums(x > gene_min_expr  ) >1)
		     fil_rows<-row.names(x)[fil_rows]
		# Normalize by cell-type (use max value for a receptor as reference)
    }else if(filter_type ==2){
    # Filter
    }
    # Jan 2021: Normalize by total pathway activity (experimental)
    # x_fil is the matrix used by the heatmap
    # Here we normalize by max gene in pathway or by sum of pathway genes
    # we apply the normalization to all profiles. But we do filter by min.expresison when running the heatmaps
    # we also return the fil row list of fitered cell types
    if(norm_by_sum){
      x_fil = apply(x, 1, function(x_r){x_r/sum(x_r)}  )
      colnames(x_fil) <- row.names(x) # apply returns transposed matrix
      x_fil = t(x_fil) # should work
    }else{
      x_fil = x # do nothing to the data and return the same matrix
    }

    # For min max normalization, outliers can bias the scaling
    # Here we saturate the counts to the 0.99 percentile such that outliers above this value do not affect the max
    # Each feature is saturated by their corresponding quantile

    if(saturate_quantile){
      max_sat_gene = apply(x, 2, quantile, quant_threshold) # starts from x
      for(s in 1:dim(x)[2])
        x[which(x[,s]>max_sat_gene[s]),s]<- max_sat_gene[s]

      x_fil = x  # this is the matrix we are using now.
      # be careful with applying multiple normalizations at the same time (not recommended)
      # since the order might not be as expected by default
      # Here I assume only this normalization will be used (followed by min.max downstream)
    }

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

pca_dist<- function(master_seurat, devel_adult = data.frame() ){

  dist_pca = dist(Embeddings(master_seurat ,'pca')) %>% as.matrix()


  row.names(dist_pca)<-devel_adult$global_cluster
  colnames(dist_pca)<-devel_adult$global_cluster

  return(dist_pca)
}


# Diversity score based on the UMAP coordinates
makeUmapStats <- function(scatter_export2 = data.frame(),
                          silh_scores = data.frame(),
                          dist_method = 'umap', user_dist_matrix = matrix() , this_pathway = c() ){


	# calculate silhouette score
	if(length(silh_scores ) ==0)
		rand_silh_scores  = cluster_silhouette(scatter_export2, this_pathway)

	#compute the distance in the umap space
  if(dist_method =='umap'){
     # make a distance matrix in the UMAP space
     umap_mat <- scatter_export2 %>% dplyr::select(UMAP_1, UMAP_2, UMAP_3)
     row.names(umap_mat) <- scatter_export2$global_cluster

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
				pull(global_cluster) %>% as.character() # we are using the global_cluster id as row.name

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
recursiveSpectral <- function(p_list = list(), k_opt = 30,n_motifs = 30,
                              silh_cutoff = 0.5,
                              master_clustered = data.frame(), n_iter = 9,
                              spectrum_kernel = 'density', this_pathway = c()){

      # 6. Compile results from recursive clustering
    recursive_res <- clusterRecursive(p_list, k_opt, min_s_default =silh_cutoff,
                                        master_clustered = master_clustered,
                                        max_iter = n_iter,  spec_kernel = spectrum_kernel, this_pathway = this_pathway )
    all_clusters_recursive = recursive_res$recursive_labels
    remaining_profiles = recursive_res$remaining

    # 7. New data frame with unique recursive labels
    recursive_df_list<-makeMasterRecursive(all_clusters_recursive, master_clustered)
    master_recursive = recursive_df_list$master
    master_recursive_tidy = recursive_df_list$tidy
    label_map = recursive_df_list$label_map #data frame with global_cluster, motif_label, pathway_clust

    # if there are less labels that the proposed number of final motifs,
    # we overwrite the number of clusters and define n_motifs = n_recursive_labels
    if(n_motifs>label_map$motif_label %>% unique() %>% length)
       n_motifs = label_map$motif_label %>% unique() %>% length

    final_motifs_res = makeFinalMotifs(master_recursive_tidy, n_motifs = n_motifs, label_map)
    scatter_export2 = final_motifs_res$scatter_export
    recursive_master_labels = final_motifs_res$recursive_label
    motif_matrix = final_motifs_res$motif_matrix # Expression matrix of average expressing within a recursive cluster
    final_motifs_df = final_motifs_res$final_motifs # map with ID for recursive clustering and final motifs ID
    p_motifs = final_motifs_res$p_motifs # heatmap object with the dendrogram corresponding to the final motif clustering (merging of recursive clusters)

    # 8. Final data.frame
    master_clustered %>%
        rename(motif_label_original = motif_label) %>%
        left_join(recursive_master_labels  %>% dplyr::select(global_cluster,merged_label), by ='global_cluster') %>%
        rename(motif_label = merged_label) -> df


    df<-df[!is.na(df$motif_label), ]
    # df now contains ALL data for those datapoint that were classified during the recursive algorithm
    # UMAP, profiles, meta_data and cluster label

    return(list(master = df, scatter_export = scatter_export2,
                final_motifs = final_motifs_df, tree = p_motifs, recursive_tidy = master_recursive_tidy))
}



## #  #  #
##  #  #  #
#### Recursive parameters

# runs as a script for now
# calls plotMotif3D which has the meta data
clusterRecursive <- function(p_list, k_opt, min_s_default =0.5,
                            master_clustered = data.frame(),
                            max_iter =  9, spec_kernel = 'density', this_pathway = c() ){
  # 1. Run recursive clustering based on the 1st round of spectral
  # this can start right away after the pipeline
  # prepare parameters for the long run
  aa = letters[seq( from = 1, to = 26 )]
  aa = c(aa, paste(aa, aa, sep=""), paste(aa, aa,aa, sep=""), paste(aa, aa,aa,aa,  sep=""))
  round_ids = aa

  #min_s_default = 0.5
  # from the first round (we use 0.5 to make it more stringent)

  # This section runs right aways after the original pipeline

  # 1.1 make master
  #master_clustered <- makeMasterClustered(p_list, k_opt = k_opt)

  # 1.2 compute the silhouette scores for the original clustering
  silh_res_sum <- cluster_silhouette(master_clustered, this_pathway)

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
  i = 2
  min_s = 0.3
  cluster_factor = 10
  #for(i in 2:max_iter){
  while(dim(remaining_profiles)[1]>cluster_factor){
  	# 1. Re-cluster profiles that did not pass the threshold in the previous iteration
  	# for i =1 the previous iteration is the main pipeline.
    # calculate k based on the number of remaining profiles:
    if(dim(remaining_profiles)[1] > 20){
      round_k = round(dim(remaining_profiles)[1]/cluster_factor) # aprox clusters of size  10 profiles
    }else{
      round_k = 2
    }

  	round2_res = silh_recluster(remaining_profiles, k_r = round_k, spec_kernel = spec_kernel)
  	# save results and new data fram with motif label
    # silhouette score is the average for that label
  	silh_res2 = round2_res$silh
    # Data.frame of individual profiles with their label for the current iteration
  	round2 = round2_res$labels
    # return individual silhouette score for each data point
    raw_silh = round2_res$raw_silh
  	# We take out the profiles that passed the score threshold
    silh_res2 %>% dplyr::filter(sm> min_s) %>% pull(label) %>% as.character( ) -> which_labels_pass
  	# create motif_label column for compatibility with the 1st round of clustering
    if(length(which_labels_pass)==0){
      # When we get a round with no high quality clusters we reduce the threshold for silhouette score and the cluster factor for next round
      min_s = 0.2
      cluster_factor = 5
      silh_res2 %>% dplyr::filter(sm> min_s) %>% pull(label) %>% as.character( ) -> which_labels_pass

      if(length(which_labels_pass)==0){ # no high quality clusters even at silh =0.2 then break
        print("no high quaility clusters found in this iteration. Finishing recursive clustering. ")
        remaining_profiles = round2 # no filtering of labels. Last round. Exit
        break

      }else{
        # reset the silhouette score for next round
        min_s = 0.3
      }
    }

  	round2 %>%
  		dplyr::filter( pathway_clust %in% which_labels_pass) %>%
  		select(global_cluster, pathway_clust) %>%
  		mutate(motif_label = paste(pathway_clust,round_ids[i], sep =""))-> profiles_out
  	# profiles_out is an annotated data frame with global_cluster and motif_label
  	# which is ALL we care about (for v1.0)

  	# Make heatmap to visualize the quality for clusters during this iteration
  	makeHeatmap(round2, which_labels_pass, this_pathway,
  							title = paste('Clustering round', toString(i)) )

  	# Save the profiles that passed the threshold
  	all_clusters_recursive[[i]]<-profiles_out
  	# save the profiles that didn't pass the threshold and iterate again
  	remaining_profiles = round2 %>% dplyr::filter(!pathway_clust %in% which_labels_pass)
  	print( paste('Run', i, 'finished. Profiles out: ', toString(dim(profiles_out)[1])))
    # print data frame for those profiles that passed the silhouette threshold
    print(  silh_res2 %>% dplyr::filter(sm> min_s))
    print(  raw_silh %>% dplyr::filter(label %in% which_labels_pass) %>% arrange(label, desc(silh))  ) # print scores for individual data.points

    i = i+1
  } # For loop

  dev.off()

  return( list(recursive_labels = all_clusters_recursive, remaining= remaining_profiles) )
} # clusterRecursive()

# Silhoutte score for recursive iterations

# This function receives a dataframe that was already filtered for low silhouette scores
# It will re-cluster the gene expression profiles (which are included in the data.frame)
# Main KEY = global_cluster

# Output: 2 dataframes
# Silh_res: results from the silhouette score for each cluster
# round2: metadata annotated with the new cluster label as pathway_clust

silh_recluster<-function(not_pass_clustered = data.frame() , k_r = 20,
                        spec_method = 1, spec_kernel = 'density'){
  # subset the clusters
  count_mat_re = not_pass_clustered[,this_pathway] %>% as.matrix()
  row.names(count_mat_re) <- not_pass_clustered$global_cluster
  # keep row names as global clusters
  meta_data_re = not_pass_clustered
  row.names(meta_data_re) <- meta_data_re$global_cluster

  # iterate the clustering pipeline
  # Here we can specify all parameters for Spectrum
  # and for the custom manualSpectral function
  round2 <- reClusterPathway(count_mat = count_mat_re,
                            meta_data = meta_data_re,
                            k = k_r,
                            spec_method = spec_method, spec_kernel = spec_kernel
                            )

  # create labels for round 2
  #round2 %>% dplyr::mutate(motif_label = paste(pathway_clust,'b', sep = "")) -> round2

  # Silhouette score requires numeric labels
  #dist_mat <- dist.cosine(round2[,this_pathway] %>% as.matrix) %>% as.matrix
  dist_mat <- dist(round2[,this_pathway] %>% as.matrix) %>% as.matrix #euclidean
  round2_labels <- round2$pathway_clust %>% as.numeric
  names(round2_labels) <- round2$global_cluster

  s2 = silhouette(round2_labels, dist_mat)
  silh_res2 = data.frame(label = s2[,1], silh = s2[,3] , neighbor = s2[,2])

  # Here we exclude data points that don't pass a min silhouette score threshold
  # They will be passed to the next round of clustering
  # The average silh for good cluster should be higher after removing those individual points
  #silh_res2$label[silh_res2$silh<0.2] = -1

  silh_res2 %>% dplyr::group_by(label) %>%
                dplyr::summarise(sm = mean(silh), n_data = n() ) %>%
                arrange(desc(sm)) -> silh_res_summary

  return(list( silh = silh_res_summary, labels = round2, raw_silh = silh_res2) )
}

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
               clustering_distance_rows = dist(as.matrix(x)),  # dist.cosine(as.matrix(x)),
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

      return(list(scatter_export = scatter_export2, recursive_label = recursive_master_labels,
                      motif_matrix = motif_matrix, final_motifs = final_motifs_df, p_motifs = p_motifs))
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

# Calculate silhouette score for a given clustering
cluster_silhouette <-function(df = data.frame() , this_pathway = c() , dist = 'euclidean', return_singles = F ){
    #assumes that the data frame contains the pathway profiles in wide format
    ## the input data.frame must contain global_cluster AND motif_label
    ## motif_label must be numeric friendly
		df %>% dplyr::select(c(this_pathway, global_cluster, motif_label)) -> df_mat

		# silhouette requires numeric cluster labels
		labs <- df_mat$motif_label %>% as.numeric()
		names(labs) <- df_mat$global_cluster

		x <- df_mat[,this_pathway] %>% as.matrix
		row.names(x) <- df_mat$global_cluster
		if(dist =='cosine'){
			s<-silhouette(labs, dist.cosine(x))
		}else if(dist=='euclidean'){
			s<-silhouette(labs, dist(x))
		}

		s_df <- data.frame(motif_label = s[,1], silh = s[,3])
		s_df_mean <- s_df %>% dplyr::group_by(motif_label) %>%
							dplyr::summarise(ms = mean(silh)) %>%
							arrange(desc(ms)) %>% as.data.frame() %>%
							left_join(df %>% dplyr::group_by(motif_label) %>% count, by="motif_label")
    if(!return_singles){
		  return(s_df_mean)
    }else{
      df$silh = s_df$silh
      return(df)
    }
}

# df now contains ALL data for those datapoint that were classified during the recursive algorithm
# UMAP, profiles, meta_data and cluster label
compareSilhouette<- function(master_clustered, df, this_pathway = c()){
    print("Silhouette scores calculated")
    # 6. Compare the silhouette scores for Spectral (default) vs recursive
    # silhoutte scores for the spectral and recursive methods
    s1<-cluster_silhouette(master_clustered, this_pathway)
    s2<-cluster_silhouette(df, this_pathway, this_pathway)

    plot(s1$ms, ylim = c(0,1),type = "o")
    lines(s2$ms, type = "o", col = "red")
    # we use the recursive clustering silhouette scores
    # scatter_expor2 has the recursive labels for each profile
    umap_stats <- makeUmapStats(scatter_export2 , s2)
    return(umap_stats)
}

# Plot Individual heatmaps
# Subset by dataset and make heatmaps for a given pathway
makeSeparateHeatmaps<-function(devel_adult, this_pathway){
  colors$dataset %>% unique -> data_classes
  for(i in 1:length(data_classes)){
    which_datasets = which(colors$dataset==data_classes[i])
    slice_indexes <- meta_master %>% dplyr::filter(dataset %in% which_datasets) %>% dplyr::pull(global_cluster)
    pheatmap(devel_adult[slice_indexes,this_pathway],
        annotation_row = devel_adult %>% dplyr::select(Tissue),
        annotation_colors = colors
      )
  }
}



# # # # # #
# # # # # #
# # # # # #
# Hierarchical clusterin on recursive sub-clusters

# Two functions to compute the silhouette score for a given proposed clustering
# optimal k for recursive sub-clusters using the average profile for each label
silhoutte_hclust_optimalk <- function(recursive_res, metric = "euclidean", clust_method ="complete", return_mat = F, max_y = 0.4,
																		min_y = 0.1, cluster_range = 2:40){
		# after running the pipeline
		# Make average matrix across recursive labels
		recursive_res$recursive_tidy  %>%
		    group_by(motif_label.y, gene) %>%
		    summarise(m_expr = mean(expr)) %>%
		    spread(key =gene, value = m_expr) %>%
		    as.data.frame -> motif_matrix

		final_motifs_df <- recursive_res$final_motifs
		p_motifs <- recursive_res$tree

		x<-motif_matrix %>% tibble::column_to_rownames("motif_label.y")
		row.names(x) <- motif_matrix$motif_label.y


		if(!return_mat){
			# inspect the number of clusters based on the recursive labels

      if(metric =="euclidean"){
        dist_mat = dist(x)
      }else{
        dist_mat = dist.cosine(as.matrix(x))
      }

        tree = hclust(dist_mat , method = clust_method )

			sapply(cluster_range, function(i){silhouette(cutree(tree, i), dist_mat )[,3] %>% mean() } ) %>%
	plot(type ="o", main= "Silhouette for sub-clusters", ylab ="silhouette",
				xlab ="N clusters", ylim =c(min_y,max_y))
			abline(h = 0.2);abline(h =0.3)
		}else{
			return(x)
		}
}

# optimal k for recursive sub-clusters but using the individual profiles  + recursive labels
silhouetteRecursive_celltypes <-function(recursive_res, this_pathway, singles = F, cluster_range = 2:40, metric = "euclidean", clust_method = "complete", filter_profiles = F, min_expr = 0.3){

		hms = c()
    s_dist = list()
		for(i in cluster_range ){
		    master_hclust<- recluster_Hclust(recursive_res ,
														data.frame()  , this_pathway,
														n_motifs = i, clust_method = clust_method, metric = metric , filter_low = filter_profiles , fil_threshold = min_expr )

		    s_hclust <- cluster_silhouette(df = master_hclust, this_pathway, return_singles = singles, dist = metric)

        if(!singles){
		        hms[i] = s_hclust$ms %>% median() # this is the mean of the mean silhouette

        }else{
            s_dist[[i]] = s_hclust$silh
        }
    }

    if(!singles){
		    return(hms)
    }else{
        #names(s_dist) <- as.character(cluster_range)
        return(s_dist)
    }

}

# Re-cluster the recursive sub-clusters
# assings labels to individual data.points

recluster_Hclust <- function(recursive_res,master_recursive, this_pathway, n_motifs,clust_method ="complete" , metric = 'euclidean', filter_low = F, fil_threshold = 0.3){
				# Make tidy data frame so we can average across clusters
				#gather(df, c(this_pathway), key = "gene", value ="expr") %>%
				#    filter(motif_label>0) %>% group_by(motif_label, gene) %>%
				#    summarise(mean_expr = mean(expr)) %>%
				#    spread(gene, mean_expr ) -> tidy_clustering

				# after running the pipeline
				# Make average matrix across recursive labels
				recursive_res$recursive_tidy  %>%
				    group_by(motif_label.y, gene) %>%
				    summarise(m_expr = mean(expr)) %>%
				    spread(key =gene, value = m_expr) %>%
				    as.data.frame -> motif_matrix

				x<-motif_matrix %>% tibble::column_to_rownames("motif_label.y")

        # Here we can filter recursive labels with no expression that might introduce noise in the clustering

        if(filter_low){

          fil_pass = names(which(rowSums(x) > fil_threshold ))
          # the names here are the recursive labels
          x <- x[fil_pass, ]

        }
        # Hierarchical clustering on the 1st round of classes

        if(metric=="euclidean"){
				      dist_mat = dist(x)
        }else{
              dist_mat = dist.cosine(x %>% as.matrix )
        }

				tree = hclust(dist_mat, method =clust_method)
				#n_motifs = 20 # based on the silhouette score
				clust_h = cutree(tree, n_motifs)
				hclust_labels = data.frame(hclust_final = clust_h , motif_label.y = names(clust_h))

				#hclust_labels$motif_label %>% as.character() -> hclust_labels$motif_label.y

				recursive_res$recursive_tidy  %>% spread(key =gene, value = expr) -> master_df

        # if we are filtering by recursive label, here we remove those labels that did not
        # pass the min expression threshold and continue to produce the final data.frame with the filtered profiles
        if(filter_low)
          master_df %>% dplyr::filter( motif_label.y %in% hclust_labels$motif_label.y ) -> master_df

			  master_df %>% left_join(hclust_labels, by ="motif_label.y") -> master_clustered_hclust
				# Make final labels with the id: motif_labels
				master_clustered_hclust %>%
					mutate(motif_label = hclust_final) -> master_clustered_hclust
        # the output data.frame contains individual profiles with their final labels after applying hclust on the recursive-labels
        # this master_clustered data frame should be used from now on for all other functions.
        # row.names won't match if compared to previous data structures
				return( master_clustered_hclust)
}


############################################################################################
############################################################################################

# # # # #
# Export for App
# 1.1 We need a new devel_adult object
exporToShiny<-function(p_list, devel_adult, scatter_export){
    # save the row.names as global_cluster id
    rownames_to_column(devel_adult, 'global_cluster') -> devel_adult_ann # this is not normalized
    norm_counts <- p_list$counts %>% as.data.frame()
    # save row.names as global cluster id
    rownames_to_column(norm_counts, 'global_cluster') -> norm_counts

    # 2. Run the plotMotif3D to extract the 3D UMAP coordinates
    #scatter_export = plotMotif3D(p_clust = p_list$heatmap , which_motif = 15, scatter.data,k_opt = k_opt, export_csv = T)
    # 3. Save data frames so the app can access them. They must match on dimensions and pathway profiles
    write.csv(scatter_export, file = "app/global_transcriptome_motifLabeled.csv", quote =F, row.names = T)
    # 4. Meta data
    write.csv(devel_adult_ann[p_list$fil_rows, ], file = "app/annotated_counts.csv", quote =F, row.names = T)
    # normalized data (by sum of pathway genes)
    write.csv(norm_counts[p_list$fil_rows,], file = "app/filtered_counts.csv", quote =F, row.names = T)
}



# Control
# Random sets of genes
# Main function:
# For a given number of clusters, it will run the Spectrum + manual_spectral pipeline
# Returns a data frame with cluster statistics + the count matrix for the specified set of genes
silh_run <- function(x){
    this_pathway = x
    devel_adult <- makeMainDataFrame(this_pathway) #pink variables go to Shiny
    k = 29# k+1 where 1 represents the outlier cluster from recursive
		#k =  70 for the first round

    out<-tryCatch(
        {
            p_list<-clusterPathway(
                devel_adult = devel_adult,
                which_pathway = this_pathway,
                max_satu =1, # not needed if using min.max
                min_expr = 2, #filte on rowsums
                k_opt = k, # number of clusters Spectral
                spec_method = 1,
                spec_kernel = 'stsc',
                unique_cell_types = F,
                min_max = T,
                sat_quantile = T, # saturate values at quantile
                sat_val = 0.99,
								min_expr_single_cell = 0.2# use this quantile as max for min.max
            )

            master_clustered = makeMasterClustered(p_list, k_opt = k)
            s1<-cluster_silhouette(master_clustered, this_pathway = this_pathway) # original clustering
						umap_stats <-    makeUmapStats(master_clustered, s1,
																					this_pathway = this_pathway,
																					dist_method ='user',
																					user_dist_matrix = dist_pca)

						umap_stats_cosine <-    makeUmapStats(master_clustered, s1,
																					this_pathway = this_pathway,
																					dist_method ='user',
																					user_dist_matrix = dist_cosine_pca)




            list(umap_stats, umap_stats_cosine, master_clustered)

        },
        error=function(cond) {
            message(cond)
            return(NA)
        },
        warning = function(cond){
            return(NULL)
        },
        finally  = {}
    )# end tryCatch

    return(out)
}


# Feb 10th
# Parse the output into data.frames:
# Remove errors (appear as NA in the results list)
filter_na_errors<-function(x, n_length =3)
    x[sapply(lapply(x, length), function(x) x==n_length)]

# execute after the large control batch run
makeControl_df_ <-function(x = list() , dist_index = 1 ){
	# this fun accepts the list generated by lapply(rand_pathways, silh_run)
	#x <- rand_silh_all_bmp

	cell_type_dist <- filter_na_errors(x)
	rand_umap_stats_df<-do.call(rbind, lapply(cell_type_dist, function(x){x[[dist_index]]}  ) )

	# Filter NA values for singletons
	rand_umap_stats_df %>% dplyr::filter(!is.na(umap_dist) & !is.na(umap_dist_sd)) -> rand_umap_stats_df

	return(rand_umap_stats_df)
}


cluster_prob <- function(x, ms=0, d=0, n = 0){
	sum(x$ms >=ms & x$umap_dist>=d & x$n >=n)/dim(x)[1]
}

# for each data point in our target data frame
# we compute the p-value based on the null distribution in null_df
# null_df was built using the random sets of genes
cluster_pvals <- function(df_stats, null_df){
	p_vals = c()
	for(i in 1:dim(df_stats)[1]){
		x = df_stats[i,]
		p_vals[i] = cluster_prob(null_df, x$ms, x$umap_dist, x$n)
	}
	df_stats$p_vals = p_vals

	return( df_stats)
}

# Feb 12th: Alternatively we can compute the ecdf (p_values distribution)
# for each individual pathway. This is probably a better comparison since
# we want to know whether out target pathway is statistically significant as a pathway
# as opposed to finding significance at the cluster level
makeControl_df <- function(x = list() , dist_index = 1){
  	cell_type_dist <- filter_na_errors(x)
    for(i in 1:length(cell_type_dist))
      cell_type_dist[[i]][[dist_index]]$batch = i
    # join by row all batches
    # but now each batch should have a unique id so we can group later
	  rand_umap_stats_df<-do.call(rbind, lapply(cell_type_dist, function(x){x[[dist_index]]}  ) )

    # Filter NA values for singletons
  	rand_umap_stats_df %>% dplyr::filter(!is.na(umap_dist) & !is.na(umap_dist_sd)) -> rand_umap_stats_df
    #rand_umap_stats_df$umap_dist[is.na(rand_umap_stats_df$umap_dist)] <- 0

	  return(rand_umap_stats_df)

}



runAllControls <- function(experiment = 'hvg'){


      # 0. Select the control
      if(experiment =='hvg'){
        control_experiment = rand_silh_hvg_bmp_28clusters
      }else{
        control_experiment = rand_silh_all_bmp_28clusters
      }

      # 1. Make df from the long list run
      control_28_clusters<-control_experiment %>% makeControl_df()
      # 1.1 we only consider clusters with positive silhouette score.
      # (not doing it can bias the p-value since we can get large clusters with low S-score but
      # the large size of the cluster can dominate the p-value)
      control_28_clusters<-control_28_clusters %>% dplyr::filter(ms >0)

      # 2. Compute the p-vals for the target pathway df
      final_stats_pathway<- cluster_pvals(umap_stats_recursive, control_28_clusters) %>% arrange(p_vals)
      # 3. Compute the p-vals for the null df (what does this actually mean)
      # This distribution considers a cluster across ALL samplings
      # it gives us the p-value for an individual cluster when compared to all possible
      # clusterings for samples of size n_genes
      control_pvals <-cluster_pvals(control_28_clusters, control_28_clusters) %>% arrange(p_vals)
      control_pvals <- control_pvals %>% mutate(adj_pval = p.adjust(p_vals, method ="BH"))

      final_stats_pathway <- final_stats_pathway %>% mutate(adj_pval = p.adjust(p_vals, method ="BH"))
      final_stats_pathway$batch = 0 # for consistency



      n_batches = length(unique(control_pvals$batch))
      dist_axis = seq(0.0001, 1, by =0.001)
      distr_mat = matrix(0, n_batches, length(dist_axis))

      for(i in 1:n_batches){
          control_pvals %>% dplyr::filter(batch == i) -> this_sample
          distr = ecdf( this_sample$adj_pval)
          distr_mat[i, ] = distr(dist_axis)

      }

      mean_control_ecdf <-data.frame( dist_axis = dist_axis, adj_pval = apply(distr_mat, 2, mean), batch = 0)




      ggplot(control_pvals , aes(adj_pval , group = batch)) + stat_ecdf(geom = "step", alpha = 0.01) +
        scale_x_continuous(trans='log10') + theme_minimal() + theme(text = element_text(size = 25)) +
        annotation_logticks(sides = 'b') + stat_ecdf(data = final_stats_pathway, aes(adj_pval), color = "hotpink3", size = 1.5) +
        geom_line(data = mean_control_ecdf, aes(x = dist_axis, y = adj_pval), size = 1.5 )  +
        ylab("Fraction of profiles") + xlab("Adj p-value") +
        ggtitle("Motif enrichment in pathway") + coord_cartesian(xlim=c(0.005,1))
}
