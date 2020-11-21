# Code for assessing differential gene expression in slingshot pseudotime
# This code uses the dataset from Zhong 2018 available as a Seurat object in ../dataset/Brain
# Aug 07 2020
#

# Pathway-Pathway correlations
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
load_path <- function(which_pathway = 'Bmp', which_datasets = c(), file_id  = '_raw_development.csv', meta_cols = 4){
    aa = read.csv(paste( which_pathway,file_id , sep=''))
		if(length(which_datasets)!=0){
				aa %>% filter(dataset %in% which_datasets) -> aa
		}
		#aa %>% select(-dataset) -> aa
    # Remove meta.data columns (not necessary right now)
    aa[,1:(dim(aa)[2]-meta_cols)] -> aa
    return(aa)
}

# Main function:
# Perform spectral clustering using Spectrum and the average count matrix
# Returns a list with the count matrix, the cluster labels and the distance matrix from Spectrum
# filter_datsets is the list of dataset names that we will consider for the analysis
spectrumClustering <- function(i =1 , test_pathways = c() , filter_datasets = c() , spec_methods = c() , files_id = '_raw_development.csv', n_meta = 4){
	# Each pathway gets an optimal number of clusters
	# WE need to compute this k from somewhere else (Spectrum maybe?)
	path_mat <- load_path(test_pathways[i], filter_datasets, files_id , n_meta) %>% as.matrix()
	# Spectrum needs row names
	row.names(path_mat) <- 1:dim(path_mat)[1]

	res_spectrum  = Spectrum(t(path_mat), maxk = 30, method = spec_methods[i], showres = F)
	# Make data.frame with all labels
	labels = data.frame(label = res_spectrum$assignments %>% as.character())
	row.names(labels) <- row.names(path_mat)

	# Take distance matrix from Spectrum results
	dist_mat = max(max(res_spectrum$similarity_matrix )) - res_spectrum$similarity_matrix

	num_clust =length(unique(labels))
	# Inspect with pheatmap (cor distance is not exactly the same as in Aff prop)

	return(list(path_mat, labels, dist_mat, res_spectrum))
}

test_pathways = c('Bmp', 'Notch', 'Wnt', 'Fgfr','Lpa', 'Srsf','Wnt_l','Bmp_l',
                  'Bmp_Tgfb','Eph_r','Eph_l')


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
                                return_list = F){
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
    pdf(paste('plots/',plot_id, '.pdf',sep=''))

  for( i in 1:length(test_pathways)){

      # this function call the loading function
      # And clusters the ACM from there. So this is where data is actually loaded:
      # We need Spectrum for the distance matrix
      list_res = spectrumClustering(i, test_pathways, filter_datasets, spec_methods_, file_id, cols_meta )
      labels = list_res[[2]]
      path_mat = list_res[[1]]
      dist_mat = as.dist(list_res[[3]])

      # For manual spectral clustering we get the spectrum results object
      spec_res = list_res[[4]]

      num_spec = length(unique(labels$label))


      # Get labels from manual spectral clustering:
      # Nov 12th: Add a new column with out manual spectral labels:
      # Makes a heatmap with the count matrix annotated by a few different labels
      # Manual clustering seems to give a better matrix that works well with ward.D2
      # looks good visually
	    res_manual = manualSpectral(A = spec_res$similarity_matrix,
	                            labels = labels, k_spectral[i])

	   # Compare the two clusterings using the count matrix
	   # Distance comes from similarity matrix computed manually Q_t
	   labels = res_manual[[1]] # updated data.frame with both spectrum and manual labels

     Q_t = res_manual[[3]]
     Q_dist = as.dist( max(max(Q_t)) -Q_t)

     # Set the clustering algorithm for the final labels:
     # Some pathways look better using the manual clustering and/or the Q_t distnace matrix
     if(dist_methods[i]){ # manual matrix
        clustering_method_heatmap = 'ward.D2'
        # Se dist_mat to the manual distance matrx
        dist_mat = Q_dist
        all_labels[[i]] = labels$manual

        num_clust = length(unique(labels$manual))

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



     #dist_mat = dist.cosine(path_mat)
     # Make heatmap normalized by quantiles with a blue color scale for Fig 1
      if(save_plot){  # Save pdfs with clustering for each pathway
        blues_pal<-colorRampPalette(brewer.pal(n = 9, name = 'BuPu'))
        # For visualization, we apply saturation such that 50% expression of Gapdh is the highest value
        # This prevent visual artifacts where highly expressed genes bias the color palette range
        path_mat[path_mat>0.5] = 0.5 # 50% of Gapdh since this is normalized data
        # Nov 19th
        # Before making the actual heatmap, we use pheatmap for h-clustering using the manual distance matrix  + ward.D2
        p = pheatmap(path_mat, silent=T, clustering_distance_rows = dist_mat, clustering_method =clustering_method_heatmap )
        # now cut the tree to get a third label set and save it on data.frame labels
        labels$ward = treeClust(p, num_clust ) %>% as.character()
        # NOTE: OVER WRITE the FINAL LABEL with the ward.D2 clustering from the distance matrix
        all_labels[[i]] = labels$ward

        # Make colors for all labels
        # Make palette of distinct colors
        cols_manual  = makeQualitativePal(length(labels$manual %>% unique()))
        names(cols_manual) <- labels$manual %>% unique() %>% sort()
        cols_spectrum = makeQualitativePal(length(labels$label %>% unique() ))
        names(cols_spectrum) <- labels$label %>% unique()
        # Since we are going to map this to the global clusters, let's choose colors in deterministic order
        tail_colors = ifelse(i==2, T, F)

        cols_ward = makeQualitativePal(length(labels$ward %>% unique() ), rand_order = F, tail_colors = tail_colors)
        names(cols_ward) <- labels$ward %>% unique() %>%  sort()

        col_annotation = list(manual = cols_manual, label = cols_spectrum, ward = cols_ward)

        if(quantile_norm){


            mat_breaks <- quantile_breaks(path_mat, n = 20)
            p2  = pheatmap(
                mat               = path_mat,
                color             = blues_pal(length(mat_breaks) - 1),
                cluster_cols      = F,
                breaks            = mat_breaks,
                border_color      = NA,
                clustering_distance_rows = dist_mat,
                show_rownames     = T,
                drop_levels       = TRUE,
                fontsize          =12,
                annotation_row = labels,
                annotation_colors = col_annotation,
                cutree_rows = num_clust
            )

        # Normal heatmap (no scaling + default palette)
        }else{
          pheatmap(path_mat, annotation_row = labels,
                   clustering_distance_rows = dist_mat,
                   cutree_rows = num_clust, show_rownames = F, fontsize =14, clustering_method = clustering_method_heatmap,
                    annotation_colors = col_annotation,border_color      = NA,drop_levels       = TRUE, cluster_cols      = F,
                    color = blues_pal(20))



      }
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




makeQualitativePal <- function(n, rand_order = T, skip = 0, tail_colors = F){

  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #pie(rep(1,n), col=sample(col_vector, n))
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
}


treeClust<-function(p,k){
  #p2 = pheatmap(input_matrix, clustering_distance_rows = dist.cosine(input_mat), silent = T)
  aa = as.hclust(p$tree_row)
  return(cutree(aa, k))
}


# Takes the data frame (motif_labels) from the previous function (clusterAllPathways)
# Makes a heatmap with the corresponding MI between all pairs of pathways
# Adjusted MI requires package aricode
# Returns the matrix of mutual information
pathwayConfusion<-function(motif_labels = data.frame(), clust_method = 'ward.D2'){

    # We are testing different metrics for clustering comparison:
    MI = matrix(0, dim(motif_labels)[2], dim(motif_labels)[2])
    N_MI = matrix(0, dim(motif_labels)[2], dim(motif_labels)[2])
    ari = matrix(0, dim(motif_labels)[2], dim(motif_labels)[2])
    Adj_MI = matrix(0, dim(motif_labels)[2], dim(motif_labels)[2])

    for( i in 1:(dim(motif_labels)[2]-1)){
        for(j in (i+1):dim(motif_labels)[2]){
            # Mutual information
            MI[i,j] = infotheo::mutinformation(motif_labels[,i], motif_labels[,j])
            # Adjusted mutual information
            Adj_MI[i,j] = AMI(motif_labels[,i], motif_labels[,j])
            # Normalized mutual information
            N_MI[i,j] = NMI(motif_labels[,i], motif_labels[,j])
            # Adjusted random index
            #ari[k] = adjustedRandIndex(motif_labels[,i], motif_labels[,j]) # Adjusted Rand Index
        }
    }
    resMI = Adj_MI
    resMI  = resMI + t(resMI)
    diag(resMI) <- 1
    row.names(resMI)<-colnames(motif_labels)
    colnames(resMI)<-colnames(motif_labels)
    pheatmap(resMI, fontsize = 14, col = magma(100), clustering_method  =clust_method)
    return(resMI)
}


# From the labels data frame we can now create a confusion matrix to visualize as a Sankey diagram:
# motif labels has pathways as columns, seurat clusters as rows and cluster-labels as values
# Returns a confusion matrix (not normalized by total sum so is not a P_x_y)
confMatrix_labels <- function(motif_labels, a,b){
	k = length(unique(motif_labels[,a]))
	kk = length(unique(motif_labels[,b]))
	conf_matrix = matrix(0,k,kk)

	a_classes = unique(motif_labels[,a])
	b_classes = unique(motif_labels[, b])
	for(i in 1:length(a_classes)){

	    # Find in b the data points that belong to the current A class
	    map_classes = motif_labels[which(motif_labels[,a]==a_classes[i]), b]

	    for(j in 1:length(b_classes)){
	        conf_matrix[i,j ] = sum(map_classes==b_classes[j])
	    }
	}

  row.names(conf_matrix) <- paste(a_classes, test_pathways[a],sep="_")
  colnames(conf_matrix) <- paste(b_classes, test_pathways[b],sep="_")

  return(conf_matrix)
}


# Sankey diagram for a pair of pathways
makeSankeyDiagram <- function(motif_labels, i, j){

      # Create confusion matrix
    conf_matrix = confMatrix_labels(motif_labels, i,j)

    data = as.data.frame(conf_matrix)
    # I need a long format
    data_long <- data %>%
      rownames_to_column %>%
      gather(key = 'key', value = 'value', -rowname) %>%
      filter(value > 0)
    colnames(data_long) <- c("source", "target", "value")
    data_long$target <- paste(data_long$target, " ", sep="")

    # From these flows we need to create a node data frame: it lists every entities involved in the flow
    nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())

    # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
    data_long$IDsource=match(data_long$source, nodes$name)-1
    data_long$IDtarget=match(data_long$target, nodes$name)-1

    # prepare colour scale
    ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

    # Make the Network
    sankeyNetwork(Links = data_long, Nodes = nodes,
                         Source = "IDsource", Target = "IDtarget",
                         Value = "value", NodeID = "name",
                         sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)

}

# Manual clustering
# Starts from the similarity matrix output by Spectrum and performs Spectral clustering using the normalized graph Laplacian
# this is the output from the spectral clustering function
# list[[1]] = data.frame with annotation and labels
# list[[2]] = Spectrum object
# list[[3]] = data matrix
manualSpectral <- function(A=c(), labels =data.frame(),k_neigh = 10){
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
	for(p in 2:5){
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



# Other functions

# Oct 16 2020 Works OK
# Seurat pipeline for analysis of early embryo atlas from https://www.nature.com/articles/s41586-020-2552-x/figures/5
# Analysis of scRNA mouse embryos at 5 different time-points
# Each time-point includes 9-11 embryos.
# Authors mentioned that they clustered data by time-point and batch correct by embryo so we are doing the same
integrate_Embryo <-function(s, matrix_list = list(),  k.filter = 100, npcs_umap = 30, npcs_clusters = 20, sample_id = c() , meta_df = data.frame(), clus.resolution = 1.5 ){
  # matrix_list: list of count matrices, created with Read10X function for each time-point (sample)
  # sample_id has all the time-points as character
  # meta_df is a master meta.data df. It partially was provided by the authors
  # the number of cells in meta_df and the count matrix do not exactly match, so we have to find shared barcodes
  # We filter embryos that do not pass a min number of cells (k.filter = 100)


  this_sample = sample_id[s]
  # Subset meta-data by time-point:
  sub_meta = meta_df %>% filter(sample==this_sample)
  #get the embryos for this time-point
  sub_meta$Embryo %>% unique() -> embryos_here
  sample_matrix_list = list()
  sample_meta_list = list()
  # Subset the list of count-matrices
  # Double check that barcodes are present in meta-data AND count matrix
  for(e in 1:length(embryos_here)){
      # subset meta_data by Embryo (we are in the same time-point here )
      embryo_cells = sub_meta %>% filter(Embryo==embryos_here[e])
      shared_cells = embryo_cells$cell_id[which(embryo_cells$cell_id %in% colnames(matrix_list[[s]]))]
      sample_matrix_list[[e]] = matrix_list[[s]][, shared_cells]
  		sample_meta_list[[e]] = embryo_cells %>% filter(cell_id %in% shared_cells)
  		row.names(sample_meta_list[[e]]) = sample_meta_list[[e]]$cell_id
  }

  # Now sample matrix contains all embryos as separated matrices for this time point
  # Second part: Seurat integration
  # in E6.5 some embryos have less that 100 cells (n =3)
  # After filtering we still get 7 embryos
  which_embryos_filter<-which(do.call(rbind,lapply(sample_matrix_list, dim))[,2]>k.filter)

  # Create seurat objects
  # Filter for those embryos that have enough cells
  tiss.embryo = list()
  for(i in 1:length(which_embryos_filter)){
      ef = which_embryos_filter[i]
      tiss.embryo[[i]] = CreateSeuratObject(counts= sample_matrix_list[[ef]], meta.data = sample_meta_list[[ef]])
  }

  # Apply SCT transform to all samples in the list
  for(i in 1:length(tiss.embryo)){
      tiss.embryo[[i]] <- SCTransform(tiss.embryo[[i]], verbose = F)
  }

  # min number of cells across datasets
  # Otherwise Seurat integration does not work


  # Change the memory parameters for parallel computing (future)
  # this dataset requires 5.6Gb so I set it to 6Gb
  options(future.globals.maxSize = 6000 * 1024^2)

  # select features for downstream analysis
  embryo.features <- SelectIntegrationFeatures(object.list = tiss.embryo, nfeatures = 3000)
  tiss.embryo <- PrepSCTIntegration(object.list = tiss.embryo, anchor.features = embryo.features,
                                    verbose = FALSE)



  embryo.anchors <- FindIntegrationAnchors(object.list = tiss.embryo, normalization.method = "SCT",
                                           anchor.features = embryo.features, verbose = FALSE, k.filter = k.filter)
  embryo.integrated <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT",
                                     verbose = FALSE)



  # From here, normal pipeline
  embryo.integrated <- RunPCA(embryo.integrated, verbose = FALSE)
  embryo.integrated <- RunUMAP(embryo.integrated, dims = 1:npcs_umap)


  embryo.integrated %>% FindNeighbors(dims = 1:npcs_umap) %>% FindClusters( resolution = clus.resolution) -> embryo.integrated

  # Set Default assay to RNA
  # And normalize
  DefaultAssay(embryo.integrated)<-'RNA'
  embryo.integrated <-NormalizeData(embryo.integrated)

  return(embryo.integrated)

}

# Tiss embryo is the output from integrate_Embryo
# create ann_heatmap manually with annotate function (from Motifs_figure1.R)
exportEmbryo <- function(tiss.embryo, ann_heatmap, pathways_file , all_pathway_file_names){
  # create the cell type field (used by the next function )
  ann_heatmap<- ann_heatmap %>% mutate(cell_type = cell_ontology_class)
  # import this function from test_Funxs.R
  export_df = exportSeuratPathways(tiss.embryo, pathway_file = pathways_file, all_pathway_file_names, cell_type_ann = ann_heatmap)
  return(export_df )
}



# Cardiomyocite dataset
# But also a more general pipeline for integration of replicates starting from a tiss.embryo list of SEurat objects

integrateTissEmbryo<-function(tiss.embryo = list() , k.filter = 200, npcs_umap = 30,  clus.resolution = 1.2, npcs_neigh = 20 ){

  # Apply SCT transform to all samples in the list
  for(i in 1:length(tiss.embryo)){
      tiss.embryo[[i]] <- SCTransform(tiss.embryo[[i]], verbose = F)
  }

  # min number of cells across datasets
  # Otherwise Seurat integration does not work


  # Change the memory parameters for parallel computing (future)
  # this dataset requires 5.6Gb so I set it to 6Gb
  options(future.globals.maxSize = 6000 * 1024^2)

  # select features for downstream analysis
  embryo.features <- SelectIntegrationFeatures(object.list = tiss.embryo, nfeatures = 3000)
  tiss.embryo <- PrepSCTIntegration(object.list = tiss.embryo, anchor.features = embryo.features,
                                    verbose = FALSE)



  embryo.anchors <- FindIntegrationAnchors(object.list = tiss.embryo, normalization.method = "SCT",
                                           anchor.features = embryo.features, verbose = FALSE, k.filter = k.filter)
  embryo.integrated <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT",
                                     verbose = FALSE)



  # From here, normal pipeline
  embryo.integrated <- RunPCA(embryo.integrated, verbose = FALSE)
  embryo.integrated <- RunUMAP(embryo.integrated, dims = 1:npcs_umap)


  embryo.integrated %>% FindNeighbors(dims = 1:npcs_neigh) %>% FindClusters( resolution = clus.resolution) -> embryo.integrated

  # Set Default assay to RNA
  # And normalize
  DefaultAssay(embryo.integrated)<-'RNA'
  embryo.integrated <-NormalizeData(embryo.integrated)

  return(embryo.integrated)
}




findK <- function(){
  set.seed(23)
  icMAt <- evaluateK(counts = lineage_counts, sds = sling_output2, k = 3:12, nGenes = 100, plot = T)
}


# Make binomial fit:
makeFit <- function(k = 5, BPPARAM){
  sce_5 <- fitGAM(counts = lineage_counts, sds = sling_output2, nknots = k,parallel=TRUE, BPPARAM = BPPARAM)
  return(sce_5)
}

# all tests from Trade Seq:
makeAlltests<-function(sce_8){

  assoRes <- associationTest(sce_8, lineage = T)

  startRes <- startVsEndTest(sce_8, lineages = T)

  endRes <- diffEndTest(sce_8, pairwise = T)

  patternRes <- patternTest(sce_8, pairwise = T)


  return( list (startRes, endRes, patternRes, assoRes ))
}



# Zhong 2018

makeAnnotatedHeatmap<-function(brain_lineage2, which.pathway =1 ){

  # Make data frame for annotation
 brain_lineage2@meta.data %>% group_by(cell_types, week, seurat_clusters) %>% count() %>%
     group_by(seurat_clusters) %>% mutate(freq= n/sum(n)) %>%
     arrange(desc(freq,seurat_clusters)) %>% top_n(1,freq) %>% arrange(seurat_clusters)-> cell_type_top_label_lineage

 #row.names(cell_type_top_label)<-0:(dim(cell_type_top_label)[1]-1)
 cell_type_top_label_lineage %>% as.data.frame() -> cell_type_top_label_lineage

 # arrange data frame by week
 cell_type_top_label_lineage %>% arrange(week) -> cell_type_top_label_lineage
 row.names(cell_type_top_label_lineage)<-cell_type_top_label_lineage$seurat_clusters


 # Make two color palettes to annotate by cell type and week
 ann_colors = list(cell_types =cell_type_colors, week=brewer.pal(cell_type_top_label_lineage$week %>% unique() %>% length(), 'Greys'))
 names(ann_colors$cell_types)<-brain_lineage2@meta.data$cell_types %>% unique()
 names(ann_colors$week)<-cell_type_top_label_lineage$week %>% unique()

 # Make heatmap
 bmp.mat = avg.matrix(brain_lineage2, pathways[[which.pathway]], by='seurat_clusters')
 p1 = pheatmap(bmp.mat,clustering_distance_rows = dist.cosine(bmp.mat), annotation_row = cell_type_top_label_lineage %>% select(cell_types, week), annotation_colors = ann_colors, cutree_rows = 6, cluster_cols = F)

 return(p1)
}



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
