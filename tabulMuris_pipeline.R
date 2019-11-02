#
#library(sva)
# this function assumes a global seurat object for tabula muris sce.seurat saved in ../tabula-muris
# this function assumes a scran_seurat object (saved in ../tabula-muris/)
batchCorrection_n_Imputation<-function( this_tissue = "Pancreas",k_magic = 4, npca = 25, t = "auto"){
  # 1. For each tissue
  #this_tissue = "Pancreas"

  # 2. Select the metadata for this tissue, including cell_type, plate (for batch correction)
  #this sce.seurat object contains SCRAN normalization + MAGIC
  sce.seurat@meta.data %>% filter(tissue ==this_tissue) %>% select(cell, plate.barcode,seurat_clusters,cell_ontology_class) -> pancreas_meta_data
  #rename so it is easy to use
  names(pancreas_meta_data)<-c("cell","plate","cluster","cell_type")
  row.names(pancreas_meta_data)<-pancreas_meta_data$cell
  # remove levels that are note relevant for this tissue
  pancreas_meta_data$plate = droplevels(pancreas_meta_data$plate)

  # scran_seurat is the SCRAN normalized Tabula Muris data BEFORE imputation with MAGIC!
  # this is saved as a separate object
  # Subset the cells only for this tissue from the SCRAN normalized counts,
  # subset genes with row.sum > 0, otherwise comBat crashes
  # batch.groups are the plates in the metadata
  pancreas_norm_batch = batch.normalise.comBat(
          as.data.frame(as.matrix(scran_seurat[rowSums(as.matrix(scran_seurat[, pancreas_meta_data$cell]))>0,pancreas_meta_data$cell])),
          batch.groups = pancreas_meta_data$plate)

  # we can now apply magic to the batch-corrected data
  # default parameters for magic
  pancreas_norm_batch_magic = Rmagic::magic(t(as.matrix(pancreas_norm_batch)),k=k_magic,n.jobs =-1, npca = npca, t = t )
  pancreas_norm_batch_magic = t(pancreas_norm_batch_magic$result)
  return(list(pancreas_norm_batch_magic,pancreas_meta_data))
}
###############################################
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

gaussianFitGeneExpression<-function(gene,k = 4, return.plot = T,expr_matrix = c()){
  # from library mixtools
  low_lim = 1
  frac_cells = 0.90




  if(length(expr_matrix)==0){
    expr = sce.seurat[['RNA']]@data[gene,]
  }else{
    expr = expr_matrix[gene,]
  }

  k = ifelse(sum(expr<low_lim)/length(expr)>frac_cells, 2,4)

  mixmdl <- normalmixEM(expr,k = k)
  df = data.frame(x = mixmdl$x)

  p = plotGaussianModel(df,mixmdl,k,gene)

  if(return.plot){
    return(p)
  }else{
    return(mixmdl)
  }

}

plotGaussianModel<-function(df,mixmdl,k,gene){

  if(k==4){
   p =  ggplot(df) +
    geom_histogram(aes(x, ..density..), binwidth = max(expr)/100, colour = "black",
                   fill = "white") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                  colour = "red", lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                  colour = "blue", lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[3], mixmdl$sigma[3], lam = mixmdl$lambda[3]),
                  colour = "green", lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[4], mixmdl$sigma[4], lam = mixmdl$lambda[4]),
                  colour = "purple", lwd = 1.5) +
    ylab("Density") +  ggtitle(gene)
  }else if(k==2){

    p =  ggplot(df) +
     geom_histogram(aes(x, ..density..), binwidth = max(expr)/100, colour = "black",
                    fill = "white") +
     stat_function(geom = "line", fun = plot_mix_comps,
                   args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                   colour = "red", lwd = 1.5) +
     stat_function(geom = "line", fun = plot_mix_comps,
                   args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                   colour = "blue", lwd = 1.5) +
     ylab("Density") +  ggtitle(gene)
  }


  return(p)
}

gaussianFitAssignCluster<-function(mixmdl,  rank_labels = c("0","L","M","H")){

  low_lim = 1
  frac_cells = 0.90


  #before assigning cluster labels let's check if the expression passes some filters
  #for VERY lowly expressed genes, we should reduce the number of clases
  if(sum(mixmdl$x<low_lim)/length(mixmdl$x) > frac_cells){
    #re do the cluster with k =2 , for ON/OFF labels
    expr= mixmdl$x
    mixmdl <- normalmixEM(expr,k = 2)
    rank_labels = c("0","L")
  }


  post.df <- as.data.frame(cbind(x = mixmdl$x, mixmdl$posterior))
  names(post.df)<-c("x",rank_labels[rank(mixmdl$mu)])


  assignLabel = function(x){names(which.max(x[2:5]))}
  post.df$group = apply(post.df,1,assignLabel)

  return(post.df$group)
}

profileDictionary<-function(profiles){

  # Take a list of profile
  # Compute the frequnecy of each one and make a table of unique profiles
  # Sort the table from most frequent to least frequent
  profile_ranks = rep(0,length(profiles))
  table(profiles) %>% sort(decreasing = T)  ->sorted_profiles

  # For each profile, we will find it's rank in the table and assign it to an array
  # Such that we have an array with the same dimension and the rank as its value (i.e, there are repeated values)
  for(i in 1:length(sorted_profiles)){

    profile_ranks[which(profiles == names(sorted_profiles)[i])] = i
  }

  return(profile_ranks)
}

plotMotifBarPlot<-function(bmp_quant){

  motif_values = data.frame(genes = bmp.receptors, value = as.numeric(str_split(bmp_quant,pattern = "")[[1]]) )

  p = ggplot(motif_values, aes(x = genes,y=value)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}

fitSingleGene<-function(gene,plot.all = F,expr_matrix = c(),G= 1:9){

  if(length(expr_matrix)==0){
    x = sce.seurat[['RNA']]@data[gene,];
  }else{
    x = expr_matrix[gene,]
  }

    fit =Mclust(x,G = G)

    if(plot.all){
      x11();
      par(mfrow=c(3,1));
      hist(x,breaks = 100)
      plot(fit,what="density");plot(fit,what="classification");

    }else{
      return(fit)
    }
}

## statistics


#plotMotifDist
#this_tissue is the name of the variable
# class is the type of variable:
#cell_ontology_class for cell types
#tissue for tissues
#seurat_clusters for clusters (based on transcriptome)
plotMotifDist<-function(this_tissue,tabula = tabula, class = "tissue",binwidth = 3){
  #this_tissue = "Brain_Non-Myeloid";



  if(class=="tissue"){
    tabula %>% group_by(tissue,pathway) %>% summarise(count = n()) %>% mutate(freq = count/sum(count)) %>% arrange(tissue,desc(count)) -> motifs_by_tissue
    motifs_by_tissue %>% filter(tissue ==this_tissue) %>% arrange(desc(freq), pathway) ->df
  }else if(class=="cell_ontology_class"){
    tabula %>% group_by(cell_ontology_class,pathway) %>% summarise(count = n()) %>% mutate(freq = count/sum(count)) %>% arrange(cell_ontology_class,desc(count)) -> motifs_by_tissue
    motifs_by_tissue %>% filter(cell_ontology_class ==this_tissue) %>% arrange(desc(freq), pathway) ->df
  }else if(class =="seurat_clusters"){
    tabula %>% group_by(seurat_clusters,pathway) %>% summarise(count = n()) %>% mutate(freq = count/sum(count)) %>% arrange(seurat_clusters,desc(count)) -> motifs_by_tissue
    motifs_by_tissue %>% filter(seurat_clusters ==this_tissue) %>% arrange(desc(freq), pathway) ->df
  }

  p =ggplot(df, aes(x = freq*100, y = ..density..)) +
      geom_histogram(,position = "identity",binwidth = binwidth) + geom_density() +
      ggtitle(this_tissue) + theme(text = element_text(size=20))  + xlab(" %cells with motif  ") +
      ylab("Density") +  xlim(-5, 100)

  return(p)
}

 # # # # # # # # # # 
# # # # # # # # # #
# global meta data object:
if(!exists("tabula")){
  sce.seurat@meta.data %>% select(cell, tissue, cell_ontology_class, seurat_clusters) -> tabula
}

# GENERAL PIPELINE, Nov 1st 2019
# This function takes a list of genes,
# Retrieves data (normalized and imputed) from the Seurat object (or a matrix, if provided)
# Applies a GMM to each gene across all cells with either k = 2 for low expressed and k =4 for all other genes
# I takes the labels from the GMM and created a word with size N= number of genes
# For the general data-frame, each single cell now has a word associated with it that represents the expression profile
# The rank of each profile: inverse frequency is also calculated
# The output data frame can then be used as input for plotting and statistics functions
pipelineNov2019<-function(gene.list = bmp.receptors,tabula){

      # Modelng: Gaussian mixture model for each gene across all cell types
      model_list = lapply(gene.list,gaussianFitGeneExpression,return.plot=F)
      # Discretization: Label each gene, for each cell
      all_labels = lapply(model_list, gaussianFitAssignCluster)
      label_matrix = do.call(cbind,all_labels)

      # Based on the discrete vales, we create a list of words.
      profiles = apply(do.call(cbind,all_labels),1,paste,sep="",collapse="")

      # Each cell has now a profile assigned
      tabula$pathway = profiles

      profile_ranks <- profileDictionary(profiles)

      tabula$pathway_rank = profile_ranks

      tabula$pathway_rank_string = as.character(tabula$pathway_rank)
      tabula$pathway_rank_string<-factor(tabula$pathway_rank_string,levels = as.character( sort(unique(tabula$pathway_rank))  ))

      # Also we use a string with numerical values (probably more useful)
      tabula$pathway_quant = unlist(lapply(tabula$pathway, str_replace_all, c( "L" = "1", "M" = "2","H" = "3")))
      #we want to return the data frame with the list of motifs, and ranks, that can be used for plotting
      return(tabula)
}

# INTERNAL 
groupMotifsBy<-function(tabula, class ="seurat_clusters"){
  if(class =="seurat_clusters"){
    tabula %>% group_by(seurat_clusters,pathway_quant) %>% summarise(count = n()) %>% mutate(freq = count/sum(count)) %>% arrange(seurat_clusters,desc(count)) -> motifs_by_cluster

  }else if(class =="tissue"){

    tabula %>% group_by(tissue,pathway_quant) %>% summarise(count = n()) %>% mutate(freq = count/sum(count)) %>% arrange(tissue,desc(count)) -> motifs_by_cluster
  }else if(class =="cell_ontology_class"){

    tabula %>% group_by(cell_ontology_class,pathway_quant) %>% summarise(count = n()) %>% mutate(freq = count/sum(count)) %>% arrange(cell_ontology_class,desc(count)) -> motifs_by_cluster
  }

  return(motifs_by_cluster)
}

# USE THIS
kMeansOnMotifs<-function(tabula, k = 100, class = "seurat_clusters"){

    motifs_by_cluster = groupMotifsBy(tabula, class)
      #extract profiles as string of numbers
    all_motifs_quant_bycluster<-motifs_by_cluster$pathway_quant

    # split string, convert to number and arrange as matrix
    all_motifs_quant_number= do.call(rbind,lapply(lapply(all_motifs_quant_bycluster,str_split,"",simplify=T),as.numeric))

    #Cluster by k-means with k = 100
    k = kmeans(all_motifs_quant_number,100)

    #save the labels as a new field in the motif data-frame
    motifs_by_cluster$motif_class = k$cluster

    motifs_by_cluster
}



## # # # # plottting moitfs
# HEATMAP after kmeans
heatmapKMeansClasses<- function(motifs_by_cluster, class = "seurat_clusters"){

  clust_vs_motifs = data.frame()

  if(class == "seurat_clusters"){
    motifs_by_cluster %>% group_by(seurat_clusters,motif_class) %>% summarise(class_freq = sum(freq)) %>% arrange(seurat_clusters, desc(class_freq)) %>%
        filter(class_freq>0.1) %>% spread(key = seurat_clusters,value= class_freq,fill =0) -> clust_vs_motifs;
  }else if( class =="tissue"){
    motifs_by_cluster %>% group_by(tissue,motif_class) %>% summarise(class_freq = sum(freq)) %>% arrange(tissue, desc(class_freq)) %>%
        filter(class_freq>0.1) %>% spread(key = tissue,value= class_freq,fill =0) -> clust_vs_motifs;
  }else if(class =="cell_ontology_class"){
    motifs_by_cluster %>% group_by(cell_ontology_class,motif_class) %>% summarise(class_freq = sum(freq)) %>% arrange(cell_ontology_class, desc(class_freq)) %>%
        filter(class_freq>0.1) %>% spread(key = cell_ontology_class,value= class_freq,fill =0) -> clust_vs_motifs;
  }

  # common for all classes
  clust_vs_motifs %>% select(-motif_class) %>% as.data.frame() -> clust_vs_motifs_mat ;

  row.names(clust_vs_motifs_mat)  = clust_vs_motifs$motif_class ;
  p= pheatmap(clust_vs_motifs_mat,fontsize = 12);
  x11(); plot(p[[4]]) ;

}

makePlots(tabula = tabula){




}

atLeastN_motifsWithPercent<-function(tabula = tabula, howMany = 1){
  n_clusters = c(); vals = seq(0,1,0.05)
  for(i in 1:length(vals)){
    motifs_by_cluster %>% mutate(is_more_than = freq>vals[i]) %>% group_by(seurat_clusters) %>%
        summarise(motifs =sum(is_more_than)) %>% mutate(motifs_here = motifs >= howMany) %>% summarise(clusters_with_motifs =sum(motifs_here)) -> n_clusters_with_motif
    n_clusters[i] =  n_clusters_with_motif$clusters_with_motifs
  }
  return(n_clusters)
}
