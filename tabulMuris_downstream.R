# analysis of motifs and downstream genes
# After meeting Jordi, Dec 17th 2019

# Alejandro Granados
source("tabulMuris_pipeline.R")



# # # # DOWNSTREAM GENES
# # # #

# We want to compare motifs from one pathway to other pathway or to downstream gene targets.
# We need a data frame where every single cell has a pathway_1 ID , pathway_2 ID, ..
# Then, we can compute MI, or some other correlation metric
pipelineSingleCells<-function(tabula = tabula, gene.list_1 = bmp.receptors, gene.list_2 = bmp.downstream ,
                                      k1 = 100, k2 = 100, expr_matrix = c(), control = F){

  # run the pipeline with the first list of genes
  aa<-pipelineNov2019(tabula = tabula, gene.list = gene.list_1, expr_matrix = expr_matrix)

  bb<-pipelineNov2019(tabula = tabula, gene.list = gene.list_2, expr_matrix = expr_matrix)

  kMeansOnMotifs(aa, class = 'cell_ontology_class') -> aa
  kMeansOnMotifs(bb, class = 'cell_ontology_class') -> bb

  # rename the column so we can merge both data frames
  bb %>% rename(target_class = motif_class) -> bb
  # from the second df, just take the relevant columns
  merge(aa, select(bb,cell,target_class)) -> df


  return(df)
}
  # count the co-ocurrence of classes between target and pathway motifs

# this function plots signaling profiles vs downstream profiles, regarless of the cell type
# it can be useful to show that the distribution is not random, as the matrix shows blocks of signaling-downstream specific activity
motifVsDownstream<-function(df,fil_min = 5){
  # from
  # https://stackoverflow.com/questions/17538830/x-axis-and-y-axis-labels-in-pheatmap-in-r/21684103
  # and from figure3.R, Tabula Muris paper

  df %>% group_by(motif_class,target_class) %>% summarise(count = n()) %>% filter(count > fil_min) %>% spread(key = motif_class, value = count, fill = 0)-> hmap_df

  hmap_mat <- as.data.frame(hmap_df %>% ungroup() %>% select(-target_class))
  row.names(hmap_mat) <- hmap_df$target_class


  # make the heatmap
  # use the grid library to label the x,y axis
  setHook('grid.newpage',function() pushViewport(viewport(x =1,y=1,width=0.9,height =0.9,name ='vp',just =c('right','top'))),action ='prepend')


  pheatmap(log(hmap_mat + 1))
  setHook("grid.newpage", NULL, "replace")

  grid.text("Motif class signaling", y=-0.07, gp=gpar(fontsize=16))
  grid.text("Motif class downstream", x=-0.07, rot=90, gp=gpar(fontsize=16))

}


# We want to analyse the profiles in separate ways:
# Some are specific for a given cell type
# Some appear in multiple cell types
# Questions:
#  1 Do cells that show the specific profiles also show specific downstream gene expression?
#  2 Do disparate cell types with the same signaling profile also show the same signaling downstream expression?
# We want to heatmaps, of profiles vs downstream.
analyseByCellType <-function(df, min_frac = 0.10){
  # We first need to classify profiles as motifs or specific, whether they appear in multiple or single cell types
  # Let's count the number of cell types a profile appears in :
  # We define a profile as being present if it appear in > min_frac of the cells for that cell type

  df %>% drop_na(cell_ontology_class) %>% mutate(anno_and_tissue = paste0(cell_ontology_class, " (", tissue, ")")) %>%
  drop_na(anno_and_tissue) -> df_anno

  df_anno %>%  group_by(anno_and_tissue,motif_class) %>% summarise(count = n()) %>% mutate(freq = count/sum(count)) %>% arrange(anno_and_tissue,desc(count)) -> motifs_df

  # this works
  # quantify in how many different cell types the profiles appear, consider only > min_frac
  # then label them as motif or specific
  motifs_df %>% filter(freq>min_frac) %>% group_by(motif_class) %>% summarise(n_tissues = n()) %>% mutate(profile_type = ifelse(n_tissues==1, 'specific' , 'motif')) -> profile_types

  # now we know the type of profile, we can go back to the original data frame and fill a new column with the profile_type
  # finally: original df with the annotation of profile type

  motifs_df %>% filter(freq >min_frac) %>% merge(profile_types) -> reference_motifs


  return(list(df_anno,reference_motifs))
}
# For each pair of motif_class and anno_tissue:
#   Get the expression data, compute average, save in matrix
# Plot 2 heatmaps:
#  1. Specific
#  2. Motifs
# Colums are downstream genes, rows are motif_class, cell_type pairs
retrieveDownstreamProfiles<-function(df_anno,reference_motif,targets = bmp.downstream,class = 'specific', expr_matrix = c()){

  df_anno %>% filter(motif_class==1, anno_and_tissue=='lung endothelial cell (Lung)') %>% select(cell) -> this_cluster

  reference_motifs %>% filter(profile_type == class) -> spec_profiles
  all_mat = matrix(0,dim(spec_profiles)[1],length(targets))

  for( i in 1:dim(spec_profiles)[1]){
    df_anno %>% filter(motif_class==spec_profiles[i,]$motif_class, anno_and_tissue==spec_profiles[i,]$anno_and_tissue) %>% select(cell) -> this_cluster
    expresion_this_cluster = expr_matrix[targets,this_cluster$cell]
    mean_expr = apply(expresion_this_cluster,1,mean)
    all_mat[i,] = mean_expr
  }

  row.names(all_mat)= paste(spec_profiles$motif_class," ",spec_profiles$anno_and_tissue, sep ="")
  colnames(all_mat) = targets
  return(all_mat)
}
# # # # END DOWNSTREAM
# # # # #


# # # # JAN 2020
# calculate distance between cell types
# based on a list of >1000 TFs
# we will just compute euclidean dist and check if cell type annotations / dendrogram make sense.
# compare to tabula muris paper, they have a dendrogram of cell types in Fig 5

# 1. calculate distance between cell types
# 2. Run the pipeline and get the motifs
# 4. Calculate the distance between motifs and create a matrix
# 5. Map the two matrices using the cell type ID
# 6. Make a 2D plot dist_pathway vs dist_cell_types

distMatCellTypes <- function(tabula){
  # get the TF names
  tfs      <- read.csv('../tabula-muris/23_tf_analysis/GO_term_summary_20171110_222852.csv')
  tf.names <- as.character( tfs %>% distinct(Symbol) %>% pull(Symbol) )
  tf.names <- make.names(tf.names)
  tf.names <- tf.names[tf.names %in% rownames(magic_norm)]
  # we want to group them by cell type ontology:
  tabula %>% drop_na(cell_ontology_class) -> tabula
  all_types = unique(tabula$cell_ontology_class)

  all_mat = matrix(0,length(all_types),length(tf.names))
  # which matrix are we using to get the expression data
  for(i in 1:length(all_types)){
    these_cells <- tabula %>% filter(cell_ontology_class ==all_types[i])
    all_mat[i, ] = apply(norm_counts[tf.names,these_cells$cell],1,mean)

  }
  row.names(all_mat) = all_types
  colnames(all_mat) = tf.names



  return(all_mat)
}

runGeneList <-function(gene_list = wnt.genes){

    result = runFullPipeline(tabula,input_matrix = sce.seurat[['RNA']]@data[gene_list,tabula$cell], gene.list = gene_list,control_type = 0,group_classes = 'cell_ontology_class')

    result %>% filterNonExpressing() %>% kMeansOnMotifs(class = 'cell_ontology_class') -> result_full
    #result_full %>% filter(n_expres>5) %>% makeAllPlots(this_pathway = "",k=50,class = 'cell_ontology_class')

    return(result_full)
}
# result_full after running the pipeline
# min_freq, minimum fraction of cells in a cell type expressin the motif to call it 'expressed'
# min_expres, min fraction of genes in the pathway that have values > 0
library(philentropy)
distMatPathways <- function(result_full, gene_list = wnt.receptors, min_freq = 0.1, min_expres=0.2,dist_method='corr'){


    result_full %>% filter(freq>=min_freq & n_expres/nchar(result_full$pathway_quant[1])>min_expres) -> high_freq_motifs


    high_freq_matrix = do.call(rbind,lapply(lapply(high_freq_motifs$pathway_quant,str_split,'',simplify=T) ,as.numeric )  )


    # let's name motifs with the cell type they appear in
    # this way we can then calculate distance and compare with distance between cell types
    row.names(high_freq_matrix)<-as.character(high_freq_motifs$cell_ontology_class)
    colnames(high_freq_matrix)<-gene_list

    if(dist_method =='corr'){
    # this has names of cell types!
    # and they could be duplicated
    dist_pathway = 1-cor(t(high_freq_matrix))
  }else if(dist_method=='euclidean'){

    dist_pathway = distance(high_freq_matrix)
    row.names(dist_pathway) = row.names(high_freq_matrix)
    colnames(dist_pathway ) = row.names(high_freq_matrix)
  }

  return(dist_pathway)
}

compareDistances<-function(dist_pathway, dist_cell_type){

  n_motifs = dim(dist_pathway)[1]
  pathway_cell_types = row.names(dist_pathway)

  all_pairs_dist_pathway = c()
  all_pairs_dist_celltype = c()
  c = 1


  for(i in 1:(n_motifs-1 )){
    # just the upper triangle matrix
    for(j in (i+1):n_motifs ){
      all_pairs_dist_pathway[c] = dist_pathway[i,j]
      all_pairs_dist_celltype[c] = dist_cell_type[ pathway_cell_types[i]  , pathway_cell_types[j]  ]
      c = c+1
    }
  }

  return(list(all_pairs_dist_pathway,all_pairs_dist_celltype))
}


# # # # #
# # # # #
# Jan 8th 2020, RANDOM ensemble
# GOAL: to create a null hypothesis for number of expected clusters
# ideally the input is already sampled, so the matrix is the dimension of length(pathway) * repeats
# this will speed up the computation and save memory crash
# the list of genes is just to get the size of the pathway


runRandomEmsemble <-function(list_genes = bmp.receptors,n_repeats = 100,input_matrix =c()){

  result_list= list()
  # the matrix is composed of a random list of genes, by the way it was constructed
  # So here we just take the row.names which is already the random set
  all_random_genes = row.names(input_matrix)
  i = 1
  # this is the fundamental N, which is the size of the original pathway (and therefore size of fake pathways)
  n_path = length(list_genes)
  for(i in 1:n_repeats){
    print(paste('Running repeat', i, '...'))
    gene_sublist = all_random_genes[  (n_path*(i-1)+1):(n_path*i)];
    # we take the random pathways in sets of N, one by one and apply the pipeline
    # subset the matrix with the current batch
    input_sub_matrix = input_matrix[gene_sublist,];
    # runFullPipeline, default parameters (No kmeans )
    result = runFullPipeline(tabula= tabula, input_matrix = input_sub_matrix,gene.list =gene_sublist,group_classes = "seurat_clusters")
    result_list[[i]] = result
  }

  # return a list of dataframes with profiles frequencies by cell type
  return(result_list)
}

# next function, take the list from the previous function and perform clustering
# We need the optimal number of clusters as the output of the following function
BIC_optimization<-function(result_list){


  quant_list = lapply(result_list,pipelineResToQuant)
  #quant_list is a matrix with all the profiles, discretized (GMM k=4)
  # they have different sizes, depending on the pathway (control, etc)
  d_clust_list <- lapply(quant_list,Mclust, G=1:100)

  return(d_clust_list)
}
