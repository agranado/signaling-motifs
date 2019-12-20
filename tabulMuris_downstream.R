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



# calculate distance between cell types
# based on a list of >1000 TFs
# we will just compute euclidean dist and check if cell type annotations / dendrogram make sense.
# compare to tabula muris paper, they have a dendrogram of cell types in Fig 5


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
