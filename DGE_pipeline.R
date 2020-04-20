# DGE analysis for group meeting

library(ggplot2)
library(pheatmap)
library(igraph)
library(dplyr)
library(tidyr)
library(viridis)

#load("globalClusters_analysis_pathwayDistanceMarch2020.Rdata")
#rm(res.list)
#load("../groupMeeting2020/dataset/Differential_exp_5pathways.rda")
#load("/Users/alejandrog/Downloads/Avg_matrix_mar28_manualPathways.rda")


seve.receptors=c("Acvr1", "Acvrl1", "Bmpr1a", "Bmpr1b", "Acvr2a", "Acvr2b",  "Bmpr2")
notch.genes = c("Dll1",   "Dll3"   ,"Dll4",   "Dtx1",   "Jag1"  , "Jag2", "Notch1", "Notch2", "Notch3", "Notch4", "Mfng",  "Rfng"  , "Lfng",   "Dlk1",   "Dlk2")


list_of_pathways=list(bmp.receptors,notch.genes, wnt.receptors,lpa_genes,nf_kb$Gene.Name)
pathway_names = c(  'BMP_R',
                    'Notch_R_lig',
                    'Wnt_R',
                    'LPA_R_L',
                    'NFkb_R_L')


dendrogramDGE<-function(i,this.res.list = res.list,
                          pathway_list = list_of_pathways){

    dist_pathway = distMatrix(this.res.list,all_clusters, which.genes=pathway_list[[i]] , p_val = 10^-6,subset_genes = T)
    h_tree = hclust(as.dist(dist_pathway))
    return(as.dendrogram(h_tree))
}


# Given a real pathway, make k lists with random genes from the control matrix
# We use the 10k variable gene matrix as default
makeRandomPathway<-function(p,k =10,
                            source_matrix = avg_Matrix_March28,
                            control_matrix = variable10k_matrix,
                            range.mean = 0.1,
                            range.var = 0.1){

 # In this way we can sample genes that have similar statistics as the real pathway
	 mean_rand300= apply(control_matrix,2,mean)
 	var_rand300= apply(control_matrix,2,sd)
	list_of_rand_pathways = list()
	rand_pathway_names = c()

	#### MAKE RANDOM PATHWAYS
	# make k random pathways from the control_matrix list
	for(i in 1:k){
    	list_of_rand_pathways[[i]] = makeRandomMean(mean_rand300,var_rand300,3,source_matrix[,list_of_pathways[[p]]],range.mean = range.mean,range.var =range.var)
		rand_pathway_names[i] = paste("rand","_",toString(i),sep="")
	}

	return(list_of_rand_pathways)
}




#
#
# doControlPlots<-function(p){
#
#   list_of_rand_pathways = list()
#   rand_pathway_names = c()
#   for(i in 1:10){
#       list_of_rand_pathways[[i]] = makeRandomMean(mean_rand300,var_rand300,3,avg_all_82Clusters[,list_of_pathways[[p]]],range.mean = 0.2,range.var = 0.3)
#       rand_pathway_names[i] = paste("rand","_",toString(i),sep="")  }
#
#   x11(); par(mfrow=c(2,5))
#   for(i in 1:length(rand_pathway_names)){
#       graphs = analysisNoPlot(i, list_of_pathways = list_of_rand_pathways , pathway_names = rand_pathway_names , rand_matrix_3000,res.listAll)
#       distributionCellTypeDist(graphs[[3]],title = rand_pathway_names[i])     }
#
# }

#Mar31
# returns the mean celltype distance for a given set of target global clusters in the celltype space
# The methods define how this distance is computed
cellTypeDist<-function(target_genes,method = "DGE",cellType_matrixDGE = dist_10k){

  if(method =="DGE"){
    cell_type_dist = cellType_matrixDGE[target_genes,target_genes]
    mean_dist_celltypes = mean(cell_type_dist[upper.tri(cell_type_dist)])
  }else if(method =="hclust"){
    cell_type_dist = dist_from_hclust_10k[target_genes,target_genes]
    mean_dist_celltypes = mean(cell_type_dist[upper.tri(cell_type_dist)])
  }

  return(mean_dist_celltypes)

}

# returns a data frame (Mar 31)
# for a given pathway
# size is the size of the cliques and cellType_dist is their distance in the celltype space according to method
# df$size, df$cellType_dist
distributionCellTypeDist<-function(this_graph, title = "",min.clique = 2,method = "DGE", whichCellMatrix = dist_10k){
  cliques_all = cliques(this_graph, min =2)
  clique_dist = unlist(lapply(cliques_all,function(x){ cellTypeDist(x$name,method, cellType_matrixDGE = whichCellMatrix) } ))
  #hist(clique_dist,breaks = 20,
      #main= title,xlab = "Avg Cell-Type Distance within block")

  df_output=data.frame(size =unlist(lapply(cliques_all,length)),celltype_dist = clique_dist,names = unlist(lapply(cliques_all, function(x){ paste(x$name,collapse = "_") })))
  return(df_output)
}

# Performs DGE analysis and
# returns a list of graphs (graphs[[3]] is the one where we filter for lowly expressed global clusters)
analysisNoPlot<-function(i, list_of_pathways_ = list(), pathway_names_ = c(),expr.matrix,res.list,min.pct.genes = 0.2){

    genes_pathway = list_of_pathways_[[i]]             #select list of genes to subset matrix

    expr.mat = expr.matrix[,genes_pathway]     # Express

    which.genes = genes_pathway;
    dist_pathway = distMatrix(res.list,all_clusters, which.genes=which.genes, p_val = 10^-6,subset_genes = T)


    adj_pathway = as.matrix(as.data.frame(dist_pathway)==0)
    graph_from_adjacency_matrix(adj_pathway,mode='undirected',diag=F) -> graph1

    singletons = apply(dist_pathway, 2, function(x){sum(x==0)})==1
    adj_pathway_connected = adj_pathway[!singletons,!singletons]

    graph_from_adjacency_matrix(adj_pathway_connected,mode='undirected',diag=F)  -> graph2

    #min.pct =min.pct.genes  # only profiles having 0.2 of their genes expressed at higher level than the mean across clusters
    min.genes.expr =  floor(length(list_of_pathways_[[i]]) * min.pct.genes)
    pass.filter <- names(which( rowSums(expr.mat>mean(expr.mat)) >=min.genes.expr))

    adj_pathway_connected = adj_pathway[pass.filter, pass.filter]
    graph_from_adjacency_matrix(adj_pathway_connected,mode='undirected',diag=F) -> graph3
    return(list(graph1,graph2,graph3))
}


runDGEonControl<-function(i, list_rand =list(), path_names = c(), control_mat = c(), control.res = list(), dist_method = "DGE" ){
  graphs = analysisNoPlot(i,
                          list_of_pathways_ = list_of_rand_pathways ,
                          pathway_names_ = rand_pathway_names , control_matrix,control_res.list)
  dist_cliques =  distributionCellTypeDist(graphs[[3]],title = "", method = dist_method)
  return(dist_cliques)
}
# path_cliques has the REAL pathways runned already
# Be sure to run with the real pathways using the same parameters, to make good comparisons

runControlDistance<-function( p=1 , source_matrix = avg_Matrix_March28,
                              list_of_pathways = list(bmp.receptors,notch.genes, wnt.receptors,lpa_genes,nf_kb$Gene.Name,seven.receptors),
                              path_cliques_real = path_cliques,
                              pathway_names = c(
                                              'BMP_TGFb_R',
                                              'Notch_R_lig',
                                              'Wnt_R',
                                              'LPA_R_L',
                                              'NFkb_R_L',
                                              'BMPr'),
                              control_matrix = variable10k_matrix, control_res.list = res.list_10k, celltype_matrix = dist_10k,
                              dist_method = "DGE",
                              k = 10,range.mean = 0.2, range.var = 0.2,
                              save.pdf = F, save.path = "../groupMeeting2020/plots/CellTypeDist/",id_file = "10k",
                              pdf.width = 20, pdf.height = 6){

  # The mean and variance of the control matrix
  # In this way we can sample genes that have similar statistics as the real pathway
  mean_rand300= apply(control_matrix,2,mean)
  var_rand300= apply(control_matrix,2,sd)
  # Filter very low expressed genes for the random sets
  keep_genes = mean_rand300>0.05
  mean_rand300 = mean_rand300[keep_genes]
  var_rand300 = var_rand300[keep_genes]

  list_of_rand_pathways = list()
  rand_pathway_names = c()

  #### MAKE RANDOM PATHWAYS
  # make k random pathways from the control_matrix list
  for(i in 1:k){
      list_of_rand_pathways[[i]] = makeRandomMean(mean_rand300,var_rand300,3,source_matrix[,list_of_pathways[[p]]],range.mean = range.mean,range.var =range.var)
      rand_pathway_names[i] = paste("rand","_",toString(i),sep="")  }

  ### RUN THE CLIQUES

  path_cliques_rand = list()

  # RETURN A DF
  for(i in 1:length(rand_pathway_names)){
      graphs = analysisNoPlot(i,
                              list_of_pathways_ = list_of_rand_pathways ,
                              pathway_names_ = rand_pathway_names , control_matrix,control_res.list)
      path_cliques_rand[[i]] = distributionCellTypeDist(graphs[[3]], method = dist_method, whichCellMatrix = celltype_matrix)
  } # ENDFOR

  # in parallel (Mar 31st)
  # path_cliques_rand = mclapply(1:8, runDGEonControl, list_rand = list_of_rand_pathways,
  #                             path_names =   rand_pathway_names, control_mat = control_matrix,
  #                             control.res = control_res.list   , mc.cores = 8)

      # MERGE ALL CONTROLS INTO ONE
      wnt_rand_control10<-do.call(rbind,path_cliques_rand)
      wnt_rand_control10$type = rep('rand',dim(wnt_rand_control10)[1])

      # MERGE WITH REAL PATHWAY
      wnt_cliques<-path_cliques_real[[p]]
      wnt_cliques$type = rep(pathway_names[p],dim(wnt_cliques)[1])

      # DONE
      wnt_all = rbind(wnt_cliques,wnt_rand_control10)
      wnt_all %>% ggplot(aes(x = as.factor(size), y = celltype_dist, fill = type)) + geom_boxplot() +
          theme_minimal() +
          theme(text = element_text(size=20)) + ggtitle(pathway_names[p]) +
          ylab("Cell type diversity") + xlab("Size of motif (N of profiles)") +
          theme(legend.title = element_blank()) + scale_fill_brewer(palette="Spectral")

      if(save.pdf){
        ggsave( paste(save.path,"cellType_dist_",dist_method,"_",pathway_names[p],"_",id_file,".pdf") ,width = pdf.width ,height = pdf.height)
      }

      return(wnt_all)

}



# This function takes a reference array with mean and standard deviation
# It then samples from this array accordingly to the target pathway in expr.mat
# Works fine BUT the 10k set of genes does not contain highly expressed genes...
makeRandonMean <-function(mean_rand300,var_rand300,i,expr.mat, range.mean =0.1,range.var =0.2){

    mean_target_pathway<-apply(expr.mat,2,mean)
    var_target_pathway<-apply(expr.mat,2,sd)
    rand_pathway  = c()
    for(g in 1:dim(expr.mat)[2]){
        subset_same_stats = mean_rand300> mean_target_pathway[g]-range.mean & mean_rand300< mean_target_pathway[g] + range.mean
        subset_same_stats = subset_same_stats & var_rand300> var_target_pathway[g]-range.var & var_rand300< var_target_pathway[g] + range.var
        rand_pathway[g] = sample(names(which(subset_same_stats)),1)
    }
    return(rand_pathway)
}





doMainPlots<-function(i,plots.path="../groupMeeting2020/plots/",
                      expr.matrix = avg_all_82Clusters,
                      res.list = list(),
                      list_of_pathways = list(),
                      pathway_names,
                      k = 8){

    #Starts here
    # Gene expression Matrix

    pdf(paste(plots.path,pathway_names[i],"_DGEplots.pdf",sep =""),width = 12,height = 15)


                                              # select pathway
    genes_pathway = list_of_pathways[[i]]             #select list of genes to subset matrix

    expr.mat = expr.matrix[,genes_pathway]     # Expression Matrix
    p = pheatmap(expr.mat,fontsize = 15)              # Plot 1


    quantileHeatmap(expr.mat,pathway_names[i])

    # Distance matrix
    which.genes = genes_pathway;
    dist_pathway = distMatrix(res.list,all_clusters, which.genes=which.genes, p_val = 10^-6,subset_genes = T)
    #zeroCount_pathway = countProfilesZeroBlocks(dist_pathway,clust_method = 'complete')
    #zeroCount_pathway %>% unique() %>% length()

    p2 = pheatmap(dist_pathway,fontsize = 15)

    # # Network
    # adj_pathway = as.matrix(as.data.frame(dist_pathway)==0)
    # graph_from_adjacency_matrix(adj_pathway,mode='undirected',diag=F) %>%     plot(vertex.colod ='gold', vertex.size=9)
    # #With no singletons in the graph
    #
    # singletons = apply(dist_pathway, 2, function(x){sum(x==0)})==1
    # adj_pathway_connected = adj_pathway[!singletons,!singletons]
    #
    # graph_from_adjacency_matrix(adj_pathway_connected,mode='undirected',diag=F) %>%
    #   plot(vertex.colod ='gold', vertex.size=9,main = paste("DGE network ",pathway_names[i], ". Only connected profiles",sep="") )

    # plot a filtered version of the distance matrix
    min.pct = 0.2 # only profiles having 0.2 of their genes expressed at higher level than the mean across clusters
    min.genes.expr =  floor(length(list_of_pathways[[i]]) * min.pct)
    pass.filter <- names(which( rowSums(expr.mat>mean(expr.mat)) >=min.genes.expr))
    p = pheatmap(dist_pathway[pass.filter,pass.filter],fontsize = 14,silent = F,main = pathway_names[i],cutree_rows = k)


    # adj_pathway_connected = adj_pathway[pass.filter, pass.filter]
    # graph_from_adjacency_matrix(adj_pathway_connected,mode='undirected',diag=F) %>%
    #     plot(vertex.colod ='gold', vertex.size=9,main = paste("DGE network ",pathway_names[i], ". Only connected profiles",sep="") )


    dev.off()



}




quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

quantileHeatmap<-function(expr.mat,pathway.name){

  mat_breaks <- quantile_breaks(expr.mat, n = 20)

  p2  = pheatmap(
    mat               = expr.mat,
    color             = inferno(length(mat_breaks) - 1),
    breaks            = mat_breaks,
    border_color      = NA,
    show_colnames     = T,
    show_rownames     = T,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = pathway.name
  )


}


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



nGenesDifferent<-function(res.list = list(),i = 1,j =2, p_val=10^-6,which.genes = bmp.receptors, subset_genes = T){
    df = res.list[[i]][[j]];
    if(is.data.frame(df)){

        if(subset_genes ==T){
            df %>% mutate(gene = row.names(df)) %>% filter(gene %in% which.genes ) %>% filter(p_val_adj<= 10^-6) %>% dim() -> d
        }else{
            df %>% mutate(gene = row.names(df)) %>% filter(p_val_adj<= 10^-6) %>% dim() -> d
        }

        return(d[1])


    }else{return(0)}


}



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


plotProfiles<-function(bmp.tidy = bmp.tidy , targets = c(),return.plots =F,ymax = 4){
    p_list = list();
    for(i in 1:length(targets)){
        bmp.tidy %>% filter(cluster==targets[i]) %>%
            ggplot(aes(x=Gene,y = avg_log_expr)) +
                geom_bar(stat = 'identity') +
                ggtitle(paste("Profile ",targets[i])) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                theme(text = element_text(size = 20)) + ylim(0,ymax)->p;
            p_list[[i]] = as.grob(p)
    }
    if(return.plots){
        return(p_list)
    }else{
        grid.arrange(grobs = p_list)
    }
}


inspectProfiles <- function(input.tidy = c(),which.clusters){

    p1 = plotProfiles(bmp.tidy = notch.tidy, targets = which.clusters) ;

    tabula_meta %>% filter(seurat_clusters %in% which.clusters) %>% mutate(Cell_type = paste(cell_ontology_class," (",tissue,")",sep="")) %>% group_by(Cell_type) %>% count() %>% filter(n>10) %>%  ggplot(aes(x=Cell_type,y = n)) + geom_bar(stat = 'identity',fill='darkcyan') + coord_flip() + theme(text =element_text(size = 25)) + xlab("Cell type") -> p2
x11();  grid.arrange(p1,p2, ncol =2)

}
