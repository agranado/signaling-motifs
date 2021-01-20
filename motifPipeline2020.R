# Main script for motifs paper
# Relates to Figure1
# A general pipeline to retrieve sets of genes for a scRNA seq datasets and perform clustering analysis using the silhouette score
library(cluster)
library(stylo)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)
# Lists of pathways
# We define the lists for some pathways baed on literatue knowledge
# For other pathways we used automated databases + manual curation

# tiss.norm is a Seurat object with the scRNA seq data
# The object has to be normalized and analysed previous to the following analysis


### PATHWAYS
get_pathway_genes<-function(seurat_obj = tiss.norm, upper = F){
  # NOTCH
  notch.genes = c("Dll1",   "Dll3"   ,"Dll4",    "Jag1"  , "Jag2", "Notch1", "Notch2", "Notch3", "Notch4", "Mfng",  "Rfng"  , "Lfng")
  # BMP / TGBFb
  bmp.receptors<-c( "Bmpr1a" ,"Bmpr1b" ,"Acvr1"  ,"Acvrl1" ,"Acvr1b" ,"Tgfbr1" ,"Acvr1c" ,"Acvr2a", "Acvr2b", "Bmpr2" ,"Tgfbr2")
  #onlyt BMP
  seven.receptors =c('Acvr1', 'Acvrl1', 'Bmpr1a','Bmpr1b', 'Acvr2a', 'Acvr2b',  'Bmpr2')
  # WNT
  wnt.receptors = c( paste('Fzd',as.character(seq(1,10)),sep=""), 'Lrp5', 'Lrp6')
  wnt.ligands = grep('Wnt.*',row.names(seurat_obj),value = T, ignore.case = T)

  #Eph-ephrin
  #  binding and activation of Eph/ephrin intracellular signaling pathways can only occur via direct cell-cell interaction
  ephrin_receptors = grep("Eph.*",row.names(seurat_obj), value = T, ignore.case = T) %>% sort( )
  ephrin_receptors<-ephrin_receptors[1:14]

  ephrin_ligands = grep("Efn.*",row.names(seurat_obj), value = T,ignore.case = T) %>% sort( )

  #Fgf signaling
  fgf_genes = grep("Fgf.*",row.names(seurat_obj),value = T,ignore.case = T)

  # A positive control
  # cell cycle genes
  cell_cycle = c("Ccnd1","Cdk4","Dkk3","Pax6","Ccnb1","Ccne2","Pcna","Ube2c")
  # sonic
   sonic_genes = c("Shh",   "Ihh" ,  "Dhh" ,  "Ptch1", "Smo" ,  "Gli1" , "Gli2" , "Gli3",  "Sufu")
  # Splicing family: by Ding
  splice_srsf = grep("Srsf.",row.names(seurat_obj),value = T,ignore.case = T)
  # LPA signaling
  LPA = grep("Lpar.",row.names(seurat_obj),value =T,ignore.case = T)
  # Make lists of pathways
  pathway_list = list(notch.genes, bmp.receptors, seven.receptors, wnt.receptors,
                    ephrin_receptors, ephrin_ligands,fgf_genes, splice_srsf, LPA, wnt.ligands)
  all_pathways_names = c("Notch", "Bmp-Tgfb receptors", "BMP 7receptors","Wnt receptors",
                          "Eph receptors", "Eph ligands" ,"Fgfr_LR","Splice_srsf",'LPA_signaling', 'Wnt_ligands')

  if(upper){
    pathway_list = lapply(pathway_list, toupper)
  }
  return(list(pathway_list,all_pathways_names))
}

# # # # # #
### FUNCTIONS

filterMinExp<-function(gene.list = c(),global_stats = c(), min.expr = 0.02){
      pathway_filtered = global_stats[gene.list,] %>% filter(mean> min.expr) %>% pull(gene) %>% as.character()
      return(pathway_filtered)
}
# Function Expression profiles for a list of genes
# First step: calculate the average expression profile for each global cluster
# For a given scRNA object AND a list of genes, it retrieves the average expression profiles for each cluster
# FetchData is using SCT normalization from Seurat
# All expression values used in the paper are SCT normalized. See the Seurat SCT vignette for additional information.
avg.matrix<-function(seurat.obj,genes.list,by= "seurat_clusters",upper.case = T,
                      min.expr = 0.02, min.genes.pr = 0.2, filter_clusters = F){
  cellnames = rownames(seurat.obj@meta.data)
  genenames = rownames(seurat.obj)


  if(upper.case) genes.list = toupper(genes.list)

  genes.list = genenames[which(toupper(genenames) %in% toupper(genes.list))]
  # Active Assay is SCT (you can in principle do it using the RNA slot but that is only Log normalized, though to be fair the values are almost identical)
  data.to.plot = FetchData(seurat.obj, c(genes.list, by))
  data.to.plot$cell = rownames(data.to.plot)
  #we added 2 new fields, get the gene names by excluding them (or get the before)...
  genes.plot = colnames(data.to.plot)[1:(length(colnames(data.to.plot))-2)]
  data.to.plot %>% gather( key =genes.plot, c(genes.plot), value = expression) -> data.to.plot
  data.to.plot %>% dplyr::group_by_at(c(by, "genes.plot")) %>% dplyr::summarize(avg.exp = mean(expression)) %>%
      spread(genes.plot,avg.exp) -> mean.expr.matrix
  mean.mat =as.matrix( mean.expr.matrix[,-1]); rownames(mean.mat)<-unlist(mean.expr.matrix[,by])
  #data.to.plot %>% dplyr::group_by_at(c(by, "genes.plot")) %>% dplyr::summarize(avg.exp = mean(expression),sd.exp = sd(expression))

  return(mean.mat)
}

# Percetage of expressing cells per cluster
# For a given list of genes, this function returns a matrix of pct.exp
# we can then use this matrix to filter data points before clustering
pct.exp.matrix <- function(seurat.obj, genes.list, by ='seurat_clusters', upper.case = F ){

  cellnames = rownames(seurat.obj@meta.data)
  genenames = rownames(seurat.obj)


  if(upper.case) genes.list = toupper(genes.list)


  genes.list = genenames[which(toupper(genenames) %in% toupper(genes.list))]
  # For percentag of expressing cells we use the counts instead of the normalized data.
  data.to.plot = FetchData(seurat.obj, c(genes.list, 'seurat_clusters'), slot = 'counts')
  data.to.plot$cell = rownames(data.to.plot)
  data.to.plot %>% gather(key = genes.plot, genes.list, value = expression) %>%
      group_by(seurat_clusters, genes.plot) %>%
      summarise(pct.exp = sum(expression>0,na.rm = T)/n(), mean = mean(expression)) -> tidy.pathway


  # Spread into matrix form
  tidy.pathway %>% select(-mean)    %>%  spread(genes.plot,pct.exp, fill = 0 ) -> mean.expr.matrix
  mean.mat =as.matrix( mean.expr.matrix[,-1]); rownames(mean.mat)<-unlist(mean.expr.matrix[,by])


  return(mean.mat)

}

# We want to remove clusters that have no significan expression in at least x genes
# This relates to whether we can define a motif with only 1 gene expressed. But also works as a noise filter, since
# profile with general low expression can't be assesed by the DGE matrix
filterClusters<-function(input_mat = c(), min.expr = 0.5,min.gene.expr =2){
  input_rank = rowSums(input_mat>min.expr)
  return(input_mat[which(input_rank >= min.gene.expr),])


}


# Function Transcriptome genes mean expression
# returns a vector with all transcriptome genes mean AND standard deviation
# This will be used to sample random pathways
# filter for low expression necessary to remove genes with almost 0 expression across all clusters (potential source of noise)
#
transcriptomeStats<-function(min.expression = 0.02){
    # Calculate means and var for all genes (for control)
    all_means=apply(tiss.norm[['RNA']]@data,1,mean)
    all_vars =apply(tiss.norm[['RNA']]@data,1,sd)

    keep_genes = all_means > min.expression
    all_means = all_means[keep_genes]
    all_vars = all_vars[keep_genes]
    return(data.frame(mean =all_means,sd = all_vars, gene=names(all_means)))
}

transcriptomeStats<-function(min.expression = 0.02){
  # Precomputed on AWS based on the average count matrix
  # Note: original function computes mean and variance from single cells
  load('/home/agranado/MEGA/Caltech/rnaseq/groupMeeting2020/globalDistance/global_stats_avg_matrix.rda')
  all_means = global_stats_avg_matrix$mean
  all_vars = global_stats_avg_matrix$sd
  all_genes = global_stats_avg_matrix$gene
  # filter min expression on mean
  keep_genes = all_means > min.expression

  all_means = all_means[keep_genes]
  all_vars = all_vars[keep_genes]
  df = data.frame(mean =all_means,sd = all_vars, gene= all_genes[keep_genes] )
  row.names(df) <- df$gene
  return(df)
}

# # # # # #
# Function Clustering scores
# Depends on avg.matrix , randomSameStats
# For a given list of genes, this function computes the silhouette score of the corresponding expresion matrix
# It also computes the silhouette score for a negative control consisting of random sets of genes samples from the transcriptom
# The negative control sets of genes is built to have the same mean and variance of the real pathway

wssPathwayControl<-function(this_pathway = bmp.receptors,scale_matrix = F){
  # before doing kmeans we can scale the input matrix by column such that each gene has mean = 0, sd = 1
  # this prevent a single highly-expressed genes to dominate the clustering and might increase signal from lowly expressed genes.
  # this, however, might create false positives


      n_reps = 100;max.k = 70

      wss = matrix(0,n_reps, max.k )

      #100 times for a random ensemble of genes with the same mean and var as the real pathway
      for(i in 1:n_reps){

          repeat{
          aa = randomSameStats(gene_list = this_pathway)
          if(length(unique(aa))==length(this_pathway)) break
          }
            #print(aa)
          avg_matrix =  if(scale_matrix) scale(avg.matrix(tiss.norm,genes.list = aa, by="seurat_clusters")) else  avg.matrix(tiss.norm,genes.list = aa, by="seurat_clusters")
          wss[i,]= wssKmeans(avg_matrix,k.max = max.k)
      }

      # do it for the real pathway
      #avg_matrix = ifelse(scale_matrix, scale(avg.matrix(tiss.norm,genes.list = this_pathway, by="seurat_clusters")),avg.matrix(tiss.norm,genes.list = this_pathway, by="seurat_clusters"))
      avg_matrix =  if(scale_matrix) scale(avg.matrix(tiss.norm,genes.list = this_pathway, by="seurat_clusters")) else  avg.matrix(tiss.norm,genes.list = this_pathway, by="seurat_clusters")

      aa =wssKmeans(avg_matrix,k.max = max.k)


    #  matplot(t(wss),type ='l',col = alpha('gray',0.4),ylab = "Total within-clusters sum of squares",xlab ="N clusters",xlim = c(0,50),main =pathways[[p]],cex.lab = 2, cex.axis =2 );
    #  lines(aa,col ="red",lwd = 2)

      return(list(wss,aa))

}

silhPathwayControl<-function(this_pathway = bmp.receptors,scale_matrix = F,max.k = 70,
      subset =F, which.clusters = c(), clustering_method = "kmeans", groups = 'seurat_clusters',  n_reps = 300){

      # if(subset){
      #   max.k = length(which.clusters)-5
      # }

      # if(subset){
      #   max.k = sum(which.clusters)-10
      # }
      # do it for the real pathway


      # Starts here:
      # 0. Get the average count matrix
      avg_matrix = avg.matrix(tiss.norm,genes.list = this_pathway, by=groups)
      # Remove cell types with zero values for all genes in pathway
      if(subset){
        which.clusters = rowSums(avg_matrix)> 0
        avg_matrix = avg_matrix[which.clusters, ]
      }
      # Distance matrix for silhouette score:
      # Compute separatley either before/after scaling or using different distance metrics e.g. cosine
      silh_distance = dist(avg_matrix %>% scale() ) # Euclidean distance on the scaled matrix (default)
      #silh_distance = dist.cosine(avg_matrix) # Use cosine distance on the unscaled data


      # scale matrix?
      if(scale_matrix){
          avg_matrix = scale(avg_matrix)
      }


      # 1. calculate the silhouette score for the real data
      silh_score =silhoutteKmeans(avg_matrix,k.max = max.k, clust_method = clustering_method, distance_matrix = silh_distance)


      wss = matrix(0,n_reps, max.k-1)
      # Return also the list of random pathways, for further analysis
      random_pathways = list()

      # 2. 100 times for a random ensemble of genes with the same mean and var as the real pathway
      for(i in 1:n_reps){
          repeat{
          aa = randomSameStats(gene_list = this_pathway)
          if(length(unique(aa))==length(this_pathway)) break
          }
            #print(aa)
          avg_matrix = avg.matrix(tiss.norm,genes.list = aa, by=groups)


          if(subset){
              avg_matrix = avg_matrix[which.clusters, ]
          }

          # Oct 6th 2020
          # Compute the distance metric for silhouette score
          silh_distance = dist(avg_matrix %>% scale()  ) # Euclidean distance on the scaled matrix (default)
          #silh_distance = dist.cosine(avg_matrix) # Use cosine distance on the unscaled data

          avg_matrix =  if(scale_matrix) scale(avg_matrix) else  avg_matrix

          # Here we can subset the matrix to filter for clusters that have enough levels of expression to call it a expression profile
          # Whatever clusters we choose here (based on the real pathway profile) we will calculate the control for the same cluster
          # This means that sometimes the selected clusters will have low expression for the random pathways (since we filtered only based on the real)

          wss[i,]= silhoutteKmeans(avg_matrix,k.max = max.k,clust_method = clustering_method, silh_distance)

          random_pathways[[i]] = aa
      }



      #matplot(t(wss),type ='l',col = alpha('gray',0.4),ylab = "Silhouette score",xlab ="N clusters",xlim = c(0,50),main =pathways[[p]],cex.lab = 2, cex.axis =2 );
      #lines(aa,col ="red",lwd = 2)
      # 3. Return both control and real calculation
      return(list(wss,silh_score, random_pathways))

}

# Spectrum analysis for random sets of genes:
# Oct 6th 2020
spectrumAnalysis <- function(this_pathway=bmp.receptors, spec_method = 2, max.k = 70, tune_kernel = F, scale_matrix = F, n_reps=100, subset = F, groups = 'seurat_clusters'){
  # Expression matrix across cell types
  x = avg.matrix(tiss.norm, this_pathway, by ='seurat_clusters')

  # 0. Remove cell types with zero values for all genes in pathway
  if(subset){
    which.clusters = rowSums(x)> 0
    x = x[which.clusters, ]
  }
  # scale matrix?
  # Scaling takes place after filtering..
  # Note: 10052020 When we filter cell types, we are removing very low expression (close to zero) Do we want the z-score to consider those cell types though?
  # Basically by filtering, we are incresing the floor of lowest expression
  if(scale_matrix){
      x = scale(x)
  }


  # 1. Calculate Spectrum result for the real pathway
  pathway_res = Spectrum(t(x), method = spec_method, showpca = F, maxk =max.k, dotsize = 2, tunekernel = tune_kernel, silent = T)


  # 2. Calculate random pathways
  random_res_list = list()
  for(i in 1:n_reps){
      repeat{
      aa = randomSameStats(gene_list = this_pathway)
      if(length(unique(aa))==length(this_pathway)) break
      }
        #print(aa)
      avg_matrix =  if(scale_matrix) scale(avg.matrix(tiss.norm,genes.list = aa, by=groups)) else  avg.matrix(tiss.norm,genes.list = aa, by=groups)

      # Here we can subset the matrix to filter for clusters that have enough levels of expression to call it a expression profile
      # Whatever clusters we choose here (based on the real pathway profile) we will calculate the control for the same cluster
      # This means that sometimes the selected clusters will have low expression for the random pathways (since we filtered only based on the real)

      if(subset){
        avg_matrix = avg_matrix[which.clusters, ]
      }

      random_res_list[[i]] = Spectrum(t(avg_matrix),method = spec_method, showpca = F, maxk =max.k, dotsize = 2, tunekernel = tune_kernel, silent = T, showres = F)
      # Spectrum fixes the random seed which creates a problem when sampling random sets of genes

      # We can not use random number to set the seed because the random.seed is locked!
      # This sets a new random seed using clock data (?) which is independent of the curret random seed
      set.seed(NULL)
      set.seed(round(runif(1,0,1)*65000))
      #wss[i,]= silhoutteKmeans(avg_matrix,k.max = max.k,clust_method = clustering_method)

  }

  return(list(pathway_res, random_res_list))


}

# Based on the Differential expression matrix, we can estimate the optimal k using spectrum
# Oct 6th 2020
estimateK_spectrum <-function(this_pathway=bmp.receptors, res.list_ = list(), p_val_threshold = 10^-6,
                              which.method = 'fold', global_stats = c() , n_reps = 100){

  # estimate DGE matrix for the real pathway
  dge.mat = distMatDEsingle(res.list_, p_val = p_val_threshold ,
                           which.genes = this_pathway ,
                           cluster_list  = 0:(n_clusters-1),
                           which.method = which.method )

  # make it a similarity matrix, instead of distance
  dge.mat_sim = if(which.method =='cosine') 1-dge.mat else (max(dge.mat) - dge.mat)/max(dge.mat)

  k_pathway = estimate_k(dge.mat_sim)


  # 2. Now we can subset the matrix and calculate k for each subset
  # Should be faster than comuting the distance matrix each time
  #
  random_res_list = list()
  for(i in 1:n_reps){
      repeat{
      aa = randomSameStats(gene_list = this_pathway)
      if(length(unique(aa))==length(this_pathway)) break
      }

      # Random pathways
      # 2. estimate DGE matrix for the transcriptome
      dge.mat.rand = distMatDEsingle(res.list_, p_val = p_val_threshold ,
                               which.genes = aa ,
                               cluster_list  = 0:(n_clusters-1),
                               which.method = which.method )

      dge.mat_sim = if(which.method =='cosine') 1-dge.mat.rand else (max(dge.mat.rand) - dge.mat.rand)/max(dge.mat.rand)


      random_res_list[[i]] = estimate_k(dge.mat_sim)

      set.seed(NULL)
      set.seed(round(runif(1,0,1)*65000))
  }
  return( list(k_pathway, random_res_list))
}






# # # # # # # Jan 14th 2020
# Function Calculate clustering scores for a given expression matrix
# clustering pathways based on average cluster expression
# Do we see a lower dimension compared with random sets of genes

wssKmeans<-function(input_matrix,k.max= 50 ){
  #set.seed(123)
  # Compute and plot wss for k = 2 to k = 15.


  wss <- sapply(1:k.max,
                function(k){kmeans(input_matrix, k, nstart=50,iter.max = 30 )$tot.withinss})
  return(wss)
  # plot(1:k.max, wss,
  #      type="b", pch = 19, frame = FALSE,
  #      xlab="Number of clusters K",
  #      ylab="Total within-clusters sum of squares")

}
silhoutteKmeans<-function(input_matrix,k.max, clust_method = 'kmeans', distance_matrix = c() ){

  if(clust_method =='kmeans'){
  wss <- sapply(2:k.max,
                function(k){res = kmeans(input_matrix, k, nstart=50,iter.max = 30 );ss =silhouette(res$cluster,distance_matrix);return(mean(ss[,3]))})
  }else if(clust_method=='pam'){
  wss <- sapply(2:k.max, function(k){res = pam(input_matrix, k );ss =silhouette(res$clustering,dist(input_matrix));return(mean(ss[,3]))})

  }else if(clust_method =="spectral"){
  wss <- sapply(2:k.max,
                  function(k){res = spectralKmeans(input_matrix, k );ss =silhouette(res$cluster,dist(input_matrix));return(mean(ss[,3]))})
              }else if(clust_method =="cosine_kmeans"){
  wss <- sapply(2:k.max,
                    function(k){res = kmeans(dist.cosine(input_matrix) %>% as.matrix(), k );ss =silhouette(res$cluster,dist(input_matrix));return(mean(ss[,3]))})

              }else if(clust_method=="cutree"){
                # This function calculates a hierarchical clustering using cosine distance
                # Then makes partitions of the resulting tree for different values of k
                # The silhouette score is calculated using Euclidean distance
                # p = pheatmap(input_matrix, clustering_distance_rows = dist.cosine(input_matrix), silent = T)
                p = pheatmap(input_matrix, ,silent = T)
                wss <- sapply(2:k.max,
                    function(k){res = treeClust(p, k );ss =silhouette(res,distance_matrix);return(mean(ss[,3]))})

              }else if(clust_method =="cutree_cosine"){
                  p = pheatmap(input_matrix, silent = T, clustering_distance_rows = dist.cosine(input_matrix))
                  wss <- sapply(2:k.max,
                    function(k){res = treeClust(p, k );ss =silhouette(res,distance_matrix);return(mean(ss[,3]))})

              }else if(clust_method=="cutree_cor"){
                  p = pheatmap(input_matrix, silent = T, clustering_distance_rows = as.dist(1-cor(t(input_matrix))))
                  wss <- sapply(2:k.max,
                  function(k){res = treeClust(p, k );ss =silhouette(res,distance_matrix);return(mean(ss[,3]))})

              }

  return(wss)
}

treeClust<-function(p,k){
  #p2 = pheatmap(input_matrix, clustering_distance_rows = dist.cosine(input_mat), silent = T)
  aa = as.hclust(p$tree_row)
  return(cutree(aa, k))
}


# Load ana analyze all pathways in PathBank
# Subset for pathways in Signaling class
# May 15th 2020
# This generates a list of list of genes that we can then run in parallel for silhouette score
loadPathBank<-function(){
  proteins = read.csv("/home/agranado/Downloads/pathbank_all_proteins.csv")
  proteins %>% filter(Species =="Mus musculus") -> proteins.mouse


  proteins.mouse %>% filter(grepl('.*Signaling.*',Pathway.Name) | grepl('.*Signalling.*',Pathway.Name)) -> mouse.signaling


  mouse.signaling %>% group_by(Pathway.Name) %>% count() %>% select(Pathway.Name) -> mouse.pathway.names
  pathway_names = mouse.pathway.names$Pathway.Name
  signaling_genes_list = list()
  for(i in 1:length(pathway_names)){
    this_pathway =pathway_names[i]
    mouse.signaling %>% filter(Pathway.Name == this_pathway) %>% pull(Gene.Name) %>% as.character() -> signaling_genes_list[[i]]
  }

  return(signaling_genes_list)
}


# Function Create negative control pathways
# For a given list of genes, this function will randomly sample the transcriptom
# Return a list with the same number of genes
randomSameStats<-function(gene_list = wnt.receptors, all_stats = global_stats){
  path_means = all_stats[gene_list,]$mean;
  path_vars = all_stats[gene_list,]$sd

  all_means = global_stats$m
  names(all_means)<-global_stats$gene
  all_sd = global_stats$sd
  names(all_sd) = global_stats$gene

  s = 0.1
  rand_set = c()
  for(i in 1:length(gene_list)){
    # the range of mean/variance to sample from, s, should be proportional to the expression of the gene
    if(path_means[i]<= 0.3){
      s = 0.05
    }else if(path_means[i]<=1.4){
      s = 0.1
    # For highly-expresed genes, there are not too many option to sample from
    # and 0.1 makes the range too narrow. So we use a range proportional to the mean of that gene
    }else{
      s = 0.2*path_means[i]
    }

    same_means = which(all_means>=path_means[i]-s & all_means<=path_means[i]+s)
    same_vars = which(all_sd>=path_vars[i]-s & all_sd<=path_vars[i]+s)
    same_stats = same_means[which(same_means %in% same_vars)]
    rand_set[i]=sample(names(same_stats),1)
  }
 return(rand_set)
}


### PLOTS
#Jan 2020 using ggplot for multiple pathways
ggplotPathway<-function(results,type = 'silh', pathway_name =''){

  pathway.data = results[[2]]
  rand.ensemble = results[[1]]


  if(type =='silh'){
    low.lim = 2
    high.lim = length(pathway.data)+1
    score.type ="Silhouette"
    legend.pos = 'none'
  }else if(type =="wss"){
    low.lim = 1
    high.lim = length(pathway.data)
    score.type = "WSS"
    legend.pos = 'none'
  }else if(type =="pca"){
    low.lim = 1
    # number of principal components
    high.lim = length(pathway.data)

    rand.ensemble = t(apply(results[[1]],1,cumsum))
    pathway.data = cumsum(pathway.data)
    score.type = "PCA"
    legend.pos = c(0.7,0.2 )
  }else if(type =='MI'){
    score.type = 'MI'

  }

  #repeats
  row.names(rand.ensemble)<-as.character(1:dim(rand.ensemble)[1])
  # K
  colnames(rand.ensemble)<-as.character(low.lim:high.lim)

  aa<-as.data.frame(rand.ensemble)
  bb= data.frame( k = low.lim:high.lim, m_score = pathway.data,sd_score = rep(0,length(pathway.data)),data = rep("Pathway",length(pathway.data)))

  control.df<-gather(aa, "k","score") %>% group_by(k) %>%
    summarise(m_score = mean(score),sd_score = sd(score)) %>% mutate(k = as.numeric(k)) %>%
      mutate(data = rep("Control",length(pathway.data)))

  control.df %>% rbind(bb) %>% ggplot(aes(x = k, y = m_score, color =data)) +
    geom_ribbon(data=control.df, aes(ymin = m_score - sd_score, ymax=m_score + sd_score),alpha = 0.2,color =NA) +
      geom_line(size = 0.5) + theme_classic() + scale_colour_manual(values=c("black", "deeppink3")) + theme(text =element_text(size =8)) +
        xlab("N clusters") + ylab(paste("Clustering Score (",score.type, ")",sep="")) +
        theme(legend.position= legend.pos) +
        theme(axis.text.x=element_text(size=10), axis.text.y = element_text(size =10)) +
        ggtitle(pathway_name)  -> p

  if(type=="wss"){
    p = p +   scale_y_continuous(trans='log10') +coord_cartesian(ylim=c(100,1000),xlim = c(2,20))
  }else if(type=="silh"){
    p = p + coord_cartesian(xlim=c(3,high.lim))
  }
  return(p)
}


ggplotWithHeatmap<-function(list_results, i, names = all_pathways_names, genes =pathway_list,
                            min.ax = 4, max.ax = 60, min.ax.y = 0, max.ax.y = 0.35, heatmap_ann = data.frame(),
                            groups = 'seurat_clusters', min.expr = 0, cut.k = 10, plot_names = c(), quantile_norm =F, palette_breaks = 10){

  # Plot in the same order as the input list
  results_1  = list_results[[1]][[i]]
  results_2 =   list_results[[2]][[i]]
  results_3 =   list_results[[3]][[i]]
  results_4 =   list_results[[4]][[i]]
  results_5 =   list_results[[5]][[i]]


  #p1 =  ggplotPathway(results_pca,type = "pca")
  p1 = ggplotPathway(results_1,type = "silh", pathway_name = paste(names[i],plot_names[1])) +  coord_cartesian(xlim= c(min.ax,max.ax), ylim = c(min.ax.y,max.ax.y))
  p2 = ggplotPathway(results_2,type = "silh",pathway_name = paste(names[i],plot_names[2])) + coord_cartesian(xlim= c(min.ax,max.ax))

  p3 =  ggplotPathway(results_3,type = "silh",pathway_name = paste(names[i],plot_names[3])) +  coord_cartesian(xlim= c(min.ax,max.ax), ylim = c(min.ax.y,max.ax.y))
  p4 =  ggplotPathway(results_4,type = "silh",pathway_name = paste(names[i],plot_names[4])) +  coord_cartesian(xlim= c(min.ax,max.ax), ylim = c(min.ax.y,max.ax.y))
  p5 =  ggplotPathway(results_5,type = "wss",pathway_name = paste(names[i],plot_names[5])) #+  coord_cartesian(xlim= c(min.ax,max.ax), ylim = c(min.ax.y,max.ax.y))

  mat = avg.matrix(tiss.norm, genes[[i]],by=groups)

  #  p = pheatmap(mat, clustering_distance_rows = dist.cosine(mat), clustering_distance_cols = as.dist(1-cor(mat)), silent = T, fontsize = 11)

  p = exprMatrixAsAnnotatedHeatmap(mat, annotation_df = heatmap_ann, min.expr = min.expr, return_plot = T, silent_plot = T, k = cut.k, quantile = quantile_norm, n_breaks = palette_breaks)

  lay =rbind(c(1,1,1,1,2,5),
             c(1,1,1,1,3,6),
             c(1,1,1,1,4,7))
  p = grid.arrange(grobs = list(p[[4]],p1,p2,p3,p4, p5), layout_matrix = lay, top =textGrob(names[i],gp=gpar(fontsize=20,font=3)))
}




### ### MAKE lists

# Our lists of genes contain manual curation of pathway components
# In the dataset, however, we  don't know if those genes are expressed at highlevels.
# We apply a filter for minimum expression.
#
global_stats = transcriptomeStats()
global_stats = transcriptomeStats2()
gc()

#Generate all pathways list
pathway_info = get_pathway_genes()
pathway_list = pathway_info[[1]]
all_pathways_names = pathway_info[[2]]

# Apply min expression filter:
pathway_list = lapply(pathway_list, filterMinExp, global_stats = global_stats, min.expr = 0.02)

# Load all signaling genes in the PathBank database
# Loads a csv file available for download at PathBank
signaling_genes_list = loadPathBank()
signaling_genes_list = lapply(signaling_genes_list, filterMinExp, global_stats = global_stats, min.expr = 0.02)
gc()
