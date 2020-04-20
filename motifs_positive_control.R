#positive control plots for motifs
# Goal: how do the silhouette and WSS scores look when we KNOW the actual structure of the data (in terms of motifs)


m = 5  # number of motifs in the dataset
c = 80 # number of cell types
g = 10 # or the same as a real pathway
input_mat<-mat.receptors <-avg_all_82Clusters[,bmp.receptors]
chooseMotifs<-function(input_mat, g = dim(input_mat)[2], c = dim(input_mat)[1],m=5, pm =rep(1/m,m)){
    #sample the real matrix for m motifs
    mat.motifs<- input_mat[sample(1:dim(input_mat)[1],m),]

    # Each cell type picks a motif with probability pm
    c_id = c()
    for(i in 1:c){
        r = runif(1)
        sum_pr = 0
        for(p in 1:length(pm)){
            sum_pr = sum_pr + pm[p]
            if(r<sum_pr) {m_id = p; break}
        }
        c_id[i] = m_id
    }
    positive.matrix = matrix(0,c,g)

    # create matrix with motifs
    for(i in 1:dim(positive.matrix)[1]){
        positive.matrix[i,] = mat.motifs[c_id[i],]
    }

    return(positive.matrix)
}

addNoiseMotifs<-function(positive.matrix =c(), reference_matrix,scale_var = 10){

  noise_matrix = do.call(cbind,lapply(apply(reference_matrix,2,sd)^2/scale_var, rnorm,n = dim(positive.matrix)[1], mean = 0))
  noisy_matrix = positive.matrix + noise_matrix
  noisy_matrix[noisy_matrix<0] = 0

  return(noisy_matrix)
}

# run it once
# run multiple times with :  replicate(100, runEnsemblePositive(gene_matrix = mat.receptors, m = 5, scale_var = 10,type="pca" ))
runEnsemblePositive <- function(gene_matrix = mat.receptors, m = 8 , scale_var = 10, type = "silh" ){



  chooseMotifs(gene_matrix,m = m) %>%
    addNoiseMotifs(reference_matrix = gene_matrix, scale_var = scale_var) -> test.mat

    if(type=="silh"){
      test.mat %>% silhoutteKmeans(k.max = m*2) -> x
    }else if(type =="wss"){
      test.mat %>% wssKmeans(k.max = m*2) -> x
    }else if(type =="pca"){
      test.mat %>% pcaStatKmeans() -> x
    }else if(type =="all"){
      x = list()
        test.mat %>% silhoutteKmeans(k.max = m*2) -> x[[1]]
        test.mat %>% wssKmeans(k.max = m*2) -> x[[2]]
        test.mat %>% pcaStatKmeans() -> x[[3]]
    }

  return(x)
}

makePlotEnsemble <-function(input_matrix, type ="silh"){


   if(type =='silh'){
     low.lim = 2
     high.lim = dim(input_matrix)[1]+1
     score.type ="Silhouette"
     legend.pos = 'none'
     x_lab="N clusters"
   }else if(type =="wss"){
     low.lim = 1
     high.lim = dim(input_matrix)[1]
     score.type = "WSS"
     legend.pos = 'none'
     x_lab  ="N clusters"
   }else if(type =="pca"){
     low.lim = 1
     # number of principal components
     high.lim =dim(input_matrix)[1]

     input_matrix = apply(t(input_matrix),1,cumsum)
     score.type = "PCA"
     legend.pos = c(0.7,0.2 )
     x_lab ="N PCs"
   }

  input_matrix %>% t() %>% as.data.frame() -> pos.df

   names(pos.df)<-low.lim:high.lim

   pos.df %>% gather("k","score") %>% group_by(k) %>%
    summarise(m_score = mean(score),sd_score = sd(score))  %>%
      mutate(k = as.numeric(k)) %>% ggplot(aes(x = k , ,y=m_score ,group=1)) + geom_line() +
        geom_ribbon( aes(ymin = m_score - sd_score, ymax=m_score + sd_score),alpha = 0.2 )  +
          theme(text = element_text(size = 20)) + xlim(1,high.lim) + xlab(x_lab) + ylab(paste("Clust Score (",score.type,")", sep="")) -> p

  return(p)
}
