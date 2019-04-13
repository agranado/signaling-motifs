

library(ccfindR)
library(pheatmap)
library(parallel)
#This factorization was calculated using the log-norm. It seems that ccfindR does not use normalized data by default
#which is good because I don't want the data to be normalized twice

load("bmp_BayesFactr_normLog_12kcells.rda")

  #for each rank value (3:17) we can calcualate the kmeans elbow Analysis
  # check which elements are non-0
  max.rank = sum(unlist(lapply(bmp.nmf_@coeff,length))>0)
  #max.rank = 3
  # kmeans parameters:
  #
  k.max = 20
  data.matrix = matrix(0,max.rank,k.max)

  no_cores= detectCores() -1

  for(r in 1:max.rank){
      H_trans = t(bmp.nmf_@coeff[[r]])
      cl <- makeCluster(no_cores,type = "FORK")

      #wss<-sapply(1:k.max,function(k){kmeans(H_trans,k,nstart=50,iter.max =15)$tot.withinss} )
      #parallel version : returns a list
      wss = parLapply(cl, 1:k.max, function(k){kmeans(H_trans,k,nstart=50,iter.max =30)$tot.withinss} )

      data.matrix[r,] = unlist(wss)
      stopCluster(cl)
  }




  
