# Correaltion between Senis and Tabula Muris datasets:

# Given two expression matrices X and Y (average across global clusters)
# Calculate the correlation of each row (cell type) in Y against all the profiles in X (the reference matrix)
# Count the number of high confidence matches ( > cor_threshold) and create a rank for the profiles in Y
# Returns a vector of length = dim(Y)[1]
compareAtlasProfileCorr<-function(expr.ref = bmp.mat.tm, expr.mat = bmp.mat.18m , cor_threshold = 0.8){

  # re order columns for consistency
  expr.mat  = expr.mat[,colnames(expr.ref)]
  matches = c()
  # for each profile in the test matrix, we calculate the correlation against ALL
  # profiles in the reference matrix and count how many reference profiles have a correspondence (cor > threshold)
  for(i in 1:dim(expr.mat)[1]){
    matches[i]  =  sum(apply(expr.ref,1,cor,expr.mat[i,]) >= cor_threshold )

  }
  return(matches)

}

# For a given pathway runs compareAtlasProfilesCor between:
# Tabula Muris and 18m
# Tabula Muris and 24m
# Tabula Muris and Tabula Muris
# CONTROL: Tabula Muris vs randomized versions of the three matrices
correlationRankPlot<-function(ref.mat =bmp.mat.tm, mat18 = bmp.mat.18m, mat24 =bmp.mat.24m , cor_min = 0.85,
                              thresholds = c(2,1.5,1)){


  # remove cell types that do not significantly express components from the pathway:
  # Sort colanmes as those in ref.mat for consistency 
  bmp.mat.tm.fil = ref.mat[ rowSums(ref.mat)>=thresholds[1] ,colnames(ref.mat)]
  bmp.mat.18m.fil = mat18[ rowSums(mat18)>=thresholds[2] ,colnames(ref.mat)]
  bmp.mat.24m.fil = mat24[ rowSums(mat24)>=thresholds[3] ,colnames(ref.mat)]

  # Create randomized versions of the 3 matrices
  rand.tm = apply(bmp.mat.tm.fil,1,sample) %>% t()
  rand.18m = apply(bmp.mat.18m.fil,1,sample) %>% t()
  rand.24m = apply(bmp.mat.24m.fil,1,sample) %>% t()

  colnames(rand.tm)<-colnames(ref.mat)
  colnames(rand.18m)<-colnames(ref.mat)
  colnames(rand.24m)<-colnames(ref.mat)

  # Calculate the correlations for each profile in the target_matrix
  # Reference matrix is always the Tabula Muris (3 month old mouse)

  compareAtlasProfileCorr(expr.ref = bmp.mat.tm.fil, expr.mat  = bmp.mat.18m.fil, cor_threshold = cor_min) %>% sort(decreasing = T) %>% plot(type = "l", lwd = 2, xlab = "Cell type rank",ylab="Number of matches in Cell Atlas",main="Correlation analysis across Cell Atlas datasets")
  compareAtlasProfileCorr(expr.ref = bmp.mat.tm.fil, expr.mat  = bmp.mat.24m.fil, cor_threshold = cor_min) %>% sort(decreasing = T) %>% lines(type = "l", lwd = 2,col="blue")
  compareAtlasProfileCorr(expr.ref = bmp.mat.tm.fil, expr.mat  = bmp.mat.tm.fil, cor_threshold = cor_min) %>% sort(decreasing = T) %>% lines(type = "l", lwd = 2,col="gold3")

  # Calculate the correlation rank for the randomized data (as a control for the expected N of matches for a random profile)
  # Plot the randomized versions
  compareAtlasProfileCorr(expr.ref = bmp.mat.tm.fil, expr.mat  = rand.18m, cor_threshold = cor_min) %>% sort(decreasing = T) %>% lines(type = "l", lwd = 2,col = "grey")
  compareAtlasProfileCorr(expr.ref = bmp.mat.tm.fil, expr.mat  = rand.24m, cor_threshold = cor_min) %>% sort(decreasing = T) %>% lines(type = "l", lwd = 2,col="lightblue3")
  compareAtlasProfileCorr(expr.ref = bmp.mat.tm.fil, expr.mat  = rand.tm, cor_threshold = cor_min) %>% sort(decreasing = T) %>% lines(type = "l", lwd = 2,col="lightgoldenrod3")


}
