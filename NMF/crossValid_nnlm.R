
# NNLM package 4Apr


# assuming data loaded already

set.seed(678);


k <- 10;
W  <- matrix(runif(n*k),  n, k);
H  <- matrix(10*runif(k*m),  k, m);

n <- dim(norm.bmp)[1];
m <- dim(norm.bmp)[2];
noise <- matrix(rnorm(n*m), n, m);
A <- norm.bmp + 1/2*noise;
A[A < 0] <- 0;
A = as.matrix(A)
max.k = 20
plot(-1, xlim = c(1,max.k), ylim = c(0.5, 1000), xlab = "Rank", ylab = "MSE")
cols <- c('deepskyblue', 'orange', 'firebrick1', 'chartreuse3',"black",'chocolate','darkred','lightpink',
          'gold','brown4','grey100','navy','yellow','darkgreen','gold2','gold3','cadeblue','gray21');


all.err  = list()
for (i in 1:100) {
    col = cols[i]
    ind <- sample(length(A), 0.3*length(A));
    A2 <- A;
    A2[ind] <- NA;
    err <- sapply(X = 1:max.k,
        FUN = function(k) {
            z <- nnmf(A2, k,check.k = F,n.threads = 8,loss = "mkl");
            c(mean((with(z, W %*% H)[ind] - A[ind])^2), tail(z$mse, 1), tail(z$mkl,1));
            }
        );
    #invisible(lines(err[1,], col = col, type = 'b'));
    #invisible(lines(err[2,], col = col, type = 'b', lty = 2));
    all.err[[i]] = err
}

all.err.mat =simplify2array(all.err)

#gain in MSE as you add higer rank to factorization:
 plot(diff(apply(all.err.mat[2,,],1,mean,na.rm = T)))



# # # # # #
library(doParallel)
registerDoParallel(cores = 6)
getDoParWorkers()

all.err  = list()
foreach (i in 1:length( cols),  .export =c('nnmf'),.packages='NNLM') {
    col = cols[i]
    ind <- sample(length(A), 0.3*length(A));
    A2 <- A;
    A2[ind] <- NA;
    err <- sapply(X = 1:max.k,
        FUN = function(k) {
            z <- nnmf(A2, k,check.k = F);
            c(mean((with(z, W %*% H)[ind] - A[ind])^2), tail(z$mse, 1));
            }
        );
    invisible(lines(err[1,], col = col, type = 'b'));
    invisible(lines(err[2,], col = col, type = 'b', lty = 2));
    all.err[[i]] = err
}
