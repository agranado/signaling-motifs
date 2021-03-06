
library(Matrix)
library(Seurat)
#data from intestine single cell rna seq
test <-read.table("table_B_scRNAseq_UMI_counts.tsv",sep="\t",header=T)
#gene names are in 1st column
test2<-test[,-1]
rownames(test2)<-test[,1]

raw.data<-Matrix(as.matrix(test2),sparse =T)

#plot the histogram of reads across cells:
hist(log10(Matrix::colSums(raw.data)))


raw.data<-raw.data[keep_feature,]

#remove mitocondrial genes
#in this dataset, cells with high mitochondrial content (>50%) have been already removed
# this matrix is basically already being filtered and QC'd

mit <- grep(pattern = "^mt-", x = rownames(x = raw.data), value = TRUE)

#percent is the ratio, for each cell, between the sum of ERCC detection divided by the total count
percent.mit <- Matrix::colSums(raw.data[mit, ])/Matrix::colSums(raw.data)
sum.mit<-Matrix::colSums(raw.data[mit, ])
mit.index <- grep(pattern = "^mt-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-mit.index,] #remove the ERCC sequences


bmp.all = pathway.genes(pathway = "bmp")
bmp.indexes=match(bmp.all,rownames(raw.data))




#----------------------
#  Using SingleCellExperiment
umi<-SingleCellExperiment(
     assays = list(counts = as.matrix(test2))
 )
keep_feature<- rowSums(counts(umi)>0) >0
umi <-umi[keep_feature,]

isSpike(umi,"MT")<-grepl("^mt-",rownames(umi))

umi<-calculateQCMetrics(umi, feature_controls = list(MT=isSpike(umi,"MT")))
# Number of reads per cell
hist(umi$total_counts,breaks=100)
# Number of genes detected per cell
hist(umi$total_features,breaks=100)
#to find outliers in the data
umi <- plotPCA(
     umi,
     size_by = "total_features",
     pca_data_input = "pdata",
     detect_outliers = TRUE,
     return_SCE = TRUE
 )
