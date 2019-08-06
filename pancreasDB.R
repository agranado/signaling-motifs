

load()

load("Downloads/SRA653146_SRS2874280.sparse-RPKM.RData")

#minimum percentage of cells in which a gene should be expresed
min.percent = 0.05
#Manually Filtering cells
sm2_filtered=sm2[,log10(Matrix::colSums(sm2>0))>3 & log10(Matrix::colSums(sm2))<6  ]
#filter genes that are expressed in less than 5% of cells in the dataset
sm2_filtered = sm2_filtered[Matrix::rowSums(sm2_filtered>0) >round(dim(sm2_filtered)[2] * min.percent)  & log10(Matrix::rowSums(sm2_filtered)) <6,]

bmp.receptors<-c( "Bmpr1a" ,"Bmpr1b" ,"Acvr1"  ,"Acvrl1" ,"Acvr1b" ,"Tgfbr1" ,"Acvr1c" ,"Acvr2a", "Acvr2b", "Bmpr2" ,"Tgfbr2")

#extract names from Panglao format:
gene.names = do.call(rbind,strsplit(row.names(sm2_filtered),"_"))[,1]
ensembl.names = do.call(rbind,strsplit(row.names(sm2_filtered),"_"))[,2] #second part of the string contains this name
ensembl.names<-do.call(rbind,strsplit(ensembl.names,"\\."))[,1]

# The counts slot can be filled with RPKM values, since it is used for dropout correction
sce<-SingleCellExperiment(list(counts = sm2_filtered))
clusters <- quickCluster(sce, min.size=100)
sce <- computeSumFactors(sce, cluster=clusters)
# Now we can normalize the data
sce<-normalize(sce)

#let's convert to Seurat:
sce.seurat <- as.Seurat(sce)

#Seurat converts _ to - in gene names. Let's split it again
gene.names<-do.call(rbind,strsplit(gene.names,"-"))[,1]



subset_full =sm2_filtered[which(gene.names %in% bmp.receptors) ,]
p = pheatmap(subset_full)



# # #
# Tabula Muris Pancreas


tabula.path ="/home/agranado/MEGA/Caltech/rnaseq/tabula-muris/"

FACS_files = list.files(paste(tabula.path,"00_data_ingest/", "00_facs_raw_data/FACS" ,sep=""), full.names = TRUE)

# Pancreas is item 15

raw.data<-read.csv(FACS_files[15],row.names = 1)
