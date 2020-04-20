
library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(readr)
library(here)


tiss_dirs = grep('senis',list.dirs(recursive = F),value = T)

readTabulaMurisSenis<-function(tiss_dir = ""){

  # it is transposed as per scanpy convention
  raw.data = read.csv(paste(tiss_dir,'/X.csv',sep=""),header = F)
  raw.data = t(raw.data)

  # Meta data for cells:
  # Includes age
  meta_data = read.csv(paste(tiss_dir,'/obs.csv',sep=""))
  colnames(raw.data) = meta_data$index

  # Meta data for genes:
  gene_data = read.csv(paste(tiss_dir,'/var.csv',sep=""))
  row.names(raw.data)= gene_data$index

  # Select just for age 3m (We want a replicate of the original Tabula Muris, this is as close as it gets)
  # meta_data %>% filter(age=='3m' | age=="1m") -> selected_cells
  # raw.data_filtered = raw.data[,selected_cells$index]

  raw.data <- Matrix(as.matrix(raw.data), sparse = TRUE)

  return(raw.data)
}

readTabulaMurisSenisMeta<-function(tiss_dir = ""){
  # Meta data for cells:
  # Includes age
  meta_data = read.csv(paste(tiss_dir,'/obs.csv',sep=""))
  # Select just for age 3m (We want a replicate of the original Tabula Muris, this is as close as it gets)
  return(meta_data)
}



#
# i = 1
# raw.data.list = list()
# meta.data.list = list()
#
# for(i in 1:length(tiss_dirs)){
#     # has no col or row names.
#     # it is transposed as per scanpy convention
#     raw.data = read.csv(paste(tiss_dirs[i],'/X.csv',sep=""))
#     raw.data = t(raw.data)
#
#     # Meta data for cells:
#     # Includes age
#     meta_data = read.csv(paste(tiss_dirs[i],'/obs.csv',sep=""))
#     colnames(aa) = meta_data$index
#
#     # Meta data for genes:
#     gene_data = read.csv(paste(tiss_dirs[i],'/var.csv',sep=""))
#     row.names(raw.data)= gene_data$index
#
#     # Select just for age 3m (We want a replicate of the original Tabula Muris, this is as close as it gets)
#     meta_data %>% filter(age=='3m') -> selected_cells
#     raw.data_filtered = raw.data[,selected_cells$index]
#
#     raw.data <- Matrix(as.matrix(raw.data), sparse = TRUE)
#
#     raw.data.list <- append(raw.data.list, raw.data)
#     meta.data.list<- append(meta.data.list, selected_cells)
# }

raw.data.list = mclapply(tiss_dirs, readTabulaMurisSenis)

# they have different N of genes for some reason
# I need to either drop some tissues with low N of genes
# OR find intersection. Seems that it is only Bladder
# Ignore Bladder for now:
# next organ with less genes: Mammary glad. We can use this as our list of genes (17,232 genes)
raw.data.list[[11]] %>% row.names() -> mammary_genes
# drop Bladder (for now )
raw.data.list = raw.data.list[-1]
# Find the intersection for gene names
# Not all organs have the 17k genes from mammary glad
all_genes<-lapply(raw.data.list, row.names)
all_genes = Reduce(intersect, all_genes)

# filter the intersect-genes for all organ samples
raw.data.list.fil = lapply(raw.data.list, function(x){ return(x[all_genes,]) } )
# 17k genes, 62k cell
raw.data = do.call(cbind,raw.data.list.fil)



# READ META DATA
meta.data.list = mclapply(tiss_dirs, readTabulaMurisSenisMeta,mc.cores = 8)
# For some reason meta data for different tissues has different N of attributes
# Find the intersection
unique_meta_names = Reduce(intersect, lapply(meta.data.list, names) )
# filter for these attributes


meta.data.list_uniqueCols<-lapply(meta.data.list, function(x){ x %>% select(unique_meta_names) })
metaData_senis = do.call(rbind,meta.data.list_uniqueCols)

# subset cells and gens for final matrix

metaData_senis %>% filter(tissue != "Bladder") %>% filter(age=="24m" | tissue == "Mammary_Gland") -> selected_cells
raw.data = raw.data[,selected_cells$index]


#Seurat
row.names(selected_cells)<-selected_cells$index

tiss<-CreateSeuratObject(counts = raw.data,meta.data = selected_cells)
