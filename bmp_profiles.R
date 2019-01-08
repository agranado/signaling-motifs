library(tidyverse)
library(stringr)
library(Seurat)
library(here)

## ------------------------------------------------------------------------
load(file=here("00_data_ingest", "11_global_robj", "FACS_all.Robj"))
#object name is tiss_FACS and it is fully annotated by the tabula muris people (after running the main script)

bmps<-grep(pattern = "Bmp", x = rownames(x = tiss_FACS@data), value = TRUE)
#extract the data for all bmp genes
bmp.counts<-tiss_FACS@data[bmps,]
#type 1
bmp.receptors<-c("Bmpr1a","Bmpr1b","Acvr1","Acvrl1","Acvr1b","Tgfbr1","Acvr1c","Acvr2a","Acvr2b","Bmpr2","Tgfbr2")
bmp.ligands<-c("Bmp2","Bmp3","Bmp4","Bmp5","Bmp6","Bmp7","Bmp10","Bmp15",
            "Bmp8a","Gdf2","Gdf1","Gdf3","Gdf5","Gdf6","Gdf7","Gdf9","Gdf10","Gdf11","Gdf15")


plots.folder = "/home/agranado/MEGA/Caltech/rnaseq/bmp-profiles"

#plot the bmp receptors on top of the tSNE created for the Nature paper
#this tSNE used ALL genes and 43 principal components
x11();FeaturePlot(tiss_FACS,bmp.receptors,cols.use = c("lightgrey","blue"),nCol=4)
ggsave((plots.folder,"bmpReceptors_tSNE.pdf",sep="")  , width = 14, height = 7, units = "in")


# # # # https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html
# # # # https://satijalab.org/seurat/conversion_vignette.html
# # # # https://satijalab.org/seurat/pbmc3k_tutorial.html
# # # # https://github.com/farrellja/URD/blob/master/Analyses/QuickStart/URD-QuickStart-AxialMesoderm.md


# Alternatively, we start from scratch and cluster only using BMP
# what I should have done from the beggining..

FACS_files = list.files(here("00_data_ingest","00_facs_raw_data","FACS"), full.names = TRUE)
#kidney data, from : https://hemberg-lab.github.io/scRNA.seq.course/tabula-muris.html

dat=read.delim(FACS_files[grep("Kidney",FACS_files)],sep=",",header=T)
rownames(dat)<-dat[,1]
dat<-dat[,-1]



cellIDs <- colnames(dat)
cell_info <- strsplit(cellIDs, "\\.")
Well <- lapply(cell_info, function(x){x[1]})
Well <- unlist(Well)
Plate <- unlist(lapply(cell_info, function(x){x[2]}))
Mouse <- unlist(lapply(cell_info, function(x){x[3]}))

#read annotations
ann=read.table("./00_data_ingest/00_facs_raw_data/annotations_FACS.csv",sep=",",header=T)
ann<-ann[match(cellIDs,ann[,1]),]
celltype<-ann[,3]

#cell type composition of tissue:
summary(factor(celltype))

#quality filters
cell_anns<-data.frame(mouse = Mouse,well=Well,type=celltype)
rownames(cell_anns)<-colnames(dat)
sce<-SingleCellExperiment(assays=list(counts = as.matrix(dat)),colData = cell_anns)

isSpike(sce, "ERCC") <- grepl("ERCC-", rownames(sce))
#remove genes that are not expressed in any cell
keep_feature<-rowSums(counts(sce)>0) >0
sce<-sce[keep_feature,]




#
