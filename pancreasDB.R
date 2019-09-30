

library(parallel)
library(scran)
library(pheatmap)
library(Seurat)
library(Rmagic)
library(RColorBrewer)
library(sva)

intestine_file = "/Users/alejandrog/Downloads/SRA653146_SRS2874280.sparse-RPKM.RData"

#not in Panglao
pancreas2.files = c("GSE77980_mouse_islets_rpkm_Xin2016.txt.gz")
conditions2 = c("Pancreas_islets_Xin")

### FROM panglao
panglao.path = "/Users/alejandrog/Downloads/"
pancreas.files = c("SRA784304_SRS3823095.sparse.RData", "SRA784304_SRS3823097.sparse.RData"  ,
                    "SRA784304_SRS3823098.sparse.RData" , "SRA784304_SRS3823096.sparse.RData","SRA784304_SRS3823099.sparse")
pancreas.files = paste(panglao.path, pancreas.files, sep = "")


conditions = c("E15_progenitors_Krentz","E18_progenitors_Krentz" ,"E18_green_Krentz","E15_yellowGreen_Krentz","E18_yellow_Krentz")
#panglao has some errors in the annotation!
corrected_conditions = c("E15_red",  "E18_green", "E18_red","E15_yellow_green")
# NOTE:  After Aug 19th. The list was saved and re-named to match the original annotation order



notch.genes = c("Dll1",   "Dll3"   ,"Dll4",   "Dtx1",   "Jag1"  , "Jag2", "Notch1", "Notch2", "Notch3", "Notch4", "Mfng",  "Rfng"  , "Lfng",   "Dlk1",   "Dlk2")
bmp.receptors<-c( "Bmpr1a" ,"Bmpr1b" ,"Acvr1"  ,"Acvrl1" ,"Acvr1b" ,"Tgfbr1" ,"Acvr1c" ,"Acvr2a", "Acvr2b", "Bmpr2" ,"Tgfbr2")


# # # # # # # # # # # # # #
# Main Functions:
#

load.PanglaoDB<-function(input_file = "", min_colSums = 3, max_colSums = 6, min.percent = 0.05, max_rowSums = 6,
                            k = 10, data_type = "10x", plot.all = T, user_matrix = F,user_defined_matrix = c()){


  #Normally we want to use this function to load UMI / RPKM data from PanglaoDB
  #However, the function also performs basic filtering, so we can use it with other source data
  if(!user_matrix){
    #If this is a PanglaoDB matrix
    #LOAD the file
      load(input_file)

      if(data_type =="10x"){
        #sm2 refers to RPKM values reported in PanglaoDB from SmartSeq data
        # 10x is not normalized in the same way and so, we only get sm
        sm2 = sm
        rm(sm)
      }
  }else{
    #User provided an input matrix:
    sm2 = user_defined_matrix

  }
  #minimum percentage of cells in which a gene should be expresed

  #Manually Filtering cells
  if(plot.all){
    x11()
    par(mfrow = c(2,2))
    #make histograms to establish threshold
    hist(log10(Matrix::colSums(sm2>0)),main = "N expressed genes per cell")
    hist(log10(Matrix::colSums(sm2)), main = "Reads per cell")
  }
  #min 10^3 genes
  #max 10^6 UMIs (10x data)
  sm2_filtered=sm2[,log10(Matrix::colSums(sm2>0)) > min_colSums & log10(Matrix::colSums(sm2)) < max_colSums  ]
  #filter genes that are expressed in less than 5% of cells in the dataset
  sm2_filtered = sm2_filtered[Matrix::rowSums(sm2_filtered>0) >round(dim(sm2_filtered)[2] * min.percent)  & log10(Matrix::rowSums(sm2_filtered)) <max_rowSums,]

  if(plot.all){
    hist(log10(Matrix::colSums(sm2_filtered>0)),main = "N expressed genes per cell (filtered)")
    hist(log10(Matrix::colSums(sm2_filtered)), main = "Reads per cell (filtered)")
  }
  #extract names from Panglao format:
  gene.names = do.call(rbind,strsplit(row.names(sm2_filtered),"_"))[,1]
  ensembl.names = do.call(rbind,strsplit(row.names(sm2_filtered),"_"))[,2] #second part of the string contains this name
  ensembl.names<-do.call(rbind,strsplit(ensembl.names,"\\."))[,1]

  row.names(sm2_filtered)<-gene.names



  return(sm2_filtered)
}
#CPU time, Mac 8min single core for pancreas.files[1]

normalizeSCRAN<-function(sm2_filtered){
  # SCRAN normalization
  # We need a SCE object so we will convert our matrix to that,
  # For MAGIC, we need a UMI matrix so we will convert to Seurat so
  # we can extract the count matrix from there (not elegant)

  # The counts slot can be filled with RPKM values, since it is used for dropout correction
  sce<-SingleCellExperiment(list(counts = sm2_filtered))
  clusters <- quickCluster(sce, min.size=100) # TAKES a long time
  sce <- computeSumFactors(sce, cluster=clusters)
  # Now we can normalize the data
  # Compute normalised expression values from count data in a SingleCellExperiment object, using the size factors stored in the object.
  sce<-scater::normalize(sce)

  #let's convert to Seurat:
  sce.seurat <- as.Seurat(sce)

  norm_count_matrix = as.matrix(sce.seurat[['RNA']]@data )

  return(norm_count_matrix)

}

normaliseAndImpute<-function(sm2_filtered = c(), k_magic = 10,plot.all = F){

  # SCRAN normalization
  # The counts slot can be filled with RPKM values, since it is used for dropout correction
  sce<-SingleCellExperiment(list(counts = sm2_filtered))
  clusters <- quickCluster(sce, min.size=100) # TAKES a long time
  sce <- computeSumFactors(sce, cluster=clusters)
  # Now we can normalize the data
  # Compute normalised expression values from count data in a SingleCellExperiment object, using the size factors stored in the object.
  sce<-scater::normalize(sce)

  #let's convert to Seurat:
  sce.seurat <- as.Seurat(sce)

  #Seurat converts _ to - in gene names. Let's split it again
  #gene.names<-do.call(rbind,strsplit(gene.names,"-"))[,1]
  gene.names<-row.names(sm2_filtered)


  # PLOT the bmp receptors after normalization but BEFORE imputation
  if(plot.all){
    subset_full =sm2_filtered[which(gene.names %in% bmp.receptors) ,]
    p1 = pheatmap(log2(1+subset_full),cluster_rows = F,show_colnames = F)
  }
  #MAGIC imputation
  sce_magic  = Rmagic::magic(t(as.matrix(sce.seurat[['RNA']]@data )),k=k_magic,n.jobs =-1)
  sce_magic_data = t(sce_magic$result)
  if(plot.all){
    p2 =pheatmap(sce_magic_data[which(gene.names %in% bmp.receptors),] ,cluster_rows = F,cutree_cols = 10,show_colnames = F)
    g <- grid.arrange(arrangeGrob(grobs= list(p1[[4]],p2[[4]]),ncol=2))
    x11()
    plot(g)
    return(list(p1,p2,sce_magic_data,g))
  }else{
    return(sce_magic_data)
  }

}

plotPathwayMagic<-function(sce_magic_data= c() , gene.list = bmp.receptors,
                      cut_k = 10, main_title = "",ann_df = c()){

    #22 distinct colors
    distinct_colors = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
                        '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
                        '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080', '#ffffff', '#000000')

    gene.names = row.names(sce_magic_data)

    #with annotation for colors and cell types
    # List of color must match a given cell type
    if(length(ann_df)>0){

      mycolors <-  distinct_colors[1:length(unique(ann_df$cell_type))]
      names(mycolors) <- unique(ann_df$cell_type)
      mycolors <- list(cell_type = mycolors) # Here use the name of the variable you want to annotate!

      # subset genes from the list of pathway-genes
      # Sort then alphabetically
      # Create a heatmap with 10 cuts and save it (will be returned by the function)
      p2 =pheatmap(sce_magic_data[sort(gene.names[which(gene.names %in% gene.list)]),] ,cluster_rows = F,cutree_cols = cut_k,
            show_colnames = F,main = main_title,annotation_col = ann_df,
              annotation_colors = mycolors)

    #No annotations
    }else{
      p2 =pheatmap(sce_magic_data[sort(gene.names[which(gene.names %in% gene.list)]),] ,cluster_rows = F,cutree_cols = cut_k,
            show_colnames = F,main = main_title)
    }

    return(p2[[4]])
}



# # # # # # # # # # # # # # # # #
# Krentz 2018 Stem Cell Reports
# Pancreas Progenitors 10x sc RNA sequ
# E15.5 18.5
# Data from Panglao and from NCBI (All datasets were cross-validated for annotation)
# Metadata: From NCBI

#we can directly download the data from paper's website to compare
# There are some errors in the annotation of PanglaoDB for this particular dataset
# WE will load the raw data from the paper and compare with Panlao to correct for errors
# There is also a table of number of cells per sample in the paper, we want to use that to compare to
krentz_paper_files=grep("normalized",list.files("/Users/alejandrog/Downloads/Krentz_2019/"),value = T)
krentz_paper_files = paste("Downloads/",krentz_paper_files,sep ="")
#Returns a list of lists: [[1]] is a list of magic_matrices [[2]] is a list of the associated metadata
loadFromPaper<-function(krentz_paper_files = krentz_paper_files,res=list()){

  for(i in 1:length(krentz_paper_files)){
      file_name = krentz_paper_files[i]
      aa<-read.table(file_name,row.names = 1,header = T)
      aa<-t(aa)
      #this has metadata
      # first 4 columns are the meta data
      aa_data = aa[5:dim(aa)[1],]
      ncell = ncol(aa_data)
      ngene = nrow(aa_data)
      # Data comes as string so we need to convert it to numeric
      aa_data<- matrix(mapply(aa_data, FUN=as.numeric),  ncol=ncell, nrow=ngene)
      row.names(aa_data)<-row.names(aa)[-c(1:4)]
      colnames(aa_data) = colnames(aa)

      # In this file, the names of cells have two underscores, and therefore the index we want is 4 instead of 3
      if(length(grep("E15_yellow_green",file_name))>0) barcode_index = 4 else barcode_index = 3

      colnames_simple = str_split(str_split(colnames(aa_data),"_",simplify=T )[,barcode_index], "-", simplify=T)[,1]
      colnames(aa_data)<-colnames_simple

      #let's compile the meta data
      meta_krentz<-aa[1:4,]
      meta_krentz<-as.data.frame(meta_krentz)



      krentz_magic[[i]] = aa_data
      krentz_meta[[i]]  = meta_data
    }

  return(aa_data)
}

#NOTE E15_red has not cell_type annotation

loadMetaKrentz<-function(file_name = krentz_paper_files[1]){
  #This function only applies for that specific file format

  aa<-read.table(file_name,row.names = 1,header = T)
  aa<-t(aa)

  if(length(grep("E15_yellow_green",file_name))>0) barcode_index = 5 else barcode_index = 4

  sample_name = paste(str_split(str_split(file_name,"normalized",simplify = T)[,1],"_",simplify=T)[,3:barcode_index],collapse = "_")

  meta_krentz<-aa[1:4,]
  meta_krentz<-as.data.frame(meta_krentz)

  meta_krentz_E18_green<-t(meta_krentz)
  meta_krentz_E18_green<-as.data.frame(meta_krentz_E18_green)
  meta_krentz_E18_green$simple_name<-str_split(str_split(colnames(meta_krentz),"_",simplify = T)[,barcode_index-1],"-",simplify = T)[,1]

  return(list(meta_krentz_E18_green))

}


# # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # #
# Byrnes 2018 Nat Comm
# Whole Pancreas 10x sc RNA sequ
# E12.5, 14.5 17.5
# Data from Panglao  17.5 (UMIs) and from FigShare 14.5 (official, as Seurat object)
# From FigShare, we have metadata

# This function will load and process the Seurat objects available in FigShare
# Seruat objects for this paper are very old versions and are NOT compatible with Seurat 3
loadByrnesFigShare<-function(file_name="Downloads/Byrnes2018/E14_allcells_seur_ob.Rdata"){

  load(file_name)
  # EXtract raw and meta data
  byrnes_E14_meta<-seurat_ob@data.info
  byrnes_E14_raw = seurat_ob@raw.data
  # My function for SCRAN  + MAGIC
  byrnes_E14_magic <- normaliseAndImpute(byrnes_E14_raw)


  # Annotation of byrnes data object:
  byrnes_E14_ann = data.frame(cell_type = byrnes_E14_meta$assigned_ident[ which(row.names(byrnes_E14_meta) %in% colnames(byrnes_E14_raw)) ]  )


  # Returns the imputed data matrix for this object
  return(list(byrnes_E14_magic,byrnes_E14_ann))

  # From here: execture the plotting function
  # with annotation
  # byrnes_E14_p<-plotPathwayMagic(byrnes_E14_magic,main_title = "Byrnes E14 (2batches)",ann_df = byrnes_E14_ann)

}


# Xin dataset:

load.additionalSamples <- function(sampleID){
  if(sampleID =="Xin"){
    Xin = read.table("Downloads/pancreas2/GSE77980_mouse_islets_rpkm_Xin2016.txt.gz",header = T)
  }

}
#


# Intestine
# Load Aviv Regev intestine dataset

# We can make a function for loading / filtering ANY UMI matrix
  #Manually Filtering cells
filterData<-function(sm2= c(),plot.all = T,min_colSums = 3, max_colSums = 6, min.percent = 0.05, max_rowSums = 6){

  if(plot.all){
    x11()
    par(mfrow = c(2,2))
    #make histograms to establish threshold
    hist(log10(Matrix::colSums(sm2>0)),main = "N expressed genes per cell")
    hist(log10(Matrix::colSums(sm2)), main = "Reads per cell")
  }
  #min 10^3 genes
  #max 10^6 UMIs (10x data)
  sm2_filtered=sm2[,log10(Matrix::colSums(sm2>0)) > min_colSums & log10(Matrix::colSums(sm2)) < max_colSums  ]
  #filter genes that are expressed in less than 5% of cells in the dataset
  sm2_filtered = sm2_filtered[Matrix::rowSums(sm2_filtered>0) >round(dim(sm2_filtered)[2] * min.percent)  & log10(Matrix::rowSums(sm2_filtered)) <max_rowSums,]

  if(plot.all){
    hist(log10(Matrix::colSums(sm2_filtered>0)),main = "N expressed genes per cell (filtered)")
    hist(log10(Matrix::colSums(sm2_filtered)), main = "Reads per cell (filtered)")
  }

  return(sm2_filtered)
}








old.code<-function(){
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
}
