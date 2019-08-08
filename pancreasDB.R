

library(parallel)
library(scran)
library(pheatmap)
library(Seurat)
library(Rmagic)

intestine_file = "/Users/alejandrog/Downloads/SRA653146_SRS2874280.sparse-RPKM.RData"

panglao.path = "/Users/alejandrog/Downloads/"
pancreas.files = c("SRA784304_SRS3823095.sparse.RData", "SRA784304_SRS3823097.sparse.RData"  ,
                    "SRA784304_SRS3823098.sparse.RData" , "SRA784304_SRS3823096.sparse.RData","SRA784304_SRS3823099.sparse")
pancreas.files = paste(panglao.path, pancreas.files, sep = "")

conditions = c("E15_progenitors_Krentz","E18_progenitors_Krentz" ,"E18_green_Krentz","E15_yellowGreen_Krentz","E18_yellow_Krentz")

notch.genes = c("Dll1",   "Dll3"   ,"Dll4",   "Dtx1",   "Jag1"  , "Jag2", "Notch1", "Notch2", "Notch3", "Notch4", "Mfng",  "Rfng"  , "Lfng",   "Dlk1",   "Dlk2")
bmp.receptors<-c( "Bmpr1a" ,"Bmpr1b" ,"Acvr1"  ,"Acvrl1" ,"Acvr1b" ,"Tgfbr1" ,"Acvr1c" ,"Acvr2a", "Acvr2b", "Bmpr2" ,"Tgfbr2")


load.PanglaoDB<-function(input_file = "", min_colSums = 3, max_colSums = 6, max_rowSums = 6,
                            k = 10, data_type = "10x", plot.all = T){

  load(input_file)

  if(data_type =="10x"){
    #sm2 refers to RPKM values reported in PanglaoDB from SmartSeq data
    # 10x is not normalized in the same way and so, we only get sm
    sm2 = sm
    rm(sm)
  }

  #minimum percentage of cells in which a gene should be expresed
  min.percent = 0.05
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
#
normaliseAndImpute<-function(sm2_filtered = c(), k_magic = 10,plot.all = F){

  # SCRAN normalization
  # The counts slot can be filled with RPKM values, since it is used for dropout correction
  sce<-SingleCellExperiment(list(counts = sm2_filtered))
  clusters <- quickCluster(sce, min.size=100) # TAKES a long time
  sce <- computeSumFactors(sce, cluster=clusters)
  # Now we can normalize the data
  sce<-normalize(sce)

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
  sce_magic  = magic(t(as.matrix(sce.seurat[['RNA']]@data )),k=k_magic,n.jobs =-1)
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

plotPathwayMagic<-function(sce_magic_data= c() , gene.list = bmp.receptors, cut_k = 10){

    gene.names = row.names(sce_magic_data)
    p2 =pheatmap(sce_magic_data[which(gene.names %in% gene.list),] ,cluster_rows = F,cutree_cols = cut_k,show_colnames = F)

    return(p2[[4]])

}

#we can directly download the data from paper's website to compare

norm_matrices_data = grep("normalized",list.files("Downloads/Krentz_2019/"),value = T)

loadFromPaper<-function(file_name = "Downloads/Krentz_2019/GSM3402516_E18_green_normalized_counts_matrix.txt.gz",res=list()){

  aa<-read.table(file_name,row.names = 1,header = T)
  aa<-t(aa)
  aa_data = aa[5:dim(aa)[1],]
  aa_data<- matrix(mapply(aa_data, FUN=as.numeric),  ncol=ncell, nrow=ngene)
  row.names(aa_data)<-row.names(aa)[-c(1:4)]
  colnames(aa_data) = colnames(aa)

  colnames_simple = str_split(str_split(colnames(aa_data),"_",simplify=T )[,3], "-", simplify=T)[,1]
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
