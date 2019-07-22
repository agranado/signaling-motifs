

load()

load("Downloads/SRA653146_SRS2874280.sparse-RPKM.RData")


sm2_filtered=sm2[,Matrix::colSums(sm2>0)>1000]

bmp.receptors<-c( "Bmpr1a" ,"Bmpr1b" ,"Acvr1"  ,"Acvrl1" ,"Acvr1b" ,"Tgfbr1" ,"Acvr1c" ,"Acvr2a", "Acvr2b", "Bmpr2" ,"Tgfbr2")

gene.names = do.call(rbind,strsplit(row.names(sm2_filtered),"_"))[,1]
gene.names[which(gene.names %in% bmp.receptors)]
