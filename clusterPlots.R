#plots for lab meeting



pathwayClustersToGlobal(seurat.pathway = seurat.obj, which.pathway = "bmp",tiss  = tiss){
    all.cells =tiss[[]]$cell
    signaling.clusters = as.numeric(levels(seurat.pathway[[]]$seurat_clusters))


    for(i in 1:length(signaling.clusters)){
        seurat.pathway[[]] %>% dplyr::filter(seurat_clusters==signaling.clusters[i]) %>% select(cell) ->cells.clus.2
        meta.clus.2 = data.frame(cell= cells.clus.2, value = rep(1,length(cells.clus.2)))

        bmp.2.global = rep(0, length(all.cells))
        bmp.2.global[which(all.cells %in% meta.clus.2$cell)] =1
        tiss = AddMetaData(tiss,bmp.2.global,col.name = paste(which.pathway,"." , toString(signaling.clusters[i]),".global",sep=""))
    }

    signaling.cluster.names = paste(which.pathway,".", as.character(0:(length(signaling.clusters)-1))  ,".global",sep="")
    return(tiss)
}

pathwayClusterID <-function(seurat.pathway= c(), which.pathway = "bmp"){

    signaling.clusters = as.numeric(levels(seurat.pathway[[]]$seurat_clusters))
    signaling.cluster.names = paste(which.pathway,".", as.character(0:(length(signaling.clusters)-1))  ,".global",sep="")

}
