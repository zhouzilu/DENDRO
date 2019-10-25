DENDRO.dist = function(X,N,Z,epi=0.01,show.progress=TRUE){
  Pg = rowSums(Z,na.rm=T)/ncol(Z)
  dist = as.dist(SNV.dist(N,X,Pg,epi,show.progress))
  dist = dist-min(dist)+1
  return(dist)
}

DENDRO.dist.v1 = function(X,N,Z,epi=0.01,show.progress=TRUE){
  Ng=rowSums(!is.na(Z))
  #Pg = cbind((rowSums(Z==0,na.rm=T)+1)/Ng,(rowSums(Z==1,na.rm=T)+1)/Ng,(rowSums(Z==2,na.rm=T)+1)/Ng)
  #Pg = Pg/rowSums(Pg)
  Pg = cbind((rowSums(N-X,na.rm=T)/rowSums(N,na.rm=T))^2,2*rowSums(X,na.rm=T)*rowSums(N-X,na.rm=T)/(rowSums(N,na.rm=T)^2),(rowSums(X,na.rm=T)/rowSums(N,na.rm=T))^2)
  dist = as.dist(SNV.dist.v1(N,X,Pg,epi,show.progress))
  dist = dist-min(dist)+1
  return(dist)
}


DENDRO.cluster = function(dist,method='ward.D',plot=TRUE,label=NULL){
  clust=hclust(dist,method=method)
  if(plot){
    dend=as.dendrogram(clust)
    if(!is.null(label)){
      labels_colors(dend) =
        colorspace::rainbow_hcl(
          length(unique(label))
          )[label][order.dendrogram(dend)]
    }
    plot(dend,main='DENDRO Result')
    if(!is.null(label)){
      cols <- colorspace::rainbow_hcl(length(unique(label)))
      legend("topright", legend = 1:length(unique(label)),
             fill = cols, border = cols, bty = "n")
    }
  }
  return(clust)
}

DENDRO.tree = function(Z_cluster,label_cluster=NULL){
  ret <- phyclust.edist(t(Z_cluster), edist.model = .edist.model[4])
  if(is.null(label_cluster)){
    label_cluster=colorspace::rainbow_hcl(ncol(Z_cluster))
  }else{
    label_cluster=
      colorspace::rainbow_hcl(length(unique(label_cluster)))[label_cluster]
  }
  # summary(ret)
  # ret=log(ret)-min(log(ret))
  (ret.tree <- nj(ret))
  ret.tree$tip.label=colnames(Z_cluster)
  # plotnj(ret.tree,c(1,1,2,2,2,3,3),  show.tip.label = TRUE)

  plotnj(ret.tree,tip.color=label_cluster,show.tip.label = TRUE)

}
