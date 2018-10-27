\name{DENDRO.tree}
\alias{DENDRO.tree}

\title{
Construct a phylogenetic tree
}
\description{
This method constructs a phylogenetic tree with Neighbor Joining method (implemented in `phangorn`), suggesting the evolutionary relationship between different subclones. 
}
\usage{
DENDRO.tree(Z_cluster, label_cluster = NULL)
}
\arguments{
  \item{Z_cluster}{
  An mutation indicator matrix (1 for mutation, 0 for normal) across subclones (column) and loci (row).
}
  \item{label_cluster}{
  Cluster labels, which will be used to color phylogenetic tree node. Better to track cells and clusters.  This is optional. Default NULL
}
}
\value{
 void. 
 A phylogenetic tree will be plotted.
}
\references{
https://cran.r-project.org/web/packages/phangorn/index.html
}
\author{
Zilu Zhou
}
\examples{
data("DENDRO_demo")
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
cluster=DENDRO.cluster(dist,label=demo_qc$label)
icd = DENDRO.icd(dist,cluster)
DENDRO_label = cutree(cluster,3)
demo_cluster = DENDRO.recalculate(demo_qc$X,demo_qc$N, demo_qc$Info, demo_qc$DENDRO_label, cluster.name=c('Cluster3','Cluster2','Cluster1'))
DENDRO.tree(demo_cluster$Z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
