\name{DENDRO.recalculate}
\alias{DENDRO.recalculate}
\title{
Recalculate mutation profile for each cluster
}
\description{
This function calculates the mutation profile for each cluster, after pooling cells within each cluster together. Such that, the cells within each cluster have same mutation profiles and the result is more robust. The estimation is based on a maximum likelihood approach. Loci with mutation observed in all or no subclones are removed.
}
\usage{
DENDRO.recalculate(X, N, Info, DENDRO_label, cluster.name = NULL, top = NULL, epi = 0.001, m = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  An matrix contains variants allele read counts across cell (column) and loci (row).
}
  \item{N}{
  An matrix contains total allele read counts across cell (column) and loci (row).}
\item{Info}{
  Mutation loci information. Matrix with row number same to number of loci and column number > 1. This is optional. Default NULL
}
  \item{DENDRO_label}{
  An integer 1D vector that decide the labeled cluster in each cell. Estimated by DENDRO.dist and DENDRO.cluster.
}
  \item{cluster.name}{
  You can label the name of the cluster after observing the cell composition. This is optional. Default NULL
}
  \item{top}{
  Since many mutations may be observed, `top` ask how many top mutation sites you want to selected for downstream analysis based on its marginal likelihood. This is optional. Default NULL
}
  \item{epi}{
  Sequencing error and rare RNA editing combined rate. Default 0.001 according to Illunima.
}
  \item{m}{
  The polidy in maximum likelihood mutation estimation. Default 2.
}
}
\value{
\item{X}{An matrix contains variants allele read counts across subclones (column) and loci (row).}
\item{N}{An matrix contains total allele read counts across subclones (column) and loci (row).}
\item{Z}{An mutation indicator matrix (1 for mutation, 0 for normal) across subclones (column) and loci (row).}
\item{Info}{Input Info after filtering}
\item{Z_cluster_lg}{An likelihood matrix of mutation inference across subclones (column) and loci (row).}
}
\references{
Li, B., et al., A likelihood-based framework for variant calling and de novo mutation detection in families. PLoS Genet, 2012. 8(10): p. e1002944.
}
\author{
Zilu Zhou
}
\examples{
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
cluster=DENDRO.cluster(dist,label=demo_qc$label)
icd = DENDRO.icd(dist,cluster)
DENDRO_label = cutree(cluster,3)
demo_cluster = DENDRO.recalculate(demo_qc$X,demo_qc$N, demo_qc$Info, demo_qc$DENDRO_label, cluster.name=c('Cluster3','Cluster2','Cluster1'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
