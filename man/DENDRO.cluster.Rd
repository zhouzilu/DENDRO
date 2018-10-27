\name{DENDRO.cluster}
\alias{DENDRO.cluster}
\title{
DENDRO specific cluster method on genetic divergence matrix
}
\description{
Cluster cells based on the genetic divergence matrix using hierachical cluster methods. 
}
\usage{
DENDRO.cluster(dist, method = "ward.D", plot = TRUE, label = NULL)
}
\arguments{
  \item{dist}{
  A dissimilarity structure as produced by `dist`
}
  \item{method}{
  Submethod for `hclust`. Default is "ward.D", because statistical integreity
}
  \item{plot}{
  TRUE or FALSE, decide whether to plot the result. Default value is TRUE.
}
  \item{label}{
  An integer 1D vector that decide the label color. Default NULL, where no color label. 
}
}
\value{
An object of class hclust which describes the tree produced by the clustering process. Check https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust for more detail.

If `plot==TRUE`, a hclust plot will also be displayed.
}
\references{
https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust
}
\author{
Zilu Zhou
%%  ~~who you are~~
}
\examples{
data("DENDRO_demo")
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
cluster=DENDRO.cluster(dist,label=demo_qc$label)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
