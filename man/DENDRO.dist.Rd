\name{DENDRO.dist}
\alias{DENDRO.dist}
\title{
DENDRO specific genetic divergence evaluation function
}
\description{
Calculate the cell-to-cell divergence matrix given variant allele read counts, total allele read counts, estimated mutations and sequencing error rate. This method accounts for transcriptional bursting and sequencing error with a Beta-Binomial framework. This function is linear with number of cells and number of mutations.
}
\usage{
DENDRO.dist(X, N, Z, epi = 0.01, show.progress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  An matrix contains variants allele read counts across cell (column) and loci (row).
}
  \item{N}{
  An matrix contains total allele read counts across cell (column) and loci (row).
}
  \item{Z}{
  An mutation indicator matrix (1 for mutation, 0 for normal) across cell (column) and loci (row).
}
  \item{epi}{
  Sequencing error and rare RNA editing combined rate. Default 0.01 according to Illunima.
}
  \item{show.progress}{
  Whether to show the divergence calculation programs. For large and diverse cell population, this function will take some time and we recommend tracking progress.
}
}
\value{
`dist` returns an object of class `"dist"`.
See https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/dist for more details
}
\references{
https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/dist
}
\author{
Zilu Zhou
}
\examples{
data("DENDRO_demo")
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
