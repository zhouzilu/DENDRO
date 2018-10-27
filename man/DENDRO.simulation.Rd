\name{DENDRO.simulation}
\alias{DENDRO.simulation}
\title{
Simulation DENDRO performance
}
\description{
DENDRO.simulation assess the DENDRO performance given an imaginary clonal tree with cell numbers, read depth etc. See argument below.
}
\usage{
DENDRO.simulation(kprob = NULL, lprob = NULL, filt = 0, m = 100, n = 1000, epi = 0.001, RD = NULL, ref, k = NULL, subtype = 1, rpt = 100, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kprob}{
  Distribution of cells in each subclone. The order of probability in kprob see tutorial. Default NULL, indicates uniform distribution of cells.
}
  \item{lprob}{
  Distribution of mutations along each branch in the tree. The order of probability in lprob see tutorial. Default NULL, indicates uniform distribution of mutations.
}
  \item{filt}{
  filt threashold indicates mutations with less than filt mutations are filtered out.
}
  \item{m}{
  m indicates number of cells in your simulation. Default 100
}
  \item{n}{
  n indicates number of genes in your simulation. Default 1000
}
  \item{epi}{
  Sequencing error and rare RNA editing combined rate. Default 0.001 according to Illunima.
}
  \item{RD}{
  RD indicates read depth of the overall sequencing. Default uses data from Deng et al. 2014, which is 45X with 10,000,000 reads per cell.
}
  \item{ref}{
  The reference dataset, with X1 and X2 matrices, indicating two individual allele read counts.
}
  \item{k}{
  When kprob is NULL, this will decide the total number of subclones. Combined with subtype, they will determine the clonal tree structure. Default NULL.
}
  \item{subtype}{
  When k=4, there are two different clonal tree structure. This will determine which one is it. Default 1
}
  \item{rpt}{
  How many round of simulation do you want to run to generate the estimation and confidence interval. Default 100.
}
  \item{plot}{
  Whether you want to plot out the simulation example tree and results.
}
}
\value{
Summary statistic matrix with mean and confidence interval.
}
\references{
Deng, Qiaolin, et al. "Single-cell RNA-seq reveals dynamic, random monoallelic gene expression in mammalian cells." Science 343.6167 (2014): 193-196.
}
\author{
Zilu Zhou
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("DENDROplan_ref")
res=DENDRO.simulation(RD=4.5,n=100,ref=ref,k=4,subtype=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
