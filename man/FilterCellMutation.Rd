\name{FilterCellMutation}
\alias{FilterCellMutation}
\title{
Filter out low expressed gene and high dropout cells based on read counts
}
\description{
This function filter out low expressed gene (cut.off.VAF) and high dropout cells (cut.off.sd) based on read counts and plot out the summary distributions.
}
\usage{
FilterCellMutation(X, N, Z, Info = NULL, label = NULL, cut.off.VAF = 0.05, cut.off.sd = 5, plot = TRUE)
}
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
  \item{Info}{
  Mutation loci information. Matrix with row number same to number of loci and column number > 1. This is optional. Default NULL
}
  \item{label}{
  An integer 1D vector that decide the label color by type. The label can be assigned due to piror information such as individuals or site. Default NULL, where no color label. 
}
  \item{cut.off.VAF}{
  Variant allele frequency (VAF, sum(Z)/length(Z) ) filter threshold. Variants with variants allele frequency less than cut.off.VAF (too rare) or greater than 1-cut.off.VAF (too common) are filtered out. Default 0.05
}
  \item{cut.off.sd}{
  Cell total expression (sum(N)) filter threshold. Cells with total expression less than mean(N)-cut.off.sd*sd(N) (too small) or greater than mean(N)+cut.off.sd*sd(N) (double or triplet) are filtered out. Default 5
}
  \item{plot}{
  TRUE or FALSE, decide whether to plot the distribution and cut off (red line). Default value is TRUE.
}
}
\value{
List of filtered matrix
\item{X}{Input X after filtering}
\item{N}{Input N after filtering}
\item{Z}{Input Z after filtering}
\item{Info}{Input Info after filtering}
\item{label}{Input label after filtering}
}
\author{
Zilu Zhou
}
\examples{
data("DENDRO_demo")
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
