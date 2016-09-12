\name{pval.bmrmr}
\alias{pval.bmrmr}
\title{Compuation of statistical significance values for genes using Bootstrap-MRMR technique}
\description{The function computes ths statisical significance values for the genes from the non-parametric test "H0: i-th gene is not informative against H1: i-th gene is informative" for selection of informative genes using Bootstrap-MRMR technique}
\usage{pval.bmrmr(x, y, s, Q, plot)}
\arguments{
  \item{x}{x is a N by p data frame of gene expression values where rows represent genes and columns represent samples/subject/time point. Each cell entry represents the expression level of a gene in a sample/subject (row names of x as gene names/gene ids).}
  \item{y}{y is a p by 1 numeric vector with entries 1 and -1 representing sample labels, where 1 and -1 represents the sample label of subjects/ samples for stress and control condition respectively.}
  \item{s}{s is a scalar representing the number of bootstraps generated, s must be sufficiently large (i.e. number of times  bootstrap samples are generated}
  \item{Q}{Q is a scalar representing the quartile value of the gene rankscores (lies within 1/N to 1), usually the second quartile (Q2), i.e. 0.5 or third quartile (Q3) i.e. 0.75 is taken.}
  \item{plot}{plot is a character string must either take logical value TRUE/FALSE representing whether to plot the statistical significance values of genes in the dataset.}
}
\value{The function returns a vector of p-values for all the genes from the given statistical test in the dataset using Bootstrap-MRMR technique.}

\author{Samarendra Das}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
s=80
Q=0.5
pval.bmrmr(x, y, s, Q, plot=FALSE)
}
\keyword{gene}
\keyword{rankscore}
\keyword{bootstrap}
\keyword{p-value}
\keyword{non-parametric test}