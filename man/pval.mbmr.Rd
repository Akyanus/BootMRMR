\name{pval.mbmr}
\alias{pval.mbmr}
\title{Computation of statistical significance values for genes using Modified Bootstrap MRMR technique for a particular trait/condition}
\description{The statisical significance values (p-values) will be computed for all the genes in the dataset from the non-parametric test "H0: i-th gene is not informative against H1: i-th gene is informative" for selection of informative genes using Modified Bootstrap MRMR technique.}
\usage{pval.mbmr(x, y, m, s, Q, plot)}
\arguments{
  \item{x}{x is a N by p data frame of gene expression values where rows represent genes and columns represent samples/subject/time point. Each cell entry represents the expression level of a gene in a sample/subject (row names of x as gene names/gene ids).}
  \item{y}{y is a p by 1 numeric vector with entries 1 and -1 representing sample labels, where 1 and -1 represents the sample label of subjects/ samples for stress and control condition respectively.}
  \item{m}{m is a scalar representing the size of the Modified Bootstrap Sample (i.e. Out of p samples/subjects, m samples/subjects are randomly drawn with replacement, which constitutes one Modified Bootstrap Sample).}
   \item{s}{s is a scalar representing the number of Modified Bootstrap samples (i.e. number of times each of the m samples/subjects will be resampled from p samples/subjects). }
  \item{Q}{Q is a scalar representing the quartile value of the gene rankscores (lies within 1/N to 1), usually the second quartile, i.e. 0.5 or third quartile i.e. 0.75.}
  \item{plot}{plot is a character string must either take logical value TRUE/FALSE representing whether to plot the statistical significance values of genes in the dataset.}
}
\value{The function returns a vector of p-values for all the genes from the given statistical test in the gene space/dataset using Modified Bootstrap MRMR technique.}

\author{Samarendra Das}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
m=36
s=80
Q=0.5
pval.mbmr(x, y, m, s, Q, plot=FALSE)
}
\keyword{gene}
\keyword{rankscore}
\keyword{bootstrap}
\keyword{p-value}