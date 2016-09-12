\name{weight.mbmr}
\alias{weight.mbmr}
\title{Computation of weights for informative gene selection using Modified Bootstrap MRMR technique}
\description{Weights associated with genes in a dataset computed from the Modified Bootstrap MRMR technique will provide a reliable measure for informative gene selection.}
\usage{weight.mbmr(x, y, m, s, plot)}
\arguments{
  \item{x}{x is a N by p dataframe of gene expression, where rows are genes and columns are as samples/subjects (gene names are taken as row names).}   
  \item{y}{y is a p by 1 numeric vector with entries 1 and -1 representing sample labels, where 1 and -1 represents the sample label of subjects/ samples for stress and control condition respectively.}
  \item{m}{m is a scalar representing the size of the Modified Bootstrap Sample (i.e. Out of p samples/subjects, m samples/subjects are randomly drawn with replacement, which constitutes one Modified Bootstrap Sample).}
 \item{s}{s is a scalar representing the number of Modified Bootstrap samples (i.e. number of times each of the m samples/subjects will be resampled from p samples/subjects). }
 \item{plot}{plot is a character string must either take logical value TRUE/FALSE representing whether to plot the weights of genes in the dataset.}
 }
\details{The function returns a vector of weights associated with each genes in the dataset using Modified Bootstrap MRMR technique.}
\author{Samarendra Das}
\references{Wang J, Chen L, Wang Y, Zhang J, Liang Y, Xu D (2013) A Computational systems biology study for understanding salt tolerance mechanism in Rice. PLoS one 8(6): e64929.}
\examples{

data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
m=36
s=80
weight.mbmr(x, y, m, s, plot=FALSE)
}
\keyword{gene}
\keyword{weights}
\keyword{modified bootstrap}