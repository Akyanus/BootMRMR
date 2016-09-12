\name{bootmr.weight}
\alias{bootmr.weight}
\title{Computation of weights for informative genes/ geneset selection using Bootstrap-MRMR technique}
\description{The function computes the weights associated with each genes for a given dataset using Bootstrap-MRMR technique.}
\usage{bootmr.weight(x, y, s, plot)}
\arguments{
  \item{x}{x is a N by p dataframe of gene expression, where rows are genes and columns are as samples/subjects (gene names are taken as row names). Each cell/entry represents the expression level of a gene in a sample/subject.}   
  \item{y}{y is a p by 1 numeric vector having elements as 1/-1 representing the sample labels of samples/subjects (for two class problems, i.e. stress/control)}
  \item{s}{s is a numeric constant representing the number of bootstrap samples drawn (s must be sufficiently large)}
  \item{plot}{plot is a character string must either take logical value TRUE/FALSE representing whether the plot of the gene weights of all genes in the dataset needs to be constructed or not.}
 }
\details{The function returns a vector of weights associated with each genes computed from Bootstrap-MRMR technique for a given dataset.}
\author{Samarendra Das}

\examples{

data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
s=80
bootmr.weight(x, y, s, plot=FALSE)
}
\keyword{gene}
\keyword{weights}
\keyword{bootstrap}