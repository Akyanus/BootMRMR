\name{Weights.mrmr}
\alias{Weights.mrmr}
\title{Computation of MRMR weights for gene selection}
\description{The MRMR weights associated with each gene in the dataset are computed by using the MRMR technique for informative gene selection.}
\usage{Weights.mrmr(x, y)}
\arguments{
   \item{x}{x is a N by p dataframe of gene expression, where, rows as genes and columns as samples (with row names as gene names/ids)}
  \item{y}{y is a p by 1 numeric vector with entries 1 and -1 representing sample labels, where 1 and -1 represents the sample label of subjects/ samples for stress and control condition respectively.}
 }
\details{This function returns a vector of MRMR weights for all genes in the dataset.}
\author{Samarendra Das}
\references{Ding, C and Peng, H (2005). Minimum redundancy feature selection from microarray gene expression data. J. Bioinformatics Comput Biol 3(2):185-205.}
\examples{

data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
Weights.mrmr(x, y)
}
\keyword{gene}
\keyword{weights}
\keyword{MRMR}