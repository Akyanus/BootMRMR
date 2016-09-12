\name{bmrmr.weight.cutoff}
\alias{bmrmr.weight.cutoff}
\title{Selection of informative geneset using gene weights obtained from the Bootstrap-MRMR technique}
\description{The function enables to find set of informative genes which are obtained based on weights computed from Bootstrap-MRMR technique.}
\usage{bmrmr.weight.cutoff(x, y, s, n)}
\arguments{
   \item{x}{x is a N by p dataframe of gene expression, where rows are genes and columns are as samples or subjects (gene names are taken as row names). Each cell or entry represents the expression level of a gene for a sample or subject.}   
  \item{y}{y is a p by 1 numeric vector having elements as 1 and -1 representing the sample labels of samples or subjects (for two class problems, i.e. stress or control respectively).}
  \item{s}{s is a numeric constant representing the number of bootstrap samples drawn (s must be sufficiently large)}
  \item{n}{n must be a numeric constant representing the number of informative genes to be selected from the large gene space.}
}
\value{The function returns a set of genes, which are highly informative to the trait or condition under consideration based on the computed weights form Bootstrap-MRMR technique.}

\author{Samarendra Das}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
s=80
n=20
bmrmr.weight.cutoff(x, y, s, n)
}
\keyword{gene}
\keyword{weights}
\keyword{informative geneset}
\keyword{Bootstrap}
\keyword{MRMR}
