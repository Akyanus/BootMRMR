\name{mrmr.cutoff}
\alias{mrmr.cutoff}
\title{Informative geneset selection using MRMR weights}
\description{The function returns the informative genes/ geneset for the particular trait/condition under investigation using Maximum Relevance and Minimum Redundancy (MRMR) technique.}
\usage{mrmr.cutoff(x, y, n)}
\arguments{
  \item{x}{x is a N by p data frame of gene expression values where rows represent genes and columns represent samples/subject/time point. Each cell entry represents the expression level of a gene in a sample/subject (row names of x as gene names/gene ids).}
  \item{y}{y is a p by 1 numeric vector with entries 1 and -1 representing sample labels, where 1 and -1 represents the sample label of subjects/ samples for stress and control condition respectively.}
  \item{n}{n is a numeric constant represents the number of informative genes to be selected.}
  }
\value{An informative geneset is obtained, which is relevanit to the particular trait/condition and the genes within the selected geneset are minimum redundant using MRMR technique.}

\references{Ding, C and Peng, H (2005). Minimum redundancy feature selection from microarray gene expression data. J. Bioinformatics Comput Biol 3(2):185-205.}
\author{Samarendra Das}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
n=20
mrmr.cutoff(x, y, n)
}
\keyword{gene}
\keyword{weights}
\keyword{informative geneset}
\keyword{MRMR}