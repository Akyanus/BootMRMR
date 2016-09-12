\name{mbmr.weight.cutoff}
\alias{mbmr.weight.cutoff}
\title{Identification of informative geneset based on weights obtained from Modified Bootstrap-MRMR technique}
\description{The function enables to find set of informative genes based on weights which are obtained by maximising the relevancy of genes with classes/condition/trait and minimising the redundancy among genes using Modified Bootstrap-MRMR technique}
\usage{mbmr.weight.cutoff(x, y, m, s, n)}
\arguments{
   \item{x}{x is a N by p data frame of gene expression values where rows represent genes and columns represent samples/subject/time point. Each cell entry represents the expression level of a gene in a sample/subject (row names of x as gene names/gene ids).}
 \item{y}{y is a p by 1 numeric vector with entries 1/-1 representing sample labels, where 1/-1 represents the sample label of subjects/ samples for stress/control condition  (for two class problems).}
  \item{m}{m is a scalar representing the size of the Modified Bootstrap Sample (i.e. Out of p samples/subjects, m samples/subjects are randomly drawn with replacement, which constitutes one Modified Bootstrap Sample).}
  \item{s}{s is a scalar representing the number of Modified Bootstrap samples (i.e. number of times each of the m samples/subjects will be resampled from p samples/subjects). }
  \item{n}{n is a numeric constant representing the number of informative genes to be selected from the large gene space.}
}
\value{The function returns a set of genes, which are highly informative to the trait or condition under consideration based Modified Bootstrap-MRMR weights.}

\author{Samarendra Das}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
m=36
s=80
n=20
mbmr.weight.cutoff(x, y, m, s, n)
}
\keyword{gene}
\keyword{weights}
\keyword{informative geneset}
\keyword{Modified Bootstrap-MRMR}