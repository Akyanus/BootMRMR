\name{bmrmr.pval.cutoff}
\alias{bmrmr.pval.cutoff}
\title{Selection of informative geneset based on statistical significance value using Bootstrap-MRMR technique}
\description{The informative geneset which has maximum relevance with target class/trait and minimum redundancy among genes based on statistical significance values computed from the Bootstrap-MRMR technique.}
\usage{bmrmr.pval.cutoff(x, y, s, Q, n)}
\arguments{
  \item{x}{x is a N by p data frame of gene expression values where rows represent genes and columns represent samples or subjects or time point. Each cell entry represents the expression level of a gene in a sample/subject (row names of x as gene names or gene ids).}
  \item{y}{y is a p by 1 numeric vector with entries 1 or -1 representing sample labels, where, 1/-1 represents the sample label of subjects/ samples for stress/control condition(for two class problems).}
  \item{s}{s is a scalar representing the number of bootstrap generated, s must be sufficiently large (i.e. number of times  bootstrap samples are generated.}
  \item{Q}{Q is a scalar representing the quartile value of the rankscores of genes (lies within 1/N to 1), usually the second quartile, i.e. 0.5 or third quartile i.e. 0.75 may be taken.}
  \item{n}{n is a scalar representing the size of the informative geneset to be obtained.}
  }
\value{The function returns a list of the genes/informative geneset which are highly relevant to the particular trait/condition under investigation and minimal redundant among themselves.}

\author{Samarendra Das}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
s=80
Q=0.5
n=20
bmrmr.pval.cutoff(x, y, s, Q, n)
}
\keyword{gene}
\keyword{rankscore}
\keyword{bootstrap}
\keyword{p-value}
\keyword{informative geneset}