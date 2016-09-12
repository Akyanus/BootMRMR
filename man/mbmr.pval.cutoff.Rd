\name{mbmr.pval.cutoff}
\alias{mbmr.pval.cutoff}
\title{Selection of informative geneset based on statistical significance value using Modified Bootstrap MRMR technique}
\description{The informative geneset which has maximum relevance with target class/trait and minimum redundancy among genes are selected based on p-values obtained from Modified Bootstrap MRMR technique.}
\usage{mbmr.pval.cutoff(x, y, m, s, Q, n)}
\arguments{
  \item{x}{x is a N by p data frame of gene expression values where rows represent genes and columns represent samples/subject/time point. Each cell entry represents the expression level of a gene in a sample/subject (row names of x as gene names/gene ids).}
  \item{y}{y is a p by 1 numeric vector with entries 1/-1 representing sample labels, where 1/-1 represents the sample label of subjects/ samples for stress/control condition  (for two class problems).}
  \item{m}{m is a scalar representing the size of the Modified Bootstrap Sample (i.e. Out of p samples/subjects, m samples/subjects are randomly drawn with replacement, which constitutes one Modified Bootstrap Sample).}
  \item{s}{s is a scalar representing the number of Modified Bootstrap samples (i.e. number of times each of the m samples/subjects will be resampled from p samples/subjects). }
  \item{Q}{Q is a scalar representing the quartile value of the gene rankscores (lies within 1/N to 1), usually the second quartile, i.e. 0.5 or third quartile i.e. 0.75}
  \item{n}{n is a scalar representing the size of the informative gene set to be obtained.}
  }
\value{The function returns a list of the genes/ geneset which are highly informative to the particular trait/condition under investigation using Modified Bootstrap MRMR technique.}

\author{Samarendra Das}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
m=36
s=80
Q=0.5
n=20
mbmr.pval.cutoff(x, y, m, s, Q, n)
}
\keyword{gene}
\keyword{rankscore}
\keyword{bootstrap}
\keyword{p-value}
\keyword{informative geneset}