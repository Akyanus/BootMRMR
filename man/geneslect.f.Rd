\name{geneslect.f}
\alias{geneslect.f}
\title{Informative gene set selection using F-score}
\description{The function returns geneset which is informative for a particular trait/condition using F-score as the gene selection criterion.}
\usage{geneslect.f(x, y, s)}
\arguments{
  \item{x}{x is a N by p dataframe of gene expression, where, rows represent as genes and columns as samples/subjects (with row names as gene names/ids).}
  \item{y}{y is a p by 1 numeric vector of 1 and -1, where 1/-1 indicates the class label of the samples/subjects either of two classes (e.g. stress and control).}
  \item{s}{s is a numeric constant (< N) representing the number of genes to be selected from the large gene space.}
 }

\author{Samarendra Das}
\details{This function identifies the genes/ geneset which is informative for the particular trait/condition using F-score as a criterion.}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
s=20
geneslect.f (x, y, s)
}
\keyword{gene}
\keyword{gene subset}
\keyword{F-score}