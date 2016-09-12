\name{topsis.meth}
\alias{topsis.meth}
\title{Selection of optimal gene selection method(s)/method(s) through multi-criteria decision analysis}
\description{The function enables to rank gene selection methods/method(s) under a multi-criteria decision making set up and further selection of optimum gene selection method using Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS) approach.}
\usage{topsis.meth(x)}
\arguments{
  \item{x}{x is a M by C data frame representing the performance value of the methods under different criteria, where rows are the methods and columns are the criteria. The performance of the methods are adjudged based on mangitude of the criteria (i.e. higher the criteria value better is the method).}
  }
\value{The function returns a data frame consisting rows as method names and different columns with distance score for positive, negative ideal solution, TOPSIS score and ranks for respective methods.}

\author{Samarendra Das }
\references{Ahn BS (2011) Compatible weighting method with rank order centroid: Maximum entropy ordered weighted averaging approach. Eur J Oper Res 212: 552-559.}
\examples{

x=matrix(runif(150), 10, 15)
rownames(x)=paste("Method",1:nrow(x), sep="")
colnames(x)=paste("C",1:ncol(x), sep="")
x=as.data.frame(x)
topsis.meth(x)
}
\keyword{gene}
\keyword{method}
\keyword{TOPSIS}
\keyword{MCDM}