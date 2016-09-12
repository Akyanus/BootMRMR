\name{rice_salt}
\alias{rice_salt}
\docType{data}
\title{
A gene expression dataset of rice under salinity stress
}
\description{
This data has gene expression values of 200 genes over 40 samples/subjects for a salinity vs. control study in rice. These 40 samples belong to either of salinity stress or control condition (two class problem). This gene expression data is balanced type as the first 20 samples are under salinity stress and the later 20 samples are under control condition.
The first row of the data contains the samples/subjects labels with entries are 1 and -1, where the label '1' and '-1' represent samples generated under salinity stress and control condition respectively.
}
\usage{data("rice_salt")}
\format{
  A data frame with 200 rows as genes with 40 columns as samples/subjects. Each column (sample) represent the gene expression values of genes. Each column as microarray samples with labels -1 or 1 represents control or salinity stress respectively. 
}
\details{
The data is created by taking 200 genes from the large number of genes from NCBI GEO database. The rows are the genes and columns are the samples/subjects. The first half of the samples/subjects are generated under salinity stress condition and other half under control condition.The first row of the data contains the samples/subjects albels with entries are 1 and -1, where th label '1' and '-1' represent samples generated under salinity stress and control condition respectively.}
\source{
Gene Expression Omnibus: NCBI gene expression and hybridization array data repository.ncbi.nlm.nih.gov/geo/.
}
\examples{
data(rice_salt)
}
\keyword{datasets}
\keyword{rice}
\keyword{salt}