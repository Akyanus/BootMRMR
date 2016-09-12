#########################################BootMRMR package##############################

##############################Dependent Packages#####################################
library(stats)
#######################Gene selection using F-score algorithm ########################
geneslect.f <- function(x, y, s)
{
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  cls <- as.vector(y) ###class information###
  genenames <- rownames(x)
  g <- as.matrix(x)
  n <- nrow(g)               ###number of genes####
  M <- ncol(x)          ###number of samples###
  idx <- which(cls==1)  # indexing of positive samples
  idy <- which(cls==-1) # indexing of negative samples
  B=vector(mode="numeric", n)
  for(i in 1:n){
    f.mes <-(((mean(g[i, idx])-mean(g[i, ]))^2)+ ((mean(g[i, idy])-mean(g[i, ]))^2))/(var(g[i, idx])+var(g[i, idy]))  #####F-Score    
    B[i] <- f.mes
  }
  ranking <- sort(-B, index.return = TRUE)$ix 
  temp <- ranking[1:s]
  selectgenes <- genenames[temp]
  class(selectgenes) <- "Informative geneset"
  return(selectgenes)
}

#######################################################################################
#                       Gene selection using MRMR algorithm                           #
#######################################################################################


Weights.mrmr <- function(x, y)
{      ###x is gene expression data matrix, y is vector of class labels###
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  
  cls <- as.numeric(y)                  ###class information###
  g <- as.matrix(x)
  genes <- rownames(x)
  n <- nrow(g)                          ###number of genes####
  m <- ncol(x)                          ###number of samples###
  GeneRankedList <- vector(length=n)
  qsi <- as.vector((apply(abs(cor(t(g), method="pearson",use="p")-diag(n)), 1, sum))/(n-1))
  idx <- which(cls==1)                  ###indexing of positive samples
  idy <- which(cls==-1)                 ###indexing of negative samples
  B <- vector(mode="numeric", n)
  for(i in 1:n){
    f.mes <-(((mean(g[i, idx])-mean(g[i, ]))^2)+ ((mean(g[i, idy])-mean(g[i, ]))^2))/(var(g[i, idx])+var(g[i, idy]))  #####F-Score    
    B[i] <- f.mes
  }
  rsi <- abs(B)
  rankingCriteria <- rsi/qsi
  names(rankingCriteria) <- genes
  class(rankingCriteria) <- "MRMR weights"
  return(rankingCriteria)
}

################################Gene subset selection##################################

mrmr.cutoff <- function (x, y, n)
{
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(n < 0 & n > nrow(x))
  {
    warning(" n must be numeric and less than total number of genes in gene space")
  }
  weights <- as.vector(Weights.mrmr (x, y))
  gene.id <- sort(-weights, index.return=TRUE)$ix
  temp <- gene.id [1:n]
  genes <- rownames(x)
  select.gene <- genes[temp]
  class(select.gene) <- "Informative geneset"
  return (select.gene)
}
####################################Ends here#########################################

######################################################################################
#                  Boot-MRMR technique of gene selection                             #
######################################################################################

######################Gene Selection using Modified Bootstrap procedure###############

weight.mbmr <- function(x, y, m, s, plot=TRUE) ### x is the gene expression data matrix, y is class vector, m is bootstrap sample size, n is number of bootstraps tob e drawn#
{            
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(m < 0 & m > ncol(x))
  {
    warning("m must be numeric and less than total number of samples in the input GE data")
  }
  if(s < 0 & s <= 50)
  {
    warning("s must be numeric and sufficiently large")
  }
  cls <- as.numeric(y)               ### class information###
  genes <- rownames(x)
  g <- as.matrix(x)
  nn <- nrow(g)                       ### number of genes####
  M <- ncol(x)                       ### number of samples###
  #if(missing(m)) 
  #m <- trunc(0.9*M)                  ### Fix the size of bootstrap sample###
  GeneRankedList <- vector(length=nn)
  M1 <- matrix(0, nn, s) 
  
  for (j in 1:s) {
    samp <- sample(M, m, replace=TRUE) ###select bootstrap sample #####
    x1 <- g[, samp]
    y1 <- cls[samp] 
    qsi <- as.vector((apply(abs(cor(t(x1), method="pearson",use="p")-diag(nrow(x1))), 1, sum))/(nrow(x1)-1))
    idx <- which(y1==1)                # indexing of positive samples
    idy <- which(y1==-1)               # indexing of negative samples
    B=vector(mode="numeric", nn)
    for(i in 1:nrow(x1)){
      f.mes <-(((mean(x1[i, idx])-mean(x1[i, ]))^2)+ ((mean(x1[i, idy])-mean(x1[i, ]))^2))/(var(x1[i, idx])+var(x1[i, idy]))  #####F-Score    
      B[i] <- f.mes
    }
    rsi <- abs(B)
    rankingCriteria <- rsi/qsi
    GeneRankedList <- sort(-rankingCriteria, index.return = TRUE)$ix
    rankvalue <- sort(GeneRankedList, index.return=TRUE)$ix
    rankscore <- (nn+1-rankvalue)/(nn)
    M1[,j] <- as.vector(rankscore)
    #print(rankscore)
  }
  rownames(M1) <- genes
  rankingcriteria <- as.vector(rowSums((M1), na.rm = FALSE, dims = 1))
  names(rankingcriteria) <- genes
  if (plot==TRUE)
  {
    wgs <- as.vector(rankingcriteria)
    ranks <- sort(wgs, decreasing=TRUE, index.return=TRUE)$ix
    wgs_sort <- wgs[ranks]
    plot(1:length(ranks), wgs_sort, main="Gene selection plot", type="p", xlab="Genes", ylab="Weights")
  }
  class( rankingcriteria) <- "MBootMRMR weights"
  return(rankingcriteria)
}

#####################################Gene set selection using Boot-MRMR weights########

mbmr.weight.cutoff <- function (x, y, m, s, n)
{
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(m < 0 & m > ncol(x))
  {
    warning("m must be numeric and less than total number of samples in the input GE data")
  }
  if(s < 0 & s <= 50)
  {
    warning("s must be numeric and sufficiently large")
  }
  if(n > nrow(x))
  {
    stop("Number of informative genes to be selected must be less than total number of genes")
  }
  weights <- weight.mbmr(x, y, m, s, plot=FALSE)
  genes <- rownames(x)
  w1 <- as.vector(weights)
  gene.id <- sort(-w1, index.return=TRUE)$ix
  temp <- gene.id [1:n]
  select.gene <- genes[temp]
  class(select.gene) <- "Informative geneset"
  return (select.gene)
}

####################Computation of Boot-MRMR statistical significance values###########

pval.mbmr <- function(x, y, m, s, Q, plot=TRUE)
{ 
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(m < 0 & m > ncol(x))
  {
    warning("m must be numeric and less than total number of samples in the input GE data")
  }
  if(s < 0 & s <= 50)
  {
    warning("s must be numeric and sufficiently large")
  }
  if(Q < 0 & Q > 1)
  {
    warning("Q is the quartile value of rank scores and must be within 0 and 1")
  }
  if (missing (Q))
  {
    Q <- 0.5
  }
  cls <- as.numeric(y)               ### class information###
  genes <- rownames(x)
  g <- as.matrix(x)
  n1 <- nrow(g)                       ### number of genes####
  M <- ncol(x)                       ### number of samples###
  if(missing(m)) 
  {
    m <- trunc(0.9*M)
  }### Fix the size of bootstrap sample###
  GeneRankedList <- vector(length=n1)
  M1 <- matrix(0, n1, s) 
  ##if(missing(s)) 
  for (j in 1:s) {
    samp <- sample(M, m, replace=TRUE) ###select bootstrap sample #####
    x1 <- g[, samp]
    y1 <- cls[samp] 
    qsi <- as.vector((apply(abs(cor(t(x1), method="pearson",use="p")-diag(nrow(x1))), 1, sum))/(nrow(x1)-1))
    idx <- which(y1==1)                # indexing of positive samples
    idy <- which(y1==-1)               # indexing of negative samples
    B=vector(mode="numeric", n1)
    for(i in 1:nrow(x1)){
      f.mes <-(((mean(x1[i, idx])-mean(x1[i, ]))^2)+ ((mean(x1[i, idy])-mean(x1[i, ]))^2))/(var(x1[i, idx])+var(x1[i, idy]))  #####F-Score    
      B[i] <- f.mes
    }
    rsi <- abs(B)
    rankingCriteria <- rsi/qsi
    GeneRankedList <- sort(-rankingCriteria, index.return = TRUE)$ix
    rankvalue <- sort(GeneRankedList, index.return=TRUE)$ix
    rankscore <- (n1+1-rankvalue)/(n1)
    M1[,j] <- as.vector(rankscore)
  }
  rankscore <- as.matrix(M1)          # Raw scores obtained from SVM-RFE
  mu <- Q                                    # value under null hypothesis
  #if (missing (Q))
  R <- rankscore - mu                       # Transformed score of MRMR under H0
  sam <- nrow (R)                           # number of genes
  pval.vec <- vector(mode="numeric", length=nrow(rankscore))
  for (i in 1:sam) {
    z <- R[i,]
    z <- z[z != 0]
    n11 <- length(z)
    r <- rank(abs(z))
    tplus <- sum(r[z > 0])
    etplus <- n11 * (n11 + 1) / 4
    vtplus <- n11 * (n11 + 1) * (2 * n11 + 1) / 24
    p.value=pnorm(tplus, etplus, sqrt(vtplus), lower.tail=FALSE)
    pval.vec[i]=p.value
  }
  names( pval.vec) <- genes
  if (plot==TRUE)
  {
    pval <- as.vector(pval.vec)
    ranks <- sort(pval, decreasing=FALSE, index.return=TRUE)$ix
    pval_sort <- pval[ranks]
    plot(1:length(ranks), pval_sort, main="Gene selection plot", type="p", xlab="Genes", ylab="p-values")
  }
  class(pval.vec) <- "p value"
  return(pval.vec)
}

######################Gene set selection using statistical significance values ########

mbmr.pval.cutoff <- function (x, y, m, s, Q, n)
{
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(m < 0 & m > ncol(x))
  {
    warning("m must be numeric and less than total number of samples in the input GE data")
  }
  if(s < 0 & s <= 50)
  {
    warning("s must be numeric and sufficiently large")
  }
  if(Q < 0 & Q > 1)
  {
    warning("Q is the quartile value of rank scores and must be within 0 and 1")
  }
  if (missing (Q))
  {
    Q <- 0.5
  }
  if(n > nrow(x))
  {
    stop("Number of informative genes to be selected must be less than total number of genes")
  }
  pvalue <- pval.mbmr (x, y, m, s, Q, plot=FALSE)
  genes <- rownames(x)
  w11 <- as.vector(pvalue)
  gene.id <- sort(w11, index.return=TRUE)$ix
  temp <- gene.id [1:n]
  select.gene <- genes[temp]
  class(select.gene) <- "Informative geneset"
  return (select.gene)
}

######################Gene Selection using Standard Bootstrap procedure###############

bootmr.weight <- function(x, y, s, plot=TRUE) ### x is the gene expression data matrix, y is class vector, m is bootstrap sample size, n is number of bootstraps tob e drawn#
{            
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(s < 0 & s <= 50)
  {
    warning("s must be numeric and sufficiently large")
  }
  cls <- as.numeric(y)               ### class information###
  genes <- rownames(x)
  g <- as.matrix(x)
  n1 <- nrow(g)                       ### number of genes####
  M <- ncol(x)                       ### number of samples###
  
  GeneRankedList <- vector(length=n1)
  M1 <- matrix(0, n1, s) 
  ##if(missing(s)) 
  for (j in 1:s) {
    samp <- sample(M, M, replace=TRUE) ###select bootstrap sample #####
    x1 <- g[, samp]
    y1 <- cls[samp] 
    qsi <- as.vector((apply(abs(cor(t(x1), method="pearson",use="p")-diag(n1)), 1, sum))/(n1-1))
    idx <- which(y1==1)                # indexing of positive samples
    idy <- which(y1==-1)               # indexing of negative samples
    B=vector(mode="numeric", n1)
    for(i in 1:nrow(x1)){
      f.mes <-(((mean(x1[i, idx])-mean(x1[i, ]))^2)+ ((mean(x1[i, idy])-mean(x1[i, ]))^2))/(var(x1[i, idx])+var(x1[i, idy]))  #####F-Score    
      B[i] <- f.mes
    }
    rsi <- abs(B)
    rankingCriteria <- rsi/qsi
    GeneRankedList <- sort(-rankingCriteria, index.return = TRUE)$ix
    rankvalue <- sort(GeneRankedList, index.return=TRUE)$ix
    rankscore <- (n1+1-rankvalue)/(n1)
    M1[,j] <- as.vector(rankscore)
  }
  
  rankingcriteria <- as.vector(rowSums((M1), na.rm = FALSE, dims = 1))
  names(rankingcriteria) <- genes
  if (plot==TRUE)
  {
    wgs <- as.vector(rankingcriteria)
    ranks <- sort(wgs, decreasing=TRUE, index.return=TRUE)$ix
    wgs_sort <- wgs[ranks]
    plot(1:length(ranks), wgs_sort, main="Gene selection plot", type="p", xlab="Genes", ylab="Weights")
  }
  class(rankingcriteria) <- "Boot-MRMR weights"
  return(rankingcriteria)
}

#####################################Gene set selection using Boot-MRMR weights########

bmrmr.weight.cutoff <- function (x, y, s, n)
{
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(s < 0 & s <= 50)
  {
    warning("s must be numeric and sufficiently large")
  }
  if(n > nrow(x))
  {
    stop("Number of informative genes to be selected must be less than total number of genes")
  }
  weights <- bootmr.weight(x, y, s, plot=FALSE)
  genes <- names(weights)
  w1 <- as.vector(weights)
  gene.id <- sort(-w1, index.return=TRUE)$ix
  temp <- gene.id [1:n]
  select.gene <- genes[temp]
  class(select.gene) <- "Informative geneset"
  return (select.gene)
}

####################Computation of Boot-MRMR statistical significance values###########

pval.bmrmr <- function(x, y, s, Q, plot=TRUE)
{ 
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(s < 0 & s <= 50)
  {
    warning("s must be numeric and sufficiently large")
  }
  if(Q < 0 & Q > 1)
  {
    warning("Q is the quartile value of rank scores and must be within 0 and 1")
  }
  if (missing (Q))
  {
    Q <- 0.5
  }
  cls <- as.numeric(y)               ### class information###
  genes <- rownames(x)
  g <- as.matrix(x)
  n1 <- nrow(g)                       ### number of genes####
  M <- ncol(x)                       ### number of samples###
  GeneRankedList <- vector(length=n1)
  M1 <- matrix(0, n1, s) 
  ##if(missing(s)) 
  for (j in 1:s) {
    samp <- sample(M, M, replace=TRUE) ###select bootstrap sample #####
    x1 <- g[, samp]
    y1 <- cls[samp] 
    qsi <- as.vector((apply(abs(cor(t(x1), method="pearson",use="p")-diag(n1)), 1, sum))/(n1-1))
    idx <- which(y1==1)                # indexing of positive samples
    idy <- which(y1==-1)               # indexing of negative samples
    B=vector(mode="numeric", n1)
    for(i in 1:nrow(x1)){
      f.mes <-(((mean(x1[i, idx])-mean(x1[i, ]))^2)+ ((mean(x1[i, idy])-mean(x1[i, ]))^2))/(var(x1[i, idx])+var(x1[i, idy]))  #####F-Score    
      B[i] <- f.mes
    }
    rsi <- abs(B)
    rankingCriteria <- rsi/qsi
    GeneRankedList <- sort(-rankingCriteria, index.return = TRUE)$ix
    rankvalue <- sort(GeneRankedList, index.return=TRUE)$ix
    rankscore <- (n1+1-rankvalue)/(n1)
    M1[,j] <- as.vector(rankscore)
  }
  rankscore <- as.matrix(M1)          # Raw scores obtained from SVM-RFE
  mu <- Q                                    # value under null hypothesis
  
  R <- rankscore - mu                       # Transformed score of MRMR under H0
  sam <- nrow (R)                           # number of genes
  pval.vec <- vector(mode="numeric", length=nrow(rankscore))
  for (i in 1:sam) {
    z <- R[i,]
    z <- z[z != 0]
    n11 <- length(z)
    r <- rank(abs(z))
    tplus <- sum(r[z > 0])
    etplus <- n11 * (n11 + 1) / 4
    vtplus <- n11 * (n11 + 1) * (2 * n11 + 1) / 24
    p.value=pnorm(tplus, etplus, sqrt(vtplus), lower.tail=FALSE)
    pval.vec[i]=p.value
  }
  names( pval.vec) <- genes
  if (plot==TRUE)
  {
    pval <- as.vector(pval.vec)
    ranks <- sort(pval, decreasing=FALSE, index.return=TRUE)$ix
    pval_sort <- pval[ranks]
    plot(1:length(ranks), pval_sort, main="Gene selection plot", type="p", xlab="Genes", ylab="p-values")
  }
  class(pval.vec) <- "p values"
  return(pval.vec)
}

######################Gene set selection using statistical significance values ########

bmrmr.pval.cutoff <- function (x, y, s, Q, n)
{
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }   
  if ((!class(y)=="numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }   
  if (!length(y)==ncol(x))
  {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if(s < 0 & s <= 50)
  {
    warning("s must be numeric and sufficiently large")
  }
  if(Q < 0 & Q > 1)
  {
    warning("Q is the quartile value of rank scores and must be within 0 and 1")
  }
  if (missing (Q))
  {
    Q <- 0.5
  }
  if(n > nrow(x))
  {
    stop("Number of informative genes to be selected must be less than total number of genes")
  }
  pvalue <- pval.bmrmr(x, y, s, Q, plot=FALSE)
  genes <- names(pvalue)
  w11 <- as.vector(pvalue)
  gene.id <- sort(w11, index.return=TRUE)$ix
  temp <- gene.id [1:n]
  select.gene <- genes[temp]
  class(select.gene) <- "Informative geneset"
  return (select.gene)
}

########################Booot-MRMR Ends Here##########################################

#####################################################################################
#          Technique for Order of Preference by Similarity to Ideal Solution        # 
#####################################################################################
topsis.meth <- function (x)
{
  if(!class(x)=="data.frame")
  {
    warning("x must be a data frame and rows as methods and columns as criteria")
  }
  dat <- as.matrix(x)
  methods <- rownames(x)
  criteria <- colnames(x)
  #####Calculation of weights##############
  dat.sum <- apply(dat, 2, sum)
  dat.norm <- t(t(dat)/(dat.sum))  #####Decission matrix to Normalized mode#####
  dat.norm1 <- dat.norm*log(dat.norm)
  m <- nrow(dat)  ###number of methods###
  n <- ncol(dat) ### number of criterias###
  a <- 1/log(m)
  dat.entr <- (-a)*apply(dat.norm1, 2, sum)
  dat.diver <- 1-dat.entr
  weights <- dat.diver/sum(dat.diver)
  ##########Topsis MCDM############
  dat.seq <- dat^2
  dat.sum <- (apply(dat.seq, 2, sum))^0.5
  R.norm <- t(t(dat)/(dat.sum))
  V.norm <- t(t(R.norm)*weights)
  V.pos <- vector(mode="numeric", length=n)
  V.neg <- vector(mode="numeric", length=n)
  for (j in 1:n){
    V.pos[j] <- max(V.norm[,j])
  }
  for (j in 1:n){
    V.neg[j] <- min(V.norm[,j])
  }
  D.p <- t(t(V.norm)-V.pos)
  
  d.pos <- (apply((D.p)^2, 1, sum))^0.5
  D.n <- t(t(V.norm)-V.neg)
  d.neg <- (apply((D.n)^2, 1, sum))^0.5
  d.tot <- d.pos+d.neg
  C.Score <-d.neg/(d.pos+d.neg)
  C.rank <- rank(-C.Score)
  D.topsis <- cbind(d.pos, d.neg, C.Score, C.rank)
  colnames(D.topsis) <- c("di+", "di-", "Topsis Score", "Rank")
  rownames(D.topsis) <- methods
  class(D.topsis) <- "TOPSIS"
  return(D.topsis)
}
###############Ends here###############################################