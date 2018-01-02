#' clusterExpr
#'
#' Performs a clustering analysis and predict which samples have an enrichment for a genomic region of interest.
#'
#' @param se A SummarizedExperiment containing the gene expression data.
#' @param cvExpr The output from cvExpr.
#' @param threshold Optional. The quantile threshold of genes to be used for clustering analaysis. Default is NULL.
#' @usage clusterExpr(se, cvExpr, threshold = NULL)
#' @return Returns a SummarizedExperiment containing the original inputted se, but where an additional column labelled Ploidy has been added into the meta data containing the classification of each sample.
#' @import SummarizedExperiment
#' @import cluster
#' @importFrom methods is
#' @author Benjamin Mayne
#' @examples
#' library(GenomicRanges)
#' data(datExprs)
#' chr21 <- GRanges("21:1-46709983")
#' cvExpr.out <- cvExpr(se = datExprs, region = chr21)
#' datExprs <- clusterExpr(se = datExprs, cvExpr = cvExpr.out, threshold = "25%")
#' colData(datExprs)$Ploidy
#' @export

clusterExpr <- function(se, cvExpr, threshold = NULL){

  ### Unit tests to see if the inputted data is in the correct format
  #### se must be a RangedSummarizedExperiment
  if(!is(se, "RangedSummarizedExperiment")){
    stop("se must be a RangedSummarizedExperiment")
  }

  ### Subset genes from the cvExpr using the input from threshold
  if(is.null(threshold)){
    #### Use all the gene in the analysis if threshold = NULL
    genes <- names(cvExpr[[1]])
  } else {
    #### Subset the genes based on defined quantile threshold
    genes <- names(which(cvExpr[[1]] > cvExpr[[2]][threshold]))
  }

  ### Subset se for genes
  datExpr <- data.frame(assay(se))[genes, ]

  ### Perform a k-means clustering using two clusters
  pam.out <- pam(x=t(datExpr), k=2)

  ### Next determine the mean level of expression of each cluster
  Cluster1_Mean <- mean(colMeans(datExpr[ ,which(pam.out$clustering == 1)]))
  Cluster2_Mean <- mean(colMeans(datExpr[ ,which(pam.out$clustering == 2)]))

  ### Extract the clustering labels for each sample and assign if ploidy state with respect to the other sample
  ### based on the clustering means
  pd <- data.frame(pam.out$clustering)
  colnames(pd) <- "Ploidy"

    if(Cluster1_Mean > Cluster2_Mean){

    pd$Ploidy <- ifelse(test = pd$Ploidy == 1, yes = "Hyperploidy", no = "Hypoploidy")

  } else {

    pd$Ploidy <- ifelse(test = pd$Ploidy == 2, yes = "Hyperploidy", no = "Hypoploidy")

  }

  ### Extract the meta data from se to combine with the clustering data
  colData(se) <- merge(colData(se), pd, by="row.names", all.x=TRUE)

  return(se)

}
