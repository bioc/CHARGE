#' cvExpr
#'
#' Calculates the coefficient of variation for each gene within a defined genomic region.
#'
#' @param se A SummarizedExperiment containing the normalised gene expression data.
#' @param region A GRanges object containing the genomic location of the region of interest. This can either be an entire length or the subset of a chromosome.
#' @usage cvExpr(se, region)
#' @return Returns a list containing the CV of each gene and the the quantile threshold of the data.
#' @import SummarizedExperiment
#' @importFrom stats quantile
#' @importFrom methods is
#' @importFrom matrixStats rowSds
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps
#' @author Benjamin Mayne
#' @examples
#' library(GenomicRanges)
#' data(datExprs)
#' chr21 <- GRanges("21:1-46709983")
#' cvExpr.out <- cvExpr(se = datExprs, region = chr21)
#' @details
#' Calculates the coefficient of variation (CV) for each gene with a genomic region of interest.
#' The CV Values can be used to determine which genes are not critical for appropiate clustering and can be filtered out prior to clustering.
#' @export

cvExpr <- function(se, region){

  ### Unit tests to see if the inputted data is in the correct format
  #### se must be a RangedSummarizedExperiment
  if(!is(se, "RangedSummarizedExperiment")){
    stop("se must be a RangedSummarizedExperiment")
  }

  #### region must be a GRanges object
  if(!is(region, "GRanges")){
    stop("region must be a GRanges")
  }

  ### Subset se for genes over the region of interest
  datCounts <- assay(subsetByOverlaps(se, region))

  ### Determine the coefficient of variation (CV) for gene expression from the chromosome of interest
  #### Firstly calculate the the standard deviation of each gene
  geneSds <- rowSds(datCounts)

  #### Calucalte the mean of each gene
  geneMeans <- rowMeans(datCounts)

  #### Determine CV = 100*sd/mean
  geneCV <- 100*geneSds/geneMeans

  ### Determine the genes quantiles
  CVQuantiles <- quantile(geneCV)

  ### Return the gene's CV and the quantiles to be inputted into the next function
  ### which will determine what genes to be used for clustering
  return(list(ExpressionCV = geneCV, Quantile_Threshold = CVQuantiles))


}
