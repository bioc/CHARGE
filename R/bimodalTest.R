#' bimodalTest
#'
#' Performs a bimodal test and calcualtes Hartigan's statistic and p-value over a genomic region of interest using gene expression data set, the output from cvExpr.
#'
#' @param se A SummarizedExperiment containing the normalised gene expression data.
#' @param cvExpr The output from cvExpr.
#' @param threshold Optional. The quantile threshold of genes to be used for clustering analaysis. Default is NULL.
#' @usage bimodalTest(se, cvExpr, threshold = NULL)
#' @return Returns a list containing the output from the bimodal test.
#' @import SummarizedExperiment
#' @importFrom stats density
#' @importFrom methods is
#' @import modes
#' @import diptest
#' @author Benjamin Mayne
#' @examples
#' library(GenomicRanges)
#' data(datExprs)
#' chr21 <- GRanges("21:1-46709983")
#' cvExpr.out <- cvExpr(se = datExprs, region = chr21)
#' bimodalTest.out <- bimodalTest(se = datExprs, cvExpr = cvExpr.out, threshold = "25%")
#' @details
#' Performs a bimodal test and calculates Hartigan's dip test statistic for unimodality for a given gene expression data set.
#' A bimodality coefficient value > 5/9 suggests bimodality and the closer the bimodality ratio is to 1, the more evenly distrubted the data set.
#' The dip statistic and p-value can be used to determine if the region of interest is statistically significant.
#'
#' The second part of the function returns the Z score means which can be used to visualise the denisty or distribution of the samples.
#' @export

bimodalTest <- function(se, cvExpr, threshold = NULL){

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
  datExpr <- data.frame(t(datExpr))

  #### Calulate the mean z score for each sample by calulating the z score for each gene in each sample
  datZscores <- scale(x = datExpr)
  ZscoreMeans <- rowMeans(x = datZscores)
  datZdensity <- density(x = ZscoreMeans)

  ### Rename the title of the plot to be more informative
  datZdensity$call <- "Density plot of Z-scores from gene expression"

  ### Calcualte the biomodal amplitude, coeffeicent and ratio
  bimodStats <- data.frame(bimodality_coefficient(t(datExpr)))
  bimodStats[,2] <- bimodality_ratio(x = ZscoreMeans, fig = FALSE)

  ### Dip statistic test
  dipResult <- dip.test(x = ZscoreMeans)
  bimodStats[,3] <- as.numeric(dipResult["statistic"])
  bimodStats[,4] <- as.numeric(dipResult["p.value"])
  colnames(bimodStats) <- c("Bimodality.Coefficient", "Bimodality.Ratio", "Dip.Statistic", "Dip.pvalue")

  ### Group the output into a list to return to the user
  bimodalTestOut <- list("Bimodal.Statistics" = bimodStats,
                         "Density.Plot" = datZdensity)

  return(bimodalTestOut)

}
