#' bimodalTest
#'
#' Performs a bimodal test over a genomic region of interest
#'
#' @param se A SummarizedExperiment containing the gene expression data.
#' @param cvExprA The output from cvExpr.
#' @param threshold Optional. The quantile threshold of genes to be used for clustering analaysis. Default is NULL.
#' @usage bimodalTest(se, cvExpr, threshold = NULL)
#' @return Returns a list containing the output from the bimodal test.
#' @import SummarizedExperiment
#' @import modes
#' @import diptest
#' @author Benjamin Mayne
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
