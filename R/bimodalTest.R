
bimodalTest <- function(se, cvExpr, threshold = NULL){
  
  library(SummarizedExperiment)
  library(modes)
  
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
  datCounts <- data.frame(t(assay(se)[genes, ]))
  
  #### Calulate the mean z score for each sample by calulating the z score for each gene in each sample
  datZscores <- scale(x = datCounts)
  ZscoreMeans <- rowMeans(x = datZscores) 
  datZdensity <- density(x = ZscoreMeans)
  
  ### Rename the title of the plot to be more informative
  datZdensity$call <- "Density plot of Z-scores from gene expression"
  
  ### Calcualte the biomodal amplitude, coeffeicent and ratio
  bimod_amp <- bimodality_amplitude(ZscoreMeans, fig = FALSE)
  bimod_coef <- bimodality_coefficient(t(datCounts))
  bimod_ratio <- bimodality_ratio(x = ZscoreMeans, fig = FALSE)
  
  ### Group the output into a list to return to the user
  bimodalTestOut <- list("Bimodality.Amplitude" = bimod_amp,
                         "Bimodality.Coefficient" = bimod_coef, 
                         "Bimodality.Ratio" = bimod_ratio,
                         "Density.Plot" = datZdensity)
  return(bimodalTestOut)
  
}  