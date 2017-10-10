
cvExpr <- function(se, chr){
  
  library(SummarizedExperiment)
  library(matrixStats)
  
  ### Unit tests to see if the inputted data is in the correct format
  #### se must be a RangedSummarizedExperiment
  if(!is(se, "RangedSummarizedExperiment")){
    stop("se must be a RangedSummarizedExperiment")
  }
  
  ### Subset se for genes on the chromosome of interest 
  datCounts <- assay(se[seqnames(se) == chr])
  
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
