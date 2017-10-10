
clusterChr <- function(se, cvExpr, threshold = NULL){
  
  library(SummarizedExperiment)
  
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
  datCounts <- assay(se)[genes, ]
  
  ### Perform a k-means clustering using two clusters (euploidy vs aneuploidy) 
  kmeans.Out <- kmeans(x = t(datCounts), centers = 2)
  
  ### Next determine the mean level of expression of each cluster
  Cluster1_Mean <- mean(colMeans(datCounts[ ,which(kmeans.Out$cluster == 1)]))
  Cluster2_Mean <- mean(colMeans(datCounts[ ,which(kmeans.Out$cluster == 2)]))
  
  ### Extract the clustering labels for each sample and assign if ploidy state with respect to the other sample
  ### based on the clustering means
  pd <- data.frame(kmeans.Out$cluster)
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
  
