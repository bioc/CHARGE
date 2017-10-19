#' pcaExprs 
#'
#' Determines the coefficient of variation (CV) for gene expression 
#'
#' @param se A SummarizedExperiment containing the gene expression data and the clustering output from clusterExpr.
#' @param cvExprA The output from cvExpr.
#' @param threshold Optional. The quantile threshold of genes to be used for clustering analaysis. Default is NULL.
#' @usage pcaExprs(se, cvExpr, threshold = NULL)
#' @return Returns a PCA plot showing the seperation of samples that were labelled hyperploidy or hypoploidy in clusterExpr.
#' @import SummarizedExperiment
#' @import factoextra
#' @import FactoMineR
#' @author Benjamin Mayne
#' @export

pcaExpr <- function(se, cvExpr, threshold = NULL){
  
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
  
  ### Transform datCounts into a matrix with a column containing the clustering IDs
  ### In order to run the pca
  datCounts$Group <- colData(se)$Ploidy
  gpCol = which(colnames(datCounts) == "Group")
  
  ## Run the pca
  pca <- PCA(datCounts, quali.sup=gpCol, graph = FALSE)
  
  ### Output the plot 
  fviz_pca_ind(X = pca, label="none", habillage=factor(colData(se)$Ploidy), title="",
               addEllipses=TRUE, ellipse.level=0.5)

}