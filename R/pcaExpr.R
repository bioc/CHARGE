#' pcaExprs
#'
#' Creates a PCA plot using genes within a defined genomic region.
#'
#' @param se A SummarizedExperiment containing the gene expression data and the clustering output from clusterExpr.
#' @param cvExpr The output from cvExpr.
#' @param threshold Optional. The quantile threshold of genes to be used for clustering analaysis. Default is NULL.
#' @usage pcaExpr(se, cvExpr, threshold = NULL)
#' @return Returns a PCA plot showing the seperation of samples that were labelled hyperploidy or hypoploidy in clusterExpr.
#' @import SummarizedExperiment
#' @import factoextra
#' @import FactoMineR
#' @importFrom methods is
#' @author Benjamin Mayne
#' @examples
#' library(CHARGE)
#' library(GenomicRanges)
#' data(datExprs)
#' chr21 <- GRanges("21:1-46709983")
#' cvExpr.out <- cvExpr(se = datExprs, region = chr21)
#' datExprs <- clusterExpr(se = datExprs, cvExpr = cvExpr.out, threshold = "25%")
#' pcaExpr(se = datExprs, cvExpr = cvExpr.out, threshold = "25%")
#' @export


pcaExpr <- function(se, cvExpr, threshold = NULL){

  ### Unit tests to see if the inputted data is in the correct format
  #### se must be a RangedSummarizedExperiment
  if(!is(se, "RangedSummarizedExperiment")){
    stop("se must be a RangedSummarizedExperiment")
  }

  ### Subset genes from the cvExpr using the input from threshold
  if(is.null(threshold)){
    #### Use all the genes in the analysis if threshold = NULL
    genes <- names(cvExpr[[1]])
  } else {
    #### Subset the genes based on defined quantile threshold
    genes <- names(which(cvExpr[[1]] > cvExpr[[2]][threshold]))
  }

  ### Subset se for genes
  datExpr <- data.frame(assay(se))[genes, ]
  datExpr <- data.frame(t(datExpr))

  ### Transform datExpr into a matrix with a column containing the clustering IDs
  ### In order to run the pca
  datExpr$Group <- colData(se)$Ploidy
  gpCol = which(colnames(datExpr) == "Group")

  ## Run the pca
  pca <- PCA(datExpr, quali.sup=gpCol, graph = FALSE)

  ### Output the plot
  fviz_pca_ind(X = pca, label="none", habillage=factor(colData(se)$Ploidy), title="",
               addEllipses=TRUE, ellipse.level=0.95)

}
