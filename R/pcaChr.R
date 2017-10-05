
pcaChr <- function(se, chr){
  
  library(SummarizedExperiment)
  library(FactoMineR)
  library(factoextra)
  
  ### Unit tests to see if the inputted data is in the correct format
  #### se must be a RangedSummarizedExperiment
  if(!is(se, "RangedSummarizedExperiment")){
    stop("se must be a RangedSummarizedExperiment")
  }
  
  ### Subset se for genes on the chromosome of interest and extract the 
  datCounts <- data.frame(t(assay(se[seqnames(se) == chr])))
  
  ### Transform datCounts into a matrix with a column containing the clustering IDs
  ### In order to run the pca
  datCounts$Group <- colData(se)$Ploidy
  gpCol = which(colnames(datCounts) == "Group")
  
  ## Run the pca
  pca <- PCA(datCounts, quali.sup=gpCol, graph = FALSE)
  
  ### Output the plot 
  fviz_pca_ind(X = pca, label="none", habillage=factor(colData(se)$Ploidy), title="",
               addEllipses=TRUE, ellipse.level=0.95)
  
}
