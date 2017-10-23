#' plotcvExpr 
#'
#' Plots the coefficient of variation (CV) for gene expression over the region of interest 
#'
#' @param cvExpr The output from cvExpr function.
#' @usage plotcvExpr(cvExpr)
#' @return Returns a barplot showing the CV for each gene identifier over the region of interest.
#' @author Benjamin Mayne
#' @export

plotcvExpr <- function(cvExpr){
  
  ### The output from cvExpr contains the CV for each gene 
  ### Which will be plot in a barplot
  cvPlot <- barplot(height = cvExpr[[1]], ylab = "CV", las = 2)
  
  ### The quantiles will be highlighted in four different colours
  abline(h = cvExpr[[2]][1:4], col = c("red", "blue", "green", "gold"), lwd = 4)
  legend("topright", cex = 1, col = c("red", "blue", "green", "gold"),
         fill =  c("red", "blue", "green", "gold"),
         legend = c("0%", "25%", "50%", "75%"))
}  
  
  
  
