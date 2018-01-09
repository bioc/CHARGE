#' plotcvExpr
#'
#' Plots the coefficient of variation or expression variation for each gene over a defined genomic region.
#'
#' @param cvExpr The output from cvExpr function.
#' @usage plotcvExpr(cvExpr)
#' @return Returns a barplot showing the CV for each gene identifier over the region of interest.
#' @import graphics
#' @author Benjamin Mayne
#' @examples
#' library(GenomicRanges)
#' data(datExprs)
#' chr21 <- GRanges("21:1-46709983")
#' cvExpr.out <- cvExpr(se = datExprs, region = chr21)
#' plotcvExpr(cvExpr = cvExpr.out)
#' @details
#' Generates a bar plot showing the coefficient of variation or expression variation for each gene on the Y axis.
#' The red, blue, green and gold horizontal lines show the 0%, 25%, 50%, 75% quartiles respectively.
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


