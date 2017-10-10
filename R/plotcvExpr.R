
plotcvExpr <- function(cvExpr){
  
  library(ggplot2)
  
  datPlot <- data.frame(Gene = names(cvExpr[[1]]), CV = cvExpr[[1]])
  
  p <-ggplot(datPlot, aes(Gene, CV))
  p + geom_bar(stat = "identity")  +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      geom_hline(yintercept = cvExpr[[2]][1:4], color = c("red", "blue", "green", "gold"), show.legend = TRUE)
  
}
