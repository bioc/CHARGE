### Sliding windows test to find additions or deletions ###

library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(parallel)
library(plyr)
library(matrixStats)

load("/home/SAHMRI.INTERNAL/benjamin.mayne/Bioconductor_Package/SRA/Counts/GSE55504_se.Rdata")

bins <- tileGenome(seqinfo(EnsDb.Hsapiens.v86), tilewidth=100000, cut.last.tile.in.chrom = T)


    
bimodalBin <- function(bin, se, threshold = NULL){
  
  ### Firstly subset se for the bin ###
  datExpr <- suppressWarnings(assay(subsetByOverlaps(se, bin)))
  
  ### If threshold is not set to NULL remove genes with low variance 
  if(!is.null(threshold)){
    
    ### Simialr to cvExpr calculate the standard deviations and means of each gene within the bin
    ### and determine the expression variation
    geneSds <- rowSds(datExpr)
    
    #### Calucalte the mean of each gene
    geneMeans <- rowMeans(datExpr)
    
    #### Calucalte CV = 100*sd/mean
    geneCV <- 100*geneSds/geneMeans
    
    ### Determine the genes quantiles 
    CVQuantiles <- quantile(geneCV)
    
    ## Subset the genes with a variation threshold greater than the input value from datExpr
    genes <- names(which(geneCV > CVQuantiles[threshold]))
    datExpr <- data.frame(datExpr)[genes,]
  }
  
  #### Calulate the mean z score for each sample using the genes within the bin
  datZscores <- scale(x = datExpr)
  ZscoreMeans <- rowMeans(x = datZscores) 

  ### If statement to end the function if there are not enough genes to compute the Z scores
  if(nrow(datExpr) == 0 || any(is.nan(ZscoreMeans))){
    
    datResult <- data.frame(bin)
    datResult$Bimodality.Amplitude <- NA
    datResult$Bimodality.Coefficient <- NA
    datResult$Bimodality.Ratio <- NA
    datResult$Dip.Statistic <- NA
    return(datResult)
    
  }
  
  ### Calcualte the biomodal amplitude, coeffeicent, ratio and dip test stat and return them with the bin
  bimod_amp <- as.numeric(suppressWarnings(bimodality_amplitude(ZscoreMeans, fig = FALSE)))
  bimod_coef <- as.numeric(suppressWarnings(bimodality_coefficient(t(datExpr))))
  bimod_ratio <- as.numeric(suppressWarnings(bimodality_ratio(x = ZscoreMeans, fig = FALSE)))
  dipResult <- as.numeric(dip.test(x = ZscoreMeans)["statistic"])
  
  ### If there were not enough genes within datExpr the bimodal test may not work
  ### Use the if else statement to return NA if the bimodal test failed
  datResult <- data.frame(bin)
  datResult$Bimodality.Amplitude <- ifelse(test = length(bimod_amp) == 0, yes = NA, no = bimod_amp)
  datResult$Bimodality.Coefficient <- ifelse(test = length(bimod_coef) == 0, yes = NA, no = bimod_coef)
  datResult$Bimodality.Ratio <- ifelse(test = length(bimod_ratio) == 0, yes = NA, no = bimod_ratio)
  datResult$Dip.Statistic <- ifelse(test = length(dipResult) == 0, yes = NA, no = dipResult)
  return(datResult)
  
  
}

dat <- mclapply(X = bins, se = se, threshold = 4, FUN=bimodalBin, mc.cores=8)
dat <- ldply(dat)

dat <- dat[which(dat$Dip.Statistic > 0.08),]

bimodalBin(bin = bins[150], se = se, threshold = 4)










