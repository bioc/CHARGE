#' exprFinder
#'
#' Performs a bimodality test at multiple defined bin sizes across the genome using a sliding window approach.
#'
#' @param se A SummarizedExperiment containing the nromalised gene expression data.
#' @param ranges A GRanges object containing the genomic regions to scan.
#' @param binWidth The length of each bin.
#' @param binStep The distance the bin will slide.
#' @param threshold Optional. The quantile threshold of expression variation of genes to be used at each bin. Default is NULL.
#' @param threads Total number of threads to be used. Default is 1.
#' @usage exprFinder(se, ranges, binWidth, binStep, threshold = NULL, threads = 1)
#' @return Returns a data frame containing the genomic locations of each bin and bimodality statistics.
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import modes
#' @importFrom IRanges subsetByOverlaps
#' @importFrom methods is
#' @importFrom matrixStats rowSds
#' @importFrom stats density complete.cases
#' @import parallel
#' @import plyr
#' @author Benjamin Mayne
#' @examples
#' library(GenomicRanges)
#' library(EnsDb.Hsapiens.v86)
#' data(datExprs)
#' chr21 <- GRanges("21:1-46709983")
#' chrLengths <- GRanges(seqinfo(EnsDb.Hsapiens.v86)[c("21", "22", "Y")])
#' exprFinder.out <- exprFinder(se = datExprs, ranges = chrLengths,
#' binWidth = 1e+9, binStep = 1e+9, threshold = "25%")
#' @details
#' Uses a sliding window approach to scan over a defined genomic region.
#' It automatically performs a bimodal test and calculates Hartigan's dip test statistic for unimodality and returns a data frame
#' listing each bin and the statistical likelihood of a duplication or deletion.
#' @export

exprFinder <- function(se, ranges, binWidth, binStep, threshold = NULL, threads = 1){

  ### Unit tests to see if the inputted data is in the correct format
  #### se must be a RangedSummarizedExperiment
  if(!is(se, "RangedSummarizedExperiment")){
    stop("se must be a RangedSummarizedExperiment")
  }

  #### ranges must be a GRanges
  if(!is(ranges, "GRanges")){
    stop("ranges must be a GRanges")
  }

  #### binWidth must be a numeric
  if(!is(binWidth, "numeric")){
    stop("binWidth must be a numeric value")
  }

  ### Threads must be a numeric value
  if(!is(threads, "numeric")){
    stop("threads must be a numeric value")
  }


  ### The function bimodalBin is used to take the genes within each bin and test for bimodality
  ### If there are no genes found to be within the bin or not enough to compute Z scores then the
  ### function will return NAs for the bimodality test
  ### The output is a data frame containing the bin genomic location and the output from the bimodality and dip test
  bimodalBin <- function(bin, se, threshold){

    ### Firstly subset se for the bin ###
    datExpr <- suppressWarnings(assay(subsetByOverlaps(se, bin)))

    ### If threshold is not set to NULL remove genes with low variance
    if(!is.null(threshold)){

      ### Similar to cvExpr calculate the standard deviations and means of each gene within the bin
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

    ### Get the gene count ###
    geneCount <- nrow(datExpr)

    ### Tranpose the data frame
    datExpr <- data.frame(t(datExpr))

    #### Calulate the mean z score for each sample using the genes within the bin
    datZscores <- scale(x = datExpr)
    ZscoreMeans <- rowMeans(x = datZscores)

    ### If statement to end the function if there are not enough genes to compute the Z scores
    if(nrow(datExpr) == 0 || any(is.nan(ZscoreMeans))){

      datResult <- data.frame(bin)
      datResult$Bimodality.Coefficient <- NA
      datResult$Bimodality.Ratio <- NA
      datResult$Dip.Statistic <- NA
      return(datResult)

    }

    ### Calcualte the biomodal amplitude, coeffeicent, ratio and dip test stat and return them with the bin
    bimod_coef <- as.numeric(suppressWarnings(bimodality_coefficient(t(datExpr))))
    bimod_ratio <- as.numeric(suppressWarnings(bimodality_ratio(x = ZscoreMeans, fig = FALSE)))
    dipResultStat <- as.numeric(dip.test(x = ZscoreMeans)["statistic"])
    dipResultPValue <- as.numeric(dip.test(x = ZscoreMeans)["p.value"])

    ### If there were not enough genes within datExpr the bimodal test may not work
    ### Use the if else statement to return NA if the bimodal test failed
    datResult <- data.frame(bin)
    datResult$Bimodality.Coefficient <- ifelse(test = length(bimod_coef) == 0, yes = NA, no = bimod_coef)
    datResult$Bimodality.Ratio <- ifelse(test = length(bimod_ratio) == 0, yes = NA, no = bimod_ratio)
    datResult$Dip.Statistic <- ifelse(test = length(dipResultStat) == 0, yes = NA, no = dipResultStat)
    datResult$Dip.pvalue <- ifelse(test = length(dipResultPValue) == 0, yes = NA, no = dipResultPValue)

    datResult$No.Genes <- geneCount
    return(datResult)

  }

  ### Divide the genome up into bins using a sliding window
  bins <- unlist(slidingWindows(x = ranges, width = binWidth, step = binStep))

  ### Use the bimodalBin function on every bin
  ### mclapply can be used to make use of multiple cores (if possible)
  bimodalBinOut <- mclapply(X = bins, se = se, threshold = threshold, FUN = bimodalBin, mc.cores = threads)

  ### unlist bimodalBinOut into a singel data frame and return it
  bimodalBinOut <- ldply(bimodalBinOut)
  ### Remove regions that did not compute a p value
  bimodalBinOut <- bimodalBinOut[complete.cases(bimodalBinOut$Dip.pvalue),]

  ### Reduce bimodalBinOut to contain the regions that had the same genes
  bimodalBinReduced <- unlist(reduce(split(GRanges(bimodalBinOut), elementMetadata(GRanges(bimodalBinOut))$Dip.pvalue)))

  bimodalBinReduced <- cbind(data.frame(bimodalBinReduced),
                          bimodalBinOut[match(names(bimodalBinReduced), bimodalBinOut$Dip.pvalue), 6:10])

  ### Reorder the data frame such that low p values are the top ###
  bimodalBinReduced <- bimodalBinReduced[order(bimodalBinReduced$Dip.pvalue, decreasing = FALSE), ]
  row.names(bimodalBinReduced) <- NULL

  return(bimodalBinReduced)

}
