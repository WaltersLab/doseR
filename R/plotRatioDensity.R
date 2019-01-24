#' @title plotRatioDensity Function to plot density of ratios between two treatments (using groupings from an annotation column) within a countDat object.
#' @description This function plots the expression of the supplied countDat object, using ratios between a pair of selected treatments.
#' @usage plotRatioDensity(cD, groupings, treatment1, treatment2, mode_mean, LOG2,...)
#' @param cD A countDat object.
#' @param groupings	A grouping (annotation column), e.g. groupings="something".
#' @param treatment1	Symbol, treatment 1.
#' @param treatment2	Symbol, treatment 2.
#' @param mode_mean 	Boolean, Calculate RowMeans or RowMedians.
#' @param LOG2	Boolean, Calculate LOG2.
#' @param ...	Passthrough arguments to boxplot (additional arguments affecting the summary produced).
#' @details This function plots expression of the supplied countDat object using ratios of treatment1/treatment2.
#' @return Returns an invisible data frame containing the x-values and corresponding density for each applicable annotation column entry.
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

plotRatioDensity <- function (cD, groupings= NULL, treatment1=NULL, treatment2=NULL, mode_mean=TRUE, LOG2=TRUE,...) {
  MyGroups<-cD@annotation[[groupings]]

  if(is.null(groupings)) {
    stop ('No groupings, e.g. groupings="something"...')
    return (NULL)
  }

  if(length(cD@RPKM) == 0) {
    stop ('No RPKM data saved in count data object... cancelling...')
    return (NULL)
  }

  if(is.null(treatment1) | is.null(treatment2) ) {
    stop ('Indicate treatments, such as  treatment1="A", treatment2="B"')
    return (NULL)
  }

  if ( is.factor(MyGroups) ) { MyGroups <- droplevels(MyGroups) }

  RM <- (
    if(mode_mean) rowMeans(cD@RPKM[,cD@replicates==treatment1]) else matrixStats::rowMedians(cD@RPKM[,cD@replicates==treatment1])
  ) / (
    if(mode_mean) rowMeans(cD@RPKM[,cD@replicates==treatment2]) else matrixStats::rowMedians(cD@RPKM[,cD@replicates==treatment2])
  )

  if (LOG2) { RM <- log2(RM) }
  RM[is.infinite(RM)] <- NA
  val <- split(RM, MyGroups)
  xrge <- range(unlist(val), na.rm = TRUE, finite=TRUE)
  myX <- density(na.omit(val[[1]]), from = xrge[1], to = xrge[2])$x
  myDens <- sapply(val, FUN = function(x) { density(na.omit(x), from = xrge[1], to = xrge[2])$y } )
  cat("Current order of levels: ", paste(unique(MyGroups), sep="\t") )
  matplot(x = myX, y = myDens, ...)
  invisible(data.frame(myX,myDens))

}# plotRatioDensity
