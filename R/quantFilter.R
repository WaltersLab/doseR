#' @title quantFilter	Function to filter expression data of a countDat object.
#' @description This function filters the expression of the supplied countDat object; quantFilter is a filtering function used to remove rows (genes) of various expression data.
#' @usage quantFilter (cD, lo.bound=.25, hi.bound=.75, MEDIAN = FALSE, na.rm = TRUE)
#' @param cD A countDat object.
#' @param lo.bound	The lower cutoff, expressed as a percentage.
#' @param hi.bound The upper cutoff, expressed as a percentage.
#' @param MEDIAN Boolean, Calculate RowMeans or RowMedians.
#' @param na.rm Boolean, NA removal.
#' @details This function filters the expression of the supplied countDat object, based on a selected percentage cutoff.
#' @return Returns a filtered countDat object.
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

quantFilter <- function(cD, lo.bound=.25, hi.bound=.75, MEDIAN = FALSE, na.rm = TRUE) {

  if(lo.bound > hi.bound) {
    stop('lo.bound greater than hi.bound... cancelling...')
  }

  if(lo.bound < 0 | lo.bound > 1) {
    stop('lo.bound outside range... cancelling...')
  }

  if(hi.bound < 0 | hi.bound > 1) {
    stop('hi.bound outside range... cancelling...')
  }

  if(MEDIAN) {
    rpkm <- log2(matrixStats::rowMedians( cD@RPKM, na.rm = na.rm))
  } else {
    rpkm <- log2(rowMeans( cD@RPKM, na.rm = na.rm ))
  }

  interq <- quantile(x = rpkm[is.finite(rpkm)], probs = c(lo.bound, hi.bound) )

  outliers <- which(rpkm < interq[1] | rpkm > interq[2] )

  percent <- function(x, digits = 2, format = "f", ...) { paste0(formatC(100 * x, format = format, digits = digits, ...), "%") }
  pct <- percent(  length(outliers)  / length(rpkm)   )
  message( "Filtering removed ", length(outliers), " (", pct, ") of ", length(rpkm), " total loci." )

  return(cD[-outliers,]) # remove all rows with "outlier" status

}
