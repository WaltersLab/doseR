#' @title simpleFilter Function to filter expression data of a countDat object.
#' @description simpleFilter is a filtering function used to remove rows (genes) of various expression data.
#' @usage simpleFilter(cD, mean_cutoff=NULL, min_cutoff=NULL, median_cutoff=NULL, counts=TRUE)
#' @param cD A countDat object.
#' @param mean_cutoff	The lower cutoff, using mean.
#' @param min_cutoff The cutoff, expressed as a minimum acceptable value.
#' @param median_cutoff	The cutoff, using a median value.
#' @param counts Boolean, use raw counts data (default) or RPKM.
#' @details This function filters the expression of the supplied countDat object, based on a selected cutoff.
#' @return Returns a filtered countDat object.
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

simpleFilter <- function(cD, mean_cutoff=NULL, min_cutoff=NULL, median_cutoff=NULL, counts=TRUE)
{

  if(length(cD@RPKM) == 0) {
    stop ('No RPKM data saved in count data object... cancelling...')
    return (NULL)
  }

  ori_len<- nrow(cD@data)

  if(counts) {
    if(!is.null(mean_cutoff)) { cD <- cD[   apply(cD@data , 1, function(x) mean(x) ) >= mean_cutoff,] }
    if(!is.null(min_cutoff)) { cD <- cD[   apply(cD@data , 1, function(x) min(x) ) >= min_cutoff,] }
    if(!is.null(median_cutoff)) { cD <- cD[   apply(cD@data , 1, function(x) median(x) ) >= median_cutoff,] }
  } else {
    if(!is.null(mean_cutoff)) { cD <- cD[   apply(cD@RPKM , 1, function(x) mean(x) ) >= mean_cutoff,] }
    if(!is.null(min_cutoff)) { cD <- cD[   apply(cD@RPKM , 1, function(x) min(x) ) >= min_cutoff,] }
    if(!is.null(median_cutoff)) { cD <- cD[   apply(cD@RPKM , 1, function(x) median(x) ) >= median_cutoff,] }
  }

  if(is.null(mean_cutoff) && is.null(min_cutoff) && is.null(median_cutoff)) {
    message ('USE SYNTAX : simpleFilter(cd.head, min_cutoff = 1.744596e-09, counts=FALSE)')
    warning ('SAMPLE SYNTAX : simpleFilter(cd.head, min_cutoff = 20, mean_cutoff = 10, median_cutoff = 30 )')
    stop ('No filtering parameters specified... cancelling...')
  } else {
    percent <- function(x, digits = 2, format = "f", ...) { paste0(formatC(100 * x, format = format, digits = digits, ...), "%") }
    pct <- percent( 1 - ( nrow(cD@data) / ori_len )  )
    message( "Filtering removed ", ori_len-nrow(cD@data), " (", pct, ") of ", ori_len, " total loci." )
    invisible(cD)
  }
}
