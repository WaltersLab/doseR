#' @title quantFilter Function to filter expression data of an
#' se object.
#' @description This function filters the expression of the supplied
#' se object; quantFilter is a filtering function used to remove
#' rows (genes) of various expression data.
#' @usage quantFilter (se, lo.bound=.25, hi.bound=.75, MEDIAN = FALSE,
#' na.rm = TRUE)
#' @param se An se object.
#' @param lo.bound The lower cutoff, expressed as a percentage.
#' @param hi.bound The upper cutoff, expressed as a percentage.
#' @param MEDIAN Boolean, Calculate RowMeans or RowMedians.
#' @param na.rm Boolean, NA removal.
#' @details This function filters the expression of the supplied se
#' object, based on a selected percentage cutoff.
#' @return Returns a filtered se object.
#' @examples
#' data(hmel.se)
#' f_se <- quantFilter(se, lo.bound=0.5)
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

quantFilter <- function(se, lo.bound=.25, hi.bound=.75, MEDIAN =
FALSE, na.rm = TRUE) {

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
rpkm <- log2(matrixStats::rowMedians( assays(se)$rpkm, na.rm = na.rm))
} else {
rpkm <- log2(rowMeans( assays(se)$rpkm, na.rm = na.rm ))
}

interq <- quantile(x = rpkm[is.finite(rpkm)], probs = c(lo.bound, hi.bound) )

outliers <- which(rpkm < interq[1] | rpkm > interq[2] )

percent <- function(x, digits = 2, format = "f", ...) { paste0(formatC(
100 * x, format = format, digits = digits, ...), "%") }
pct <- percent(  length(outliers)  / length(rpkm)   )
message( "Filtering removed ", length(outliers), " (", pct, ") of ",
length(rpkm), " total loci." )

return(se[-outliers,]) # remove all rows with "outlier" status

}
