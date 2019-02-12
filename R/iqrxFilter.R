#' @title iqrxFilter Function to filter expression data within a countDat
#' object.
#' @description This function filters the expression of the supplied countDat
#' object. iqrxFilter is a filtering function used to remove rows (genes) of
#' various expression data.
#' @usage iqrxFilter(cD, iqr_multi = 1.5, MEDIAN = FALSE, na.rm = TRUE)
#' @param cD A countDat object.
#' @param iqr_multi Numeric multiplier; removes any outliers that are iqr_multi
#'  times the mid-50 percentile distance greater or less than the 25th and 75th
#'  percentiles, by default
#' @param MEDIAN Boolean, Calculate RowMeans or RowMedians.
#' @param na.rm Boolean, NA removal.
#' @details This function filters the expression of the supplied cD object,
#' based on a selected percentage cutoff and selected interquartile range
#' multiplier. The function iqrxFilter will: 1) log-base two transform all RPKM
#' values (obligatory); (2) remove any outliers that were 1.5 times the mid-50
#' percentile distance greater or less than the 75th and 25th percentiles
#' (by default), respectively; and (3) uses mean values and instead of median
#' values (by default).
#' @return Returns a filtered countDat object.
#'
#' @examples
#' data(hmel.data.doser)
#' reps <- c("Male", "Male", "Male", "Female", "Female", "Female")
#' annotxn <- data.frame("Chromosome" = factor(hmel.dat$chromosome,
#' levels = 1:21))
#' hm.tr<-hmel.dat$trxLength
#' hm<-new("countDat",data=hmel.dat$readcounts,seglens=hm.tr,
#' annotation=annotxn)
#' replicates(hm) <- reps
#' libsizes(hm) <- getLibsizes2(hm, estimationType = "total")
#' rpkm(hm) <- make_RPKM(hm)
#' f_hm <- iqrxFilter(hm)
#'
#' @author AJ Vaestermark, JR Walters.
#' @references Jue et al. BMC Genomics 2013 14:150

iqrxFilter <- function(cD, iqr_multi = 1.5, MEDIAN = FALSE, na.rm = TRUE) {

if(MEDIAN) {

rpkm <- log2(matrixStats::rowMedians( cD@RPKM, na.rm = na.rm))
} else {
rpkm <- log2(rowMeans( cD@RPKM, na.rm = na.rm ))
}

interq <- quantile(x = rpkm[is.finite(rpkm)], probs = c(.25, .75) )
# get 25 & 75 quartiles, excluding true 0 values
iqr <- interq[2] - interq[1] # get interquartile range
lo.bound <- interq[1] - (iqr * iqr_multi)
hi.bound <- interq[2] + (iqr * iqr_multi)

outliers <- which(rpkm < lo.bound | rpkm > hi.bound)

percent <- function(x, digits = 2, format = "f", ...) {
paste0(formatC(100 * x, format = format, digits = digits, ...), "%") }
pct <- percent(  length(outliers)  / length(rpkm)   )
message( "Filtering removed ", length(outliers), " (", pct, ") of ",
length(rpkm), " total loci." )

return(cD[-outliers,]) # remove all rows with "outlier" status

}
