#' @title simpleFilter Function to filter expression data of an se object.
#' @description simpleFilter is a filtering function used to remove rows
#' (genes) of various expression data.
#' @usage simpleFilter(se, mean_cutoff=NULL, min_cutoff=NULL, median_cutoff=
#' NULL, counts=TRUE)
#' @param se An se object.
#' @param mean_cutoff The lower cutoff, using mean.
#' @param min_cutoff The cutoff, expressed as a minimum acceptable value.
#' @param median_cutoff The cutoff, using a median value.
#' @param counts Boolean, use raw counts data (default) or RPKM.
#' @details This function filters the expression of the supplied se obj.
#' @return Returns a filtered se object.
#' @examples
#' data(hmel.se) ; f_se <- simpleFilter(se, mean_cutoff=0.5)
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).
simpleFilter <- function(se, mean_cutoff=NULL, min_cutoff=NULL,
median_cutoff=NULL, counts=TRUE)
{ if(length(assays(se)$rpkm) == 0) {
return (NULL) }
ori_len<- nrow(assays(se)$counts)
if(counts) { if(!is.null(mean_cutoff)) {
se <- se[  apply(assays(se)$counts , 1, function(x) mean(x) ) >= mean_cutoff,]}
if(!is.null(min_cutoff)) {
se <- se[   apply(assays(se)$counts , 1, function(x) min(x) ) >= min_cutoff,]}
if(!is.null(median_cutoff)) {
se <- se[apply(assays(se)$counts , 1, function(x) median(x) )>=median_cutoff,]}
} else {if(!is.null(mean_cutoff)) {
se <- se[   apply(assays(se)$rpkm , 1, function(x) mean(x) ) >= mean_cutoff,]}
if(!is.null(min_cutoff)) {
se <- se[   apply(assays(se)$rpkm , 1, function(x) min(x) ) >= min_cutoff,]}
if(!is.null(median_cutoff)) {
se <- se[ apply(assays(se)$rpkm , 1, function(x) median(x) )>= median_cutoff,]
}} ; percent <- function(x, digits = 2, format = "f", ...) { paste0(
formatC(100 * x, format = format, digits = digits, ...), "%") }
pct <- percent( 1 - ( nrow(assays(se)$counts) / ori_len )  )
message( "Filtering removed ",ori_len-nrow(assays(se)$counts)," (",pct,")." )
invisible(se) }
