#' @title make_RPKM Make RPKM.
#' @description make_RPKM populates RPKM slot of se
#' S4 object.
#' @usage make_RPKM(se)
#' @param se An se object.
#' @return RPKM populated object
#'
#' @examples
#' data(hmel.se)

#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

#require(SummarizedExperiment)

make_RPKM <- function(se) {
return(
sweep(sweep(assays(se)$counts, 2, as.vector(colData(se)$Libsizes/1000000),
FUN = '/'), 1, metadata(se)$seglens/1000, FUN = '/')
)
}
