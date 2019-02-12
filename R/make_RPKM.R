#' @title make_RPKM Make RPKM.
#' @description make_RPKM populates RPKM slot of modified countData (countDat)
#' S4 object.
#' @usage make_RPKM(cd_object)
#' @param cd_object A countDat object.
#' @return RPKM populated object
#'
#' @examples

#' data(hmel.data.doser)
#' reps <- c("Male", "Male", "Male", "Female", "Female", "Female")
#' annotxn <- data.frame("Chromosome" = factor(hmel.dat$chromosome,
#' levels = 1:21))
#' hmtr <- hmel.dat$trxLength
#' hm<-new("countDat",data=hmel.dat$readcounts,replicates=reps,seglens=hmtr,
#' annotation=annotxn)
#' libsizes(hm) <- getLibsizes2(hm, estimationType = "total")
#' rpkm(hm) <- make_RPKM(hm)

#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

make_RPKM <- function(cd_object) {
return(
sweep(sweep(cd_object@data, 2, as.vector(libsizes(cd_object)/1000000),
FUN = '/'), 1, cd_object@rowObservables$seglens/1000, FUN = '/')
)
}
