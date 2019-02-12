#' @title generateStats Statistics function to summarize expression
#' data of a countDat object.
#' @description generateStats is a summary function used on various
#' expression data.
#' @usage generateStats(cD, groupings= NULL, mode_mean=TRUE, LOG2=TRUE)
#' @param cD A countDat object.
#' @param groupings A grouping (annotation column), e.g. groupings="something".
#' @param mode_mean Boolean, Calculate RowMeans or RowMedians.
#' @param LOG2 Boolean, Calculate LOG2.
#' @details This function completes summary statistics of the expression
#' of the supplied countDat object.
#' @return Returns an invisible list of summary statistics, kruskal test
#' and raw data of a countDat object.
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
#' generateStats(hm, groupings='Chromosome')
#'
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

generateStats <- function (cD, groupings= NULL, mode_mean=TRUE, LOG2=TRUE) {

if(is.null(groupings)) {
stop ('No groupings, e.g. groupings="something"...')
return (NULL)
}

MyGroups<-cD@annotation[[groupings]]

RM<- (
if(mode_mean) rowMeans(cD@RPKM) else matrixStats::rowMedians(cD@RPKM)
)

if(LOG2) { RM<-log2(RM) }
RM[is.infinite(RM)]<-NA

outsize <- 6
if(length(RM[is.na(RM)])>0) { outsize<-7 }

tmp <-  split(RM, MyGroups)
val.summary <-  vapply(tmp, summary, double(outsize) )
val.k <- kruskal.test(tmp)

names(tmp) <- levels(as.factor(MyGroups))
message("View output: outlist$kruskal, outlist$summary")
invisible( list("summary" = val.summary, "kruskal" = val.k, "data" = tmp) )

}# generateStats
