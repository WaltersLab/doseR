#' @title test_diffs Statistics function to summarize expression
#' data of a countDat object, using ratios between selected treatments.
#' @description generateStats is a summary function used on various
#' expression data, using ratios between selected treatments.
#' @usage test_diffs(cD, groupings= NULL, treatment1=NULL, treatment2=NULL,
#' mode_mean=TRUE, LOG2=TRUE)
#' @param cD A countDat object.
#' @param groupings A grouping (annotation column), e.g. groupings="something".
#' @param treatment1 Symbol, treatment 1.
#' @param treatment2 Symbol, treatment 2.
#' @param mode_mean Boolean, Calculate RowMeans or RowMedians.
#' @param LOG2 Boolean, Calculate LOG2.
#' @details This function completes summary statistics of the expression of
#' the supplied countDat object.
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
#' test_diffs(hm, groupings='Chromosome',treatment1="Male",
#' treatment2="Female" )
#'
#' @return Returns an invisible list of summary statistics,
#' kruskal test and raw data of a countDat object, using ratios
#' between selected treatments.
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

test_diffs <- function (cD, groupings= NULL, treatment1=NULL,
treatment2=NULL, mode_mean=TRUE, LOG2=TRUE) {

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

MyGroups<-cD@annotation[[groupings]]

RM<- (
    if(mode_mean) rowMeans(cD@RPKM[,cD@replicates==treatment1]) else
    matrixStats::rowMedians(cD@RPKM[,cD@replicates==treatment1])
) / (
    if(mode_mean) rowMeans(cD@RPKM[,cD@replicates==treatment2]) else
    matrixStats::rowMedians(cD@RPKM[,cD@replicates==treatment2])
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

}# test_diffs
