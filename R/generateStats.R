#' @title generateStats Statistics function to summarize expression
#' data of an s.e. object.
#' @description generateStats is a summary function used on various
#' expression data.
#' @usage generateStats(se, groupings= NULL, mode_mean=TRUE, LOG2=TRUE)
#' @param se An s.e. object.
#' @param groupings A grouping (annotation column), e.g. groupings="something".
#' @param mode_mean Boolean, Calculate RowMeans or RowMedians.
#' @param LOG2 Boolean, Calculate LOG2.
#' @details This function completes summary statistics of the expression
#' of the supplied s.e. object.
#' @return Returns an invisible list of summary statistics, kruskal test
#' and raw data of an s.e. object.
#' @examples
#' data(hmel.se)
#' generateStats(se, groupings="annotation.ZA")
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

generateStats <- function (se, groupings= NULL, mode_mean=TRUE, LOG2=TRUE) {

if(is.null(groupings)) {
stop ('No groupings, e.g. groupings="something"...')
return (NULL)
}

MyGroups<-rowData(se)[[groupings]]

RM<- (
if(mode_mean) rowMeans(assays(se)$rpkm)
else matrixStats::rowMedians(assays(se)$rpkm)
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
