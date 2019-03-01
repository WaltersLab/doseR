#' @title plotRatioDensity Function to plot density of ratios between two
#' treatments (using groupings from an annotation column) within an se
#' object.
#' @description This function plots the expression of the supplied se
#' object, using ratios between a pair of selected treatments.
#' @usage plotRatioDensity(se, groupings= NULL, treatment1=NULL,
#' treatment2=NULL, mode_mean=TRUE, LOG2=TRUE,...)
#' @param se An se object.
#' @param groupings A grouping (annotation); groupings="annotation.ZA"
#' @param treatment1 Symbol, treatment 1.
#' @param treatment2 Symbol, treatment 2.
#' @param mode_mean Boolean, Calculate RowMeans or RowMedians.
#' @param LOG2 Boolean, Calculate LOG2.
#' @param ... Passthrough arguments to boxplot (additional arguments affecting
#' the summary produced).
#' @details This function plots expression of the supplied se object
#' using ratios of treatment1/treatment2.
#' @return Returns an invisible data frame containing the x-values and
#' corresponding density for each applicable annotation column entry.
#' @examples
#' data(hmel.se)
#' plotRatioDensity(se, groupings='annotation.ZA', treatment1 = 'Male',
#' treatment2 = 'Female',lty=1,type="l")
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

plotRatioDensity <- function (se, groupings= NULL, treatment1=NULL,
treatment2=NULL, mode_mean=TRUE, LOG2=TRUE,...) {
MyGroups<-rowData(se)[[groupings]]

if(is.null(groupings)) {
stop ('No groupings, e.g. groupings="something"...')
return (NULL)
}

if(length(assays(se)$rpkm) == 0) {
stop ('No RPKM data saved in se object... cancelling...')
return (NULL)
}

if(is.null(treatment1) | is.null(treatment2) ) {
stop ('Indicate treatments, such as  treatment1="A", treatment2="B"')
return (NULL)
}

if ( is.factor(MyGroups) ) { MyGroups <- droplevels(MyGroups) }

RM <- (
if(mode_mean) rowMeans(assays(se)$rpkm[,colData(se)$Treatment==treatment1])
else matrixStats::rowMedians(
assays(se)$rpkm[,colData(se)$Treatment==treatment1])
) / (
if(mode_mean) rowMeans(assays(se)$rpkm[,colData(se)$Treatment==treatment2])
else matrixStats::rowMedians(
assays(se)$rpkm[,colData(se)$Treatment==treatment2])
)

if (LOG2) { RM <- log2(RM) }
RM[is.infinite(RM)] <- NA
val <- split(RM, MyGroups)
xrge <- range(unlist(val), na.rm = TRUE, finite=TRUE)
myX <- density(na.omit(val[[1]]), from = xrge[1], to = xrge[2])$x
myDens <- vapply(val, FUN = function(x) { density(na.omit(x),
from = xrge[1], to = xrge[2])$y }, numeric(512) )
cat("Current order of levels: ", paste(unique(MyGroups), sep="\t") )
matplot(x = myX, y = myDens, ...)
invisible(data.frame(myX,myDens))

}# plotRatioDensity
