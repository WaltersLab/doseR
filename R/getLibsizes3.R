#' @title getLibsizes3 Get Lib Sizes.
#' @description getLibsizes3 method for "SE" class object, derived
#' from getLibsizes2
#' @usage getLibsizes3(se, subset = NULL,
#' estimationType = c("quantile", "total",
#' "edgeR"),quantile = 0.75, ...)
#' @param se An s.e. object.
#' @param subset Value
#' @param estimationType e.g. quantile, total, edgeR.
#' @param quantile A quantile, expressed as e.g. 0.75.
#' @param ... Passthrough arguments.
#' @return Libsize value
#' @examples
#' data(hmel.se)
#' getLibsizes3(se)
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

getLibsizes3 <- function(se, subset = NULL,
estimationType = c("quantile", "total",   "edgeR"),
quantile = 0.75, ...) {
data<-assays(se)$counts
replicates <- colData(se)$Treatment
if(missing(subset)) subset <- NULL
if(is.null(subset)) subset <- seq_len(nrow(data))
estimationType = match.arg(estimationType)
if(is.na(estimationType)) stop("'estimationType' not known")
estLibs <- function(data, replicates)
{
libsizes <- switch(estimationType,
total = colSums(data[subset,,drop = FALSE], na.rm = TRUE),
quantile = apply(data[subset,, drop = FALSE], 2, function(z) {
x <- z[z > 0]
sum(x[x <= quantile(x, quantile, na.rm = TRUE)], na.rm = TRUE)
}),
edgeR = {
if(!("edgeR" %in% loadedNamespaces()))
requireNamespace("edgeR", quietly = TRUE)
d <- edgeR::DGEList(counts = data[subset,, drop = FALSE],
lib.size = colSums(data, na.rm = TRUE))
d <- edgeR::calcNormFactors(d, ...)
d$samples$norm.factors * d$samples$lib.size
})
names(libsizes) <- colnames(data)
libsizes
}
if(length(dim(data)) == 2) estLibsizes <- estLibs(data, replicates)
if(length(dim(data)) == 3) {
combData <- do.call("cbind", lapply(seq_len(dim(data)[3]),
function(kk) data[,,kk]))
combReps <- paste(as.character(rep(replicates, dim(data)[3])),
rep(c("a", "b"), each = ncol(data)), sep = "")
estLibsizes <- estLibs(combData, combReps)
estLibsizes <- do.call("cbind",
split(estLibsizes, cut(seq_len(length(estLibsizes)), breaks =
dim(data)[3], labels = FALSE)))
}
if(!missing(se))
if(inherits(se, what = "pairedData")) return(list(
estLibsizes[seq_len(ncol(se))], estLibsizes[seq_len(ncol(se)) + ncol(se)]))
if(length(dim(data)) > 2) estLibsizes <- array(estLibsizes,
dim = dim(assays(se)$counts)[-1])
return(estLibsizes)
}
