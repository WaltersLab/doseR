#' @title plotMA.se Function to make MA plot.
#' @description This function generates MA plot.
#' @usage plotMA.se(se, samplesA, samplesB, scale = NULL,
#' xlab = 'A', ylab = 'M', ...)
#' @param se An se object.
#' @param samplesA Either a character vector, identifying sample
#' set A by either replicate name or sample name, or a numerical
#' vector giving the columns of data in the 'se' object that
#' forms sample set A.
#' @param samplesB Either a character vector, identifying sample
#' set B by either replicate name or sample name, or a numerical
#' vector giving the columns of data in the 'se' object that
#' forms sample set B.
#' @param scale If given, defines the scale on which the log-ratios
#' will be plotted.
#' @param xlab Label for the X-axis. Defaults to 'A'.
#' @param ylab Label for the Y-axis. Defaults to 'M'.
#' @param ... Any other parameters to be passed to the plot
#' function.
#' @details This function makes MA plot from se object.
#' @return Returns MA plot.
#' @examples
#' data(hmel.se)
#' plotMA.se(se, samplesA = 'Male', samplesB = 'Female')
#' @author AJ Vaestermark, JR Walters.
#' @references The 'doseR' package, 2018 (in press).

plotMA.se <- function(se, samplesA, samplesB, scale = NULL, xlab = "A",
ylab = "M", ...) {
    if (is.character(samplesA)) {
        Asamps <- which(as.character(colData(se)$Treatment) %in% samplesA)
        if (!all(samplesA %in% colData(se)$Treatment))
            Asamps <- c(Asamps, which(colnames(assays(se)$rpkm) %in% samplesA[!(
            samplesA %in% as.character(colData(se)$Treatment))]))
        if (!all(samplesA %in% c(colnames(assays(se)$rpkm),
        as.character(colData(se)$Treatment))))
            warning("Some members of 'samplesA' were not found!")
        samplesA <- Asamps
    }
    if (length(samplesA) == 0) stop("Can't find any data for sample set A!")
    if (is.character(samplesB)) {
        Bsamps <- which(as.character(colData(se)$Treatment) %in% samplesB)
        if (!all(samplesB %in% colData(se)$Treatment))
            Bsamps <- c(Bsamps, which(colnames(assays(se)$rpkm) %in%
            samplesB[!(samplesB %in% as.character(colData(se)$Treatment))]))
        if (!all(samplesB %in% c(colnames(assays(se)$rpkm),
        as.character(colData(se)$Treatment))))
            warning("Some members of 'samplesB' were not found!")
        samplesB <- Bsamps
    }
    if (length(samplesB) == 0)
        stop("Can't find any data for sample set B!")
    if (!inherits(se, what = "SummarizedExperiment"))
        stop("variable 'se' must be of or descend from class 'se'")
    Adata <- assays(se)$rpkm[, samplesA]; Bdata <- assays(se)$rpkm[, samplesB]
    Adata <- colSums(t(Adata))/length(samplesA)
    Bdata <- colSums(t(Bdata))/length(samplesB)
    Azeros <- which(Adata == 0);     Bzeros <- which(Bdata == 0)
    nonzeros <- which(Adata != 0 & Bdata != 0)
    infRatio <- ceiling(max(abs((log2(Adata) -
    log2(Bdata))[nonzeros]), na.rm = TRUE))
    if (!is.null(scale) && scale > infRatio)
        infRatio <- scale
    M <- log2(Adata) - log2(Bdata)
    M[Azeros] <- -infRatio - 2; M[Bzeros] <- infRatio + 2
    A <- (log2(Adata) + log2(Bdata))/2
    A[Azeros] <- log2(Bdata[Azeros]); A[Bzeros] <- log2(Adata[Bzeros])
    plot(y = M, x = A, ylim = c(-infRatio - 3, infRatio + 3),
    axes = FALSE, xlab = xlab, ylab = ylab, ...); axis(side = 1)
    maxis <- pretty((-infRatio + 1):(infRatio - 1), min.n = 3,
                    n = length(axTicks(side = 2)))
    maxis <- maxis[maxis < infRatio & maxis > -infRatio]
    axis(side = 2, at = c(-infRatio - 1, maxis, infRatio + 1),
    labels = c(-Inf, maxis, Inf))
    abline(h = c(-1, 1) * (1 + infRatio), col = "orange", lty = 3)
}
