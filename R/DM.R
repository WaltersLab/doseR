#' @title se.DM Function to convert se object to LME4 input.
#' @description This function generates LME4 input.
#' @usage se.DM(se, weightByLL = TRUE)
#' @param se An se object containing FPKM values and at least
#' one annotation column.
#' @param weightByLL Logical, weigh by log likelihood score.
#' @details This function converts se object.
#' @return Returns LME4 input.
#' @examples
#' data(hmel.se)
#' f_se <- quantFilter(se, lo.bound = 0.4, hi.bound = 0.5)
#' dm <- se.DM(f_se)
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).
se.DM <- function (se, weightByLL = TRUE){
sel <- se
adjustScaling <- function(SE) {
colData(SE)$Libsizes <- colData(SE)$Libsizes/sum(colData(SE)$Libsizes) *
sum(getLibsizes3(SE,estimationType = "total"))
SE
}
sel <- adjustScaling(sel)
dat <- do.call("rbind.data.frame", lapply(seq_len(ncol(assays(sel)$counts)),
function(ii) {
dat <- data.frame(rowID = seq_len(nrow(assays(sel)$counts)), data = round(
assays(sel)$counts[,
ii]), sample = factor(colnames(assays(sel)$counts)[ii], levels = colnames(
assays(sel)$counts)),
replicate = colData(sel)$Treatment[ii])
if (length(list(rowData(sel)$seglens)) > 0)
dat <- cbind(dat, do.call("cbind", list(rowData(sel)$seglens)))
if (length(list()) > 0)
dat <- cbind(dat, do.call("cbind", lapply(list(),
function(x) x[, ii])))
dat <- cbind(dat, as.data.frame(rowData(sel)@listData))
dat
}))
for (ss in seq_len(length(list(sel@colData$Libsizes)))) dat[,
names(list(sel@colData$Libsizes))[ss]] <-
rep((list(sel@colData$Libsizes)[[ss]]),
each = nrow(assays(sel)$counts))
dat <- dat[!is.na(dat[, 1]), , drop = FALSE]
names(dat)[5]<-"seglens"

rep_ch<-NULL
for (ch in colData(sel)$Libsizes) {
rep_ch <- c( rep_ch, rep( ch, nrow(se) )  ) }
dat$libsizes <- rep_ch

invisible(dat) }

###########################

#' @title glSeq LME4 wrapper.
#' @description This function is an LME4 wrapper for dosage analysis.
#' @usage glSeq(dm, model, ...)
#' @param dm The dm data. Generally a reformatted se object from se.DM.
#' @param model The model. Expressed using standard LME4 syntax, see vignette.
#' @param ... passthrough arguments.
#' @details This function is an lme4 wrapper.
#' @return Returns LME4 output.
#' @examples
#' data(hmel.se)
#' f_se <- quantFilter(se, lo.bound = 0.4, hi.bound = 0.5)
#' dm <- se.DM(f_se)
#' glSeq(dm, "-1 + replicate")
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

glSeq <- function (dm, model, ...)
lme4::glmer(formula(paste(  "data ~ ", model, "+ offset(log( libsizes *
seglens / 1e9)) + (replicate|rowID)")),
family = poisson, data = dm, ...)
