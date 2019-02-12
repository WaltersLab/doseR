require(lme4)

#' @title cD.DM Function to convert count dat object to LME4 input.
#' @description This function generates LME4 input.
#' @usage cD.DM(cD, weightByLL = TRUE)
#' @param cD A countDat object containing FPKM values and at least
#' one annotation column.
#' @param weightByLL Logical, weigh by log likelihood score.
#' @details This function converts count dat object.
#' @return Returns LME4 input.
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
#' f_hm <- quantFilter(hm, lo.bound = 0.4, hi.bound = 0.5)
#' dm <- cD.DM(f_hm)

#'
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

cD.DM <- function (cD, weightByLL = TRUE)
{
cdl <- cD
adjustScaling <- function(CD) {
libsizes(CD) <- libsizes(CD)/sum(libsizes(CD)) * sum(getLibsizes2(CD,
estimationType = "total"))
CD
}
cdl <- adjustScaling(cdl)
dat <- do.call("rbind.data.frame", lapply(seq_len(ncol(cdl@data)),
function(ii) {

dat <- data.frame(rowID = seq_len(nrow(cdl@data)), data = round(cdl@data[,
ii]), sample = factor(colnames(cdl@data)[ii], levels = colnames(cdl@data)),
replicate = cdl@replicates[ii])
if (length(cdl@rowObservables) > 0)
dat <- cbind(dat, do.call("cbind", cdl@rowObservables))
if (length(cdl@cellObservables) > 0)
dat <- cbind(dat, do.call("cbind", lapply(cdl@cellObservables,
function(x) x[, ii])))

dat <- cbind(dat, cdl@annotation)

dat
}))
if (length(cdl@sampleObservables) > 0)
for (ss in seq_len(length(cdl@sampleObservables))) dat[,
names(cdl@sampleObservables)[ss]] <- rep((cdl@sampleObservables[[ss]]),
each = nrow(cdl@data))
dat <- dat[!is.na(dat[, 1]), , drop = FALSE]
invisible(dat)
}

###########################

#' @title glSeq LME4 wrapper.
#' @description This function is an LME4 wrapper for dosage analysis.
#' @usage glSeq(dm, model, ...)
#' @param dm The dm data. Generally a reformatted countDat object from cD.DM.
#' @param model The model. Expressed using standard LME4 syntax, see vignette.
#' @param ... passthrough arguments.
#' @details This function is an lme4 wrapper.
#' @return Returns LME4 output.

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
#' f_hm <- quantFilter(hm, lo.bound = 0.4, hi.bound = 0.5)
#' dm <- cD.DM(f_hm)
#' glSeq(dm, "-1 + replicate")

#'
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

glSeq <- function (dm, model, ...)
lme4::glmer(formula(paste("data ~ ", model, "+ offset(log(libsizes *
seglens / 1e9)) + (replicate|rowID)")),
family = poisson, data = dm, ...)
