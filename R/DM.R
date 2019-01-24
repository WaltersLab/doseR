require(lme4)

#' @title cD.DM Function to convert count dat object to LME4 input.
#' @description This function generates LME4 input.
#' @usage cD.DM(cD, weightByLL = TRUE, includeGroupStructure = FALSE)
#' @param cD A countDat object containing FPKM values and at least one annotation column.
#' @param weightByLL Logical, weigh by log likelihood score.
#' @param includeGroupStructure Logical, include group structure.
#' @details This function converts count dat object.
#' @return Returns LME4 input.
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

cD.DM <- function (cD, weightByLL = TRUE, includeGroupStructure = FALSE)
{
  cdl <- cD
  adjustScaling <- function(CD) {
    libsizes(CD) <- libsizes(CD)/sum(libsizes(CD)) * sum(getLibsizes2(CD,
                                                                     estimationType = "total"))
    CD
  }
  cdl <- adjustScaling(cdl)
  dat <- do.call("rbind.data.frame", lapply(1:ncol(cdl@data), function(ii) {

    dat <- data.frame(rowID = 1:nrow(cdl@data), data = round(cdl@data[,
                                                                 ii]), sample = factor(colnames(cdl@data)[ii], levels = colnames(cdl@data)),
                      replicate = cdl@replicates[ii])
    if (length(cdl@rowObservables) > 0)
      dat <- cbind(dat, do.call("cbind", cdl@rowObservables))
    if (length(cdl@cellObservables) > 0)
      dat <- cbind(dat, do.call("cbind", lapply(cdl@cellObservables,
                                                function(x) x[, ii])))
    if (includeGroupStructure & length(groups(cdl)) > 0)
      for (gg in 1:length(groups(cdl))) dat[, paste("Group:",
                                                    names(groups(cdl))[gg], sep = "")] <- groups(cdl)[[gg]][ii]
      dat <- cbind(dat, cdl@annotation)

      dat
  }))
  if (length(cdl@sampleObservables) > 0)
    for (ss in 1:length(cdl@sampleObservables)) dat[, names(cdl@sampleObservables)[ss]] <- rep((cdl@sampleObservables[[ss]]),
                                                                                               each = nrow(cdl@data))
  dat <- dat[!is.na(dat[, 1]), , drop = FALSE]
  message("Checking for duplicate columns...", appendLF = FALSE)
  dupCols <- c(FALSE, sapply(2:ncol(dat), function(ii) {

    message(".", appendLF = FALSE)
    any(sapply(1:(ii - 1), function(jj) sum(diag(table(dat[,
                                                           jj], dat[, ii]))) == nrow(dat)))
  }))
  message("done!")
  if (any(dupCols)) {
    message(paste("Column: ", colnames(dat)[dupCols], "is a duplicate; removing this column \n"))
    dat <- dat[, !dupCols]
  }
  message("Column names; these and only these may be used to define formulae for the glSeq function.")
  message(paste(colnames(dat), collapse = "\n"))
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
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

glSeq <- function (dm, model, ...)
  lme4::glmer(formula(paste("data ~ ", model, "+ offset(log(libsizes * seglens / 1e9)) + (replicate|rowID)")),
        family = poisson, data = dm, ...)
