require(edgeR)

#' An S4 class that stores count data.
#' @slot data Contains counts data.
#' @slot RPKM Contains RPKM data.
#' @slot annotation Contains annotation column data.
#' @slot groups Contains groups data.
#' @slot replicates Contains replicates data.
#' @slot rowObservables Contains rowObservables data, (including seglens).
#' @slot sampleObservables Contains sampleObservables data.
#' @slot orderings Contains orderings data.
#' @slot nullPosts Contains nullPosts data.
#' @slot cellObservables Contains cellObservables data.

setClass("countDat", representation(data = "array", RPKM = "array", replicates = "factor", groups = "list", rowObservables = "list", sampleObservables = "list", annotation = "data.frame" , orderings = "data.frame", nullPosts = "matrix" , cellObservables = "list" ))

#' libsizes method for testClass
#'
#' @docType methods
#' @rdname libsizes-methods
#' @param x Value
#' @param value Value

setGeneric("libsizes<-", function(x, value) standardGeneric("libsizes<-"))

#' libsizes method for testClass
#'
#' @docType methods
#' @rdname libsizes-methods
#' @keywords internal

setMethod("libsizes<-", signature = "countDat", function(x, value) {
  x@sampleObservables$libsizes <- value
  x
})

#' libsizes method for testClass
#'
#' @docType methods
#' @rdname libsizes-methods

setGeneric("libsizes", function(x) standardGeneric("libsizes"))

#' libsizes method for testClass
#'
#' @docType methods
#' @rdname libsizes-methods

setMethod("libsizes", signature = "countDat", function(x) {
  x@sampleObservables$libsizes
})

###############################

#' getLibsizes method for count Dat class, derived from getLibsizes1 method (countData class).
#'
#' @docType methods
#' @rdname getLibsizes-methods
#' @param cD A count Dat object.
#' @param subset Value
#' @param estimationType e.g. quantile, total, edgeR.
#' @param quantile A quantile, expressed as e.g. 0.75.
#' @param ... Passthrough arguments.

'getLibsizes2' <- function(cD, subset = NULL, estimationType = c("quantile", "total", "edgeR"), quantile = 0.75, ...)
{

  data<-cD@data 				# NEED THIS
  replicates <- cD@replicates		# NEED THIS

  if(missing(subset)) subset <- NULL
  if(is.null(subset)) subset <- 1:nrow(data)
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
                         d <- edgeR::DGEList(counts = data[subset,, drop = FALSE], group = replicates, lib.size = colSums(data, na.rm = TRUE))
                         d <- edgeR::calcNormFactors(d, ...)
                         d$samples$norm.factors * d$samples$lib.size
                       })
    names(libsizes) <- colnames(data)
    libsizes
  }

  if(length(dim(data)) == 2) estLibsizes <- estLibs(data, replicates)
  if(length(dim(data)) == 3) {
    combData <- do.call("cbind", lapply(1:dim(data)[3], function(kk) data[,,kk]))
    combReps <- paste(as.character(rep(replicates, dim(data)[3])), rep(c("a", "b"), each = ncol(data)), sep = "")
    estLibsizes <- estLibs(combData, combReps)
    estLibsizes <- do.call("cbind",
                           split(estLibsizes, cut(1:length(estLibsizes), breaks = dim(data)[3], labels = FALSE)))
  }

  if(!missing(cD))
    if(inherits(cD, what = "pairedData")) return(list(estLibsizes[1:ncol(cD)], estLibsizes[1:ncol(cD) + ncol(cD)]))

  if(length(dim(data)) > 2) estLibsizes <- array(estLibsizes, dim = dim(cD@data)[-1])

  return(estLibsizes)
}

###############################

#' subsetting method for countDat class
#'
#' @docType methods
#' @rdname extract-methods
#' @param x countDat object Value
#' @param i first dimension, subsetting parameter
#' @param j second dimension, subsetting parameter
#' @param ... Passthrough arguments.
#' @param drop Value, Logical.

setMethod("[", "countDat", function(x, i, j, ..., drop = FALSE) {
  if(missing(j)) {
    j <- 1:ncol(x@data)
  } else {
    if(is.logical(j)) j <- which(j)

    if(!all(1:ncol(x@data) %in% j))
    {
      replicates(x) <- as.character(x@replicates[j])

      if(length(x@groups) > 0)
      {
        newgroups <- list()
        newgroups <- lapply(x@groups, function(x) {
          x[j]
          rep(1:length(unique(x[j])), sapply(unique(x[j]), function(z) sum(x[j] == z)))[unlist(sapply(unique(x[j]), function(z) which(x[j] == z)))]
        })
        x@groups <- newgroups[!duplicated(newgroups) | duplicated(x@groups)]
      }

      if(length(x@orderings) > 0)
      {
        warning("Selection of samples (columns) will invalidate the values calculated in slot 'orderings', and so these will be discarded.")
        x@orderings <- data.frame()
      }

    }
  }

  if(missing(i))
    i <- 1:nrow(x@data)
  if(is.logical(i)) i <- which(i)

  if(nrow(x@data) > 0)
    x@data <- .sliceArray2(list(i, j), x@data)
  x@RPKM <- .sliceArray2(list(i, j), x@RPKM)

  x@annotation <- x@annotation[i,, drop = FALSE]
  if(nrow(x@orderings) > 0)
    x@orderings <- x@orderings[i,, drop = FALSE]
  if(length(x@nullPosts) > 0)
    x@nullPosts <- x@nullPosts[i,,drop = FALSE]

  x@rowObservables <- lapply(x@rowObservables, function(z) .sliceArray2(list(i),z, drop = FALSE))
  x@sampleObservables <- lapply(x@sampleObservables, function(z) .sliceArray2(list(j), z, drop = FALSE))
  x@cellObservables <- lapply(x@cellObservables, function(z) .sliceArray2(list(i,j), z, drop = FALSE))

  x
})

###############################

.sliceArray2 <- function(slices, array, drop = FALSE) {
  if((is.vector(array) & sum(!sapply(slices, is.null)) > 1) || (is.array(array) & length(slices) > length(dim(array)))) warning("dimensions of slice exceed dimensions of array")
  sarray <- abind::asub(array, slices, dims = 1:length(slices), drop = drop)
  sarray
}

#' replicates method for testClass
#'
#' @docType methods
#' @rdname replicates-methods
#' @param x Value
#' @param value Value

setGeneric("replicates<-", function(x, value) standardGeneric("replicates<-"))

#' replicates method for testClass
#'
#' @docType methods
#' @rdname replicates-methods
#' @keywords internal

setMethod("replicates<-", signature = "countDat", function(x, value) {
  x@replicates <- as.factor(value)
  x
})

######################

setMethod("show", "countDat", function(object) {

  cat(paste('An object x of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))

  cat('\nSlot "replicates"\n')
  cat(as.character(object@replicates))

  cat('\nSlot "groups":\n')
  print(object@groups)

  cat('\nSlot "data":\n')

  if(nrow(object@data) > 5)
  {
    print(.showData(.sliceArray2(list(1:5), object@data)), quote = FALSE)
    cat(paste(nrow(object) - 5), "more rows...\n")
  } else print(.showData(object@data))

  cat('\nSlot "RPKM":\n')

  if(nrow(object@RPKM) > 5)
  {
    print(.showData(.sliceArray2(list(1:5), object@RPKM)), quote = FALSE)
    cat(paste(nrow(object) - 5), "more rows...\n")
  } else print(.showData(object@RPKM))

  cat('\nSlot "annotation":\n')
  if(nrow(object@annotation) > 5 & ncol(object@annotation) > 0)
  {
    print(object@annotation[1:5,])
    cat(paste(nrow(object) - 5), "more rows...\n")
  } else print(object@annotation)

})

########################

.showData <- function(data)
{
  if(is.vector(data) || length(dim(data)) <= 2) return(data)
  dimsep <- c(":", "|")
  dimlen <- length(dim(data))
  if(length(dim(data)) > 4)
    dimsep <- c(dimsep, sapply(2:(dimlen - 2), function(x) paste(rep("|", x), collapse = "")))
  dimsep <- c("", dimsep)
  dimsep <- dimsep[1:(dimlen - 1)]

  dimsep <- rev(dimsep)
  dimdat <- data

  pasteDat <- function(x, dimnum) {
    if(length(dim(x)) > 2) {
      padat <- t(apply(x, 1, function(xx) paste(pasteDat(xx, dimnum = dimnum + 1), collapse = dimsep[dimnum])))
    } else {
      padat <- (apply(x, 1, function(z) paste(z, collapse = ":")))
    }
    return(padat)
  }
  pastemat <- t(apply(data, 1, pasteDat, dimnum = 1))
  pastemat
}

######################

setGeneric(".seglens<-", function(x, value) standardGeneric(".seglens<-"))
setMethod(".seglens<-", signature = "countDat", function(x, value) {
  if(!is.numeric(value)) stop("All members of seglens for a countData object must be numeric.")

  if(inherits(value, "numeric")) {
    if(length(value) != ncol(x)) stop("Length of seglens must be identical to the number of columns of the countData object.")
    value <- matrix(value, ncol = 1)
  } else if(is.array(value))
    if(any(dim(x@data)[-1] != dim(value))) stop("Dimension of seglens must be identical to the dimension of the countData object (after dropping the first dimension).")

  if(any(value <= 0)) stop("Library sizes less than or equal to zero make no sense to me!")
  x@rowObservables$seglens <- value
  x
})

setGeneric(".seglens", function(x) standardGeneric(".seglens"))
setMethod(".seglens", signature = "countDat", function(x) {
  if(".seglens" %in% names(x@rowObservables)) return(x@rowObservables$seglens)
  if(".seglens" %in% names(x@cellObservables)) return(x@cellObservables$seglens)
  return(matrix(rep(1, nrow(x)), ncol = 1))
})

#####################

setMethod("initialize", "countDat", function(.Object, ..., data, replicates, libsizes, seglens) {

  .Object <- callNextMethod(.Object, ...)

  if(!missing(data) && is.array(data)) .Object@data <- data
  if(!missing(data) && is.list(data)) .Object@data <- array(do.call("c", data), c(dim(data[[1]]), length(data)))
  if(missing(replicates)) replicates <- .Object@replicates
  .Object@replicates <- as.factor(replicates)

  if(length(dim(.Object@data)) == 1) .Object@data <- array(.Object@data, dim = c(dim(.Object@data), max(c(0, length(replicates), length(.Object@replicates)))))

  if(length(colnames(.Object@data)) == 0) colnames(.Object@data) <- make.unique(c(as.character(unique(.Object@replicates)), as.character(.Object@replicates)))[-(1:(length(unique(.Object@replicates))))]

  if(nrow(.Object@annotation) > 0 & nrow(.Object@annotation) != nrow(.Object@data))
    warning("Number of rows of '@annotation' slot not same as '@data' slot.")

  if(any(lapply(.Object@groups, length) != ncol(.Object@data)))
    stop("All vectors in '@groups' slot must equal number of columns of '@data' slot.")

  if(length(.Object@nullPosts) != 0) {
    if(nrow(.Object@nullPosts) != nrow(.Object@data) & nrow((.Object@nullPosts) != 0))
      stop("Number of rows in '@data' slot must equal number of rows of '@nullPosts' slot.")
  } else nullPosts <- matrix(ncol = 0, nrow = nrow(.Object@data))

  .Object@groups <- lapply(.Object@groups, as.factor)

  if(!missing(libsizes)) {
    if(is.array(libsizes) && (any(dim(libsizes) != dim(.Object@data)[-1])) || (is.vector(libsizes) & length(libsizes) != ncol(.Object@data)))
      stop("If provided, the 'libsizes' variable must be a vector of equal length to the columns of the `@data' array or an array of equal dimension to a row of the `@data' array")
    if(is.array(libsizes) && is.null(colnames(libsizes))) colnames(libsizes) <- colnames(.Object@data)
    if(is.vector(libsizes) && is.null(names(libsizes))) names(libsizes) <- colnames(.Object@data)
    .Object@sampleObservables$libsizes <- libsizes
  }

  if(!missing(seglens))
  {
    if(is.vector(seglens)) {
      if(length(seglens) != nrow(.Object@data)) stop("If 'seglens' specified, and is a vector, the length of this variable must equal the number of rows of '@data' slot.")
      .Object@rowObservables$seglens <- seglens
    }
    if(is.array(seglens)) {
      if(length(dim(.Object@data)) != length(dim(seglens)) || (any(dim(.Object@data) != dim(seglens)))) stop("If 'seglens' specified, and is an array, the dimensions of this variable must equal the dimensions of the '@data' slot.")
      .Object@cellObservables$seglens <- seglens
    }
  }

  if(length(.Object@rowObservables) > 0) {
    notRow <- sapply(.Object@rowObservables, length) != nrow(.Object@data)
    if(any(notRow)) stop(paste("The following '@rowObservables' elements have an incorrect length:", paste(names(notRow)[notRow], collapse = ",")))
  }
  if(length(.Object@sampleObservables) > 0) {
    notSample <- sapply(.Object@sampleObservables, function(x)
      (is.vector(x) && length(x) != ncol(.Object@data)) | (is.array(x) && ((length(dim(x)) != length(dim(.Object@data)) - 1) | any(dim(x) != dim(.Object@data)[-1]))))

    if(any(notSample)) stop(paste("The following '@sampleObservables' elements have an incorrect length:", paste(names(notSample)[notSample], collapse = ",")))
  }
  if(length(.Object@cellObservables) > 0) {
    notCell <- sapply(.Object@cellObservables, function(oco) any(dim(oco)[1:2] != dim(.Object@data)[1:2]))
    if(any(notCell)) stop(paste("The following '@cellObservables' elements have incorrect dimensions:", paste(names(notCell)[notCell], collapse = ",")))
  }

  if(length(replicates) != 0 && length(replicates) != ncol(.Object@data))
    stop("The length of the '@replicates' slot must equal number of columns of '@data' slot.")

  .Object
})

#' replicates method for testClass
#'
#' @docType methods
#' @rdname replicates-methods
#' @param x Value
#' @param value Value

setGeneric("groups<-", function(x, value) standardGeneric("groups<-"))

#' replicates method for testClass
#'
#' @docType methods
#' @rdname replicates-methods
#' @keywords internal

setMethod("groups<-", signature = "countDat", function(x, value) {
  if(any(sapply(value, length) != ncol(x))) stop(paste(sum(sapply(value, length) != ncol(x)), "vector(s) in the groups structure are the wrong length."))
  x@groups <- lapply(value, as.factor)
  x
})

#' replicates method for testClass
#'
#' @docType methods
#' @rdname replicates-methods
#' @param x Value
#' @param value Value

setGeneric("groups", function(x) standardGeneric("groups"))

#' replicates method for testClass
#'
#' @docType methods
#' @rdname replicates-methods
#' @keywords internal

setMethod("groups", signature = "countDat", function(x) {
  x@groups
})
