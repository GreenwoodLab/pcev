#' Principal Components of Explained Variance
#' 
#' \code{computePCEV} computes the first PCEV and tests its significance.
#' 
#' This is the main function. It computes the PCEV using either the classical 
#' method or the block approach. A p-value is also computed, testing the 
#' significance of the PCEV. Note that the classical method is currently only 
#' implemented for use with a single covariate.
#' 
#' @seealso \code{\link{estimatePcev}}
#' @param response A matrix of response variables.
#' @param covariate A matrix or a data frame of covariates.
#' @param confounder A matrix or data frame of confounders
#' @param estimation Character string specifying which estimation method to use:
#'   \code{"all"} or \code{"block"}. Default value is \code{"all"}.
#' @param inference Character string specifying which inference method to use: 
#'   \code{"exact"} or \code{"permutation"}. Default value is \code{"exact"}.
#' @param index If \code{estimation = "block"}, index is a vector describing the block 
#'   to which individual response variables correspond.
#' @param shrink Should we use a shrinkage estimate of the residual variance? 
#'   Default value is \code{FALSE}..
#' @param nperm The number of permutations to perform if \code{inference = 
#'   "permutation"}
#' @return A list containing the first PCEV, the p-value, the estimate of the 
#'   shrinkage factor, etc.


computePCEV <- function(response, covariate, confounder = NULL, 
                        estimation = "all", inference = "exact", 
                        index = NULL, shrink = FALSE, nperm = 1000) {
  # Check input
  if (!estimation %in% c("all", "block")){
    stop("Estimation method should be \"all\" or \"block\"")
  }
  if (!inference %in% c("exact", "permutation")){
    stop("Inference method should be \"exact\" or \"permutation\"")
  }
  if (!is.numeric(index)) index <- NULL
  if (!is.logical(shrink)) shrink <- FALSE
  
  # Create pcev objects
  if (estimation == "all") {
    pcevObj <- PcevClassical(response,
                             covariate,
                             confounder)
  }
  if (estimation == "block") {
    pcevObj <- PcevBlock(response, 
                         covariate,
                         confounder)
  }
  
  # Perform estimation and inference
  if (inference == "permutation") {
    pcevRes <- permutePval(pcevObj, shrink, index, nperm)
  } else {
    if (ncol(covariate) == 1) {
      # This condition is maybe too loose
      # What about categorical variables?
      pcevRes <- wilksPval(pcevObj, shrink, index)
    } else {
      pcevRes <- roysPval(pcevObj, shrink, index)
      # Need to implement this
    }
  }

  # Method by block
  if (method == "block") {
    if (is.null(index) || ncol(Y) != length(index)) {
      stop("index should have length equal to ncol(Y)")
    }

    d <- length(unique(index))
    Y.PCH <- matrix(NA, nrow=nrow(Y), ncol=d)
    weights <- rep_len(0,ncol(Y))
    for (i in 1:d) {
      Y.red <- Y[, index==i, drop = FALSE]
      result <- WilksLambda(Y.red, x, shrink)
      weights[index==i] <- result$PCEV[,1]
      Y.PCH[,i] <- Y.red %*% weights[index==i]
    }
    result <- WilksLambda(Y.PCH, x, shrink)
    weight.step2 <- result$PCEV[,1]
    for (i in 1:d) {
      weights[index==i] <- weights[index==i]*weight.step2[i]
    }
    Y.PCH <- Y %*% weights
    pvalue <- result$pvalue
    rho <- result$rho
  }

  # Classical method
  if (method == "all") {
    result <- WilksLambda(Y, x, shrink)
    weights <- result$PCEV[,1]
    Y.PCH <- Y %*% weights
    pvalue <- result$pvalue
    rho <- result$rho
  }

  # return results
  return(list("PCHvalues" = Y.PCH,
              "weights" = weights,
              "Pval"= pvalue,
              "rho"= rho))
}

# Constructor functions----

#' Constructor functions for the different pcev objects
#' 
#' \code{PcevClassical} and \code{PcevBlock} create the pcev objects from the
#' provided data that are necessary to compute the PCEV according to the user's
#' parameters.
#' 
#' @seealso \code{\link{estimatePcev}}, \code{\link{computePCEV}}
#' @param response A matrix of response variables.
#' @param covariate A matrix or a data frame of covariates.
#' @param confounder A matrix or data frame of confounders
#' @return A pcev object, of the class that corresponds to the estimation 
#'   method. These objects are lists that essentially contain the data necessary
#'   for computation.
#' @name PcevObj
NULL

#' @rdname PcevObj
PcevClassical <- function(response, covariate, confounder) {
  structure(list(Y = response, 
                 X = model.matrix(~., covariate), 
                 Z = model.matrix(~., confounder)), 
            class = "PcevClassical")
}

#' @rdname PcevObj
PcevBlock <- function(response, covariate, confounder) {
  structure(list(Y = response, 
                 X = model.matrix(~., covariate), 
                 Z = model.matrix(~., confounder)), 
            class = "PcevBlock")
}