#' Principal Component of Explained Variance
#'
#' \code{computePCEV} computes the first PCEV and tests its significance.
#'
#' This is the main function. It computes the PCEV using either the classical
#' method, block approach or singular. A p-value is also computed, testing the
#' significance of the PCEV.
#'
#' The p-value is computed using either a permutation approach or an exact test.
#' The implemented exact tests use Wilks' Lambda (only for a single covariate)
#' or Roy's Largest Root. The latter uses Johnstone's approximation to the null
#' distribution. Note that for the block approach, only p-values obtained from a
#' permutation procedure are available.
#'
#' When \code{estimation = "singular"}, the p-value is computed using a
#' heuristic: using the method of moments and a small number of permutations
#' (i.e. 25), a location-scale family of the Tracy-Widom distribution of order 1
#' is fitted to the null distribution. This fitted distribution is then used to
#' compute p-values.
#'
#' When \code{estimation = "block"}, there are three different ways of
#' specifying the blocks: 1) if \code{index} is a vector of the same length as
#' the number of columns in \code{response}, then it is used to match each
#' response to a block. 2) If \code{index} is a single positive integer, it is
#' understood as the number of blocks, and each response is matched to a block
#' randomly. 3) If \code{index = "adaptive"} (the default), the number of blocks
#' is chosen so that there are about n/2 responses per block, and each response
#' is match to a block randomly. All other values of \code{index} should result
#' in an error.
#'
#' By default, missing values are not allowed. This can be relaxed with
#' \code{na_action}. If \code{na_action = "omit"}, then all rows with at least
#' one missing value will be removed from \code{response} before computation. If
#' \code{na_action = "column"}, then the estimation of the linear model
#' parameters is done column-wise with the non-missing value. This approach
#' maximises the information. Note that missing values are still not allowed in
#' \code{covariate} and \code{confounder}.
#'
#' @seealso \code{\link{estimatePcev}}
#' @param response A matrix of response variables.
#' @param covariate An array or a data frame of covariates.
#' @param confounder An array or data frame of confounders.
#' @param estimation Character string specifying which estimation method to use:
#'   \code{"all"}, \code{"block"} or \code{"singular"}. Default value is
#'   \code{"all"}.
#' @param inference Character string specifying which inference method to use:
#'   \code{"exact"} or \code{"permutation"}. Default value is \code{"exact"}.
#' @param index Only used if \code{estimation = "block"}. Default value is
#'   \code{"adapative"}. See details.
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#'   Default value is \code{FALSE}.
#' @param nperm The number of permutations to perform if \code{inference =
#'   "permutation"} or for the Tracy-Widom empirical estimate (if
#'   \code{estimation = "singular"}).
#' @param na_action how NAs are treated. The default is to raise an error. See
#'   details.
#' @param Wilks Should we use a Wilks test instead of Roy's largest test? This
#'   is only implemented for a single covariate and with \code{estimation =
#'   "all"}.
#' @return An object of class \code{Pcev} containing the first PCEV, the
#'   p-value, the estimate of the shrinkage factor, etc.
#' @examples
#' set.seed(12345)
#' Y <- matrix(rnorm(100*20), nrow=100)
#' X <- rnorm(100)
#' pcev_out <- computePCEV(Y, X)
#' pcev_out2 <- computePCEV(Y, X, shrink = TRUE)
#' @export
#' @importFrom stats coefficients cor lm.fit model.matrix optim pf
computePCEV <- function(response, covariate, confounder, 
                        estimation = c("all", "block", "singular"), 
                        inference = c("exact", "permutation"), 
                        index = "adaptive", shrink = FALSE, nperm = 1000, 
                        na_action = "fail", Wilks = FALSE) {
  # Check input
  estimation <- tryCatch(match.arg(estimation),
                         error = function(c) {
                           stop("Estimation method should be \"all\", \"block\" or \"singular\"", 
                                call. = FALSE)
                         })
  
  inference <- tryCatch(match.arg(inference),
                        error = function(c) {
                          stop("Inference method should be \"exact\" or \"permutation\"", 
                               call. = FALSE)
                        })
  
  na_action <- tryCatch(match.arg(na_action, choices = c("fail", "omit", "column")),
                        error = function(c) {
                          stop("NA action should be \"fail\", \"omit\" or \"column\"", 
                               call. = FALSE)
                        })
  
  if (!is.matrix(response)) {
    stop("The response variables should be passed as a matrix.", call. = FALSE)
  }
  if (missing(confounder)) confounder <- NULL
  
  # We don't allow for missing values in covariate and confounder
  if (anyNA(covariate) || anyNA(confounder)) {
    stop("Missing values in covariate and confounder are not allowed", 
         call. = FALSE)
  }
  if (anyNA(response) && na_action == "fail") {
    stop("Missing values in covariate and confounder are not allowed", 
         call. = FALSE)
  }
  
  # If user specifies na_action = "omit", remove incomplete cases
  if (anyNA(response) && na_action == "omit") {
    response <- na.omit(response)
    covariate <- as.matrix(covariate)[-attr(response, "na.action"), ,
                                      drop = FALSE]
    confounder <- if (is.null(confounder)) confounder else confounder[-attr(response, "na.action"), ,
                                                                      drop = FALSE]
  }
  
  # If user gives index, we should do block estimation
  if (!identical(index, "adaptive")) {
    estimation <- "block"
    message("Selecting estimation by block.")
  }
  
  # Adapative selection of blocks
  if (estimation == "block") {
    if (identical(index, "adaptive")) {
      b <- round(2*ncol(response)/nrow(response))
    }
    if (is.numeric(index) && length(index) == 1) {
      b <- index
      if (b < 0 || b > ncol(response)) {
        stop("Invalid number of blocks", call. = FALSE)
      }
    }
    if (!is.numeric(index) || length(index) == 1) {
      index <- sample(1:b, ncol(response), replace = TRUE)
    }
  }
  
  if (!is.logical(shrink)) shrink <- FALSE
  if (!is.logical(Wilks)) Wilks <- FALSE
  
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
  if (estimation == "singular") {
    pcevObj <- PcevSingular(response,
                            covariate,
                            confounder)
  }
  
  if (Wilks && ncol(pcevObj$X) != 2) {
    warning("Wilks can only be applied with a single covariate.")
    Wilks <- FALSE
  }
  
  # Add information about missing values to pcevObj
  pcevObj$overall <- (na_action != "column")
  
  # Perform estimation and inference
  if (inference == "permutation") {
    pcevRes <- permutePval(pcevObj, shrink, index, nperm)
  } else {
    if (Wilks) {
      pcevRes <- wilksPval(pcevObj, shrink, index)
    } else {
      pcevRes <- roysPval(pcevObj, shrink, index, nperm)
    }
  }
  
  # Compute PCEV
  pcevRes$PCEV <- pcevObj$Y %*% pcevRes$weights
  
  # Compute variable importance
  pcevRes$VIMP <- computeVIMP(pcevObj, pcevRes)
  
  pcevRes$pcevObj <- pcevObj 
  pcevRes$methods <- c(estimation, inference)
  names(pcevRes$methods) <- c("Estimation", "Inference")
  pcevRes$nperm <- nperm 
  pcevRes$Wilks <- Wilks
  pcevRes$shrink <- shrink
  class(pcevRes) <- "Pcev"
  
  # return results
  return(pcevRes)
}
