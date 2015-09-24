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
#' @export
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
    # 2 columns mean intercept + one other regression coefficient
    if (ncol(pcevObj$X) == 2) {
      pcevRes <- wilksPval(pcevObj, shrink, index)
    } else {
      pcevRes <- roysPval(pcevObj, shrink, index)
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
  class(pcevRes) <- "Pcev"

  # return results
  return(pcevRes)
}

computeVIMP <- function(pcevObj, list, signed=FALSE) {
  
  VIMP <- cor(pcevObj$Y, list$PCEV)[,1]
  if(!signed) {
    VIMP <- abs(VIMP)
  }
  
  return(VIMP)
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
#' @export
PcevClassical <- function(response, covariate, confounder) {
  if(is.null(confounder)) {
    structure(list(Y = response, 
                   X = model.matrix(~., as.data.frame(covariate)),
                   Z = c()), 
              class = "PcevClassical")
  } else {
    structure(list(Y = response, 
                   X = model.matrix(~., as.data.frame(covariate)), 
                   Z = model.matrix(~., as.data.frame(confounder))[,-1]), 
              class = "PcevClassical")
  }
  
}

#' @rdname PcevObj
#' @export 
PcevBlock <- function(response, covariate, confounder) {
  if(is.null(confounder)) {
    structure(list(Y = response, 
                   X = model.matrix(~., as.data.frame(covariate)), 
                   Z = c()), 
              class = "PcevBlock")
  } else {
    structure(list(Y = response, 
                   X = model.matrix(~., as.data.frame(covariate)), 
                   Z = model.matrix(~., as.data.frame(confounder))[,-1]), 
              class = "PcevBlock")
  }
  
}

# Documentation for the datasets----

#' Methylation values around BLK gene
#' 
#' A dataset containing methylation values for cell-separated samples. The 
#' methylation was measured using bisulfite sequencing. The data also contains 
#' the genomic position of these CpG sites, as well as a binary phenotype 
#' (whether the sample comes from a B cell).
#' 
#' Methylation was first measured at 24,068 sites, on 40 samples. Filtering was
#' performed to keep the 25\% most variable sites. See the vignette for more detail.
#' 
#' @format The data comes in three objects:
#' \describe{
#' \item{methylation}{Matrix of methylation values at 6,000 sites measured on 40 samples}
#' \item{pheno}{Vector of phenotype, indicating whether the sample comes from a B cell}
#' \item{position2}{Data frame recording the position of each CpG site along the BLK region}
#' }
#' @source Tomi Pastinen, McGill University
"methylation"

#' @rdname methylation
"pheno"

#' @rdname methylation
"position2"