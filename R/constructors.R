###########################
# Constructor functions----

#' Constructor functions for the different pcev objects
#' 
#' \code{PcevClassical}, \code{PcevBlock} and \code{PcevSingular} create the pcev objects from the 
#' provided data that are necessary to compute the PCEV according to the user's 
#' parameters.
#' 
#' @seealso \code{\link{estimatePcev}}, \code{\link{computePCEV}}
#' @param response A matrix of response variables.
#' @param covariate A matrix or a data frame of covariates.
#' @param confounder A matrix or data frame of confounders
#' @return A pcev object, of the class that corresponds to the estimation 
#'   method. These objects are lists that contain the data necessary for
#'   computation.
#' @name PcevObj
NULL

#' @rdname PcevObj
#' @export
PcevClassical <- function(response, covariate, confounder) {
  if (is.null(confounder)) {
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
  if (is.null(confounder)) {
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

#' @rdname PcevObj
#' @export
PcevSingular <- function(response, covariate, confounder) {
  if (is.null(confounder)) {
    structure(list(Y = response, 
                   X = model.matrix(~., as.data.frame(covariate)),
                   Z = c()), 
              class = "PcevSingular")
  } else {
    structure(list(Y = response, 
                   X = model.matrix(~., as.data.frame(covariate)), 
                   Z = model.matrix(~., as.data.frame(confounder))[,-1]), 
              class = "PcevSingular")
  }
  
}

