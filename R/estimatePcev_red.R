#' Estimation of PCEV (reduced computation)
#' 
#' \code{estimatePcev_red} estimates the PCEV by reusing previous estimate of 
#' the residual variance. This is used to speed up computations of the 
#' permutation procedures.
#' 
#' @seealso \code{\link{computePCEV}}
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param root_Vr Previous estimate of the square root of the inverse of the 
#'   residual variance
#' @param root_Vr_list List of previous estimates of the square root of the 
#'   inverse of the residual variances. Requires one for each block and one for 
#'   the second stage of the estimation.
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, index is a 
#'   vector describing the block to which individual response variables 
#'   correspond.
#' @param ... Extra parameters.
#' @return A list containing the first PCEV and the largest eigenvalue of
#'   \eqn{V_R^{-1}V_G}.
#' @export
estimatePcev_red <- function(pcevObj, ...) UseMethod("estimatePcev_red")

#' @describeIn  estimatePcev_red
estimatePcev_red.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevBlock"))
}

#' @describeIn estimatePcev_red
estimatePcev_red.PcevClassical <- function(pcevObj, index, root_Vr, ...) {
  #initializing parameters
  Y <- pcevObj$Y
  N <- nrow(Y)
  p <- ncol(Y)
  
  # Variance decomposition
  fit <- lm.fit(cbind(pcevObj$X, pcevObj$Z), Y)
  Yfit <- fit$fitted.values
  res <- Y - Yfit
  fit_confounder <- lm.fit(cbind(rep_len(1, N), pcevObj$Z), Y)
  Yfit_confounder <- fit_confounder$fitted.values
  
  Vm <- crossprod(Yfit - Yfit_confounder, Y)
  
  # Reuse previous estimate of root_Vr
  mainMatrix <- root_Vr %*% Vm %*% root_Vr
  temp1 <- eigen(mainMatrix, symmetric=TRUE)
  weights <- root_Vr %*% temp1$vectors
  d <- temp1$values
  
  return(list("weights" = weights[,1, drop=FALSE],
              "largestRoot" = d[1]))
}


#' @describeIn estimatePcev_red
estimatePcev_red.PcevBlock <- function(pcevObj, index, root_Vr_list, ...) {
  p <- ncol(pcevObj$Y)
  N <- nrow(pcevObj$Y)
  
  if (is.null(index) || p != length(index)) {
    stop("index should have length equal to number of response variables")
  }
  
  d <- length(unique(index))
  if(d > N && ncol(pcevObj$X) != 2) {
    warning("It is recommended to have a number of blocks smaller than the number of observations")
  }
  Ypcev <- matrix(NA, nrow = N, ncol = d)
  weights <- rep_len(0, p)
  
  counter <- 0
  for (i in unique(index)) {
    counter <- counter + 1
    pcevObj_red <- pcevObj 
    pcevObj_red$Y <- pcevObj$Y[, index == i, drop = FALSE]
    class(pcevObj_red) <- "PcevClassical"
    result <- estimatePcev_red(pcevObj_red, root_Vr_list$first[[counter]])
    weights[index==i] <- result$weights
    Ypcev[,counter] <- pcevObj_red$Y %*% weights[index==i]
  }
  
  pcevObj_total <- pcevObj
  pcevObj_total$Y <- Ypcev
  class(pcevObj_total) <- "PcevClassical"
  
  if (ncol(pcevObj_total$X) == 2) {
    fit_total <- lm.fit(pcevObj_total$X, pcevObj_total$Y)
    beta_total <- coefficients(fit_total)[2,]
    weight_step2 <- beta_total/crossprod(beta_total)
    
  } else {
    result <- estimatePcev(pcevObj_total, root_Vr_list$second)
    weight_step2 <- result$weights
  }
  
  counter <- 0
  for (i in unique(index)) {
    counter <- counter + 1
    weights[index==i] <- weights[index==i]*weight_step2[counter]
  }
  
  return(list("weights" = weights))
}