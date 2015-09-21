#' Estimation of PCEV
#' 
#' \code{estimatePcev} estimates the PCEV.
#' 
#' @seealso \code{\link{computePCEV}}
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance? 
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, index is a vector
#'   describing the block to which individual response variables correspond.
#' @return A list containing the variance components, the first PCEV, the 
#'   eigenvalues of \eqn{V_R^{-1}V_G} and the estimate of the shrinkage 
#'   parameter \eqn{\rho}
estimatePcev <- function(pcevObj) UseMethod("estimatePcev")

#' @describeIn  estimatePcev
estimatePcev.default <- function(pcevObj) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevBlock"))
}

#' @describeIn estimatePcev
estimatePcev.PcevClassical <- function(pcevObj, shrink) {
  #initializing parameters
  rho <- NULL
  Y <- pcevObj$Y
  N <- nrow(Y)
  p <- ncol(Y)
  bar.Y <- colMeans(Y)
  
  # Variance decomposition
  fit <- lm.fit(cbind(pcevObj$X, pcevObj$Z[,-1]), Y)
  Y.fit <- fit$fitted.values
  res <- Y - Y.fit
  Vr <- crossprod(res, Y)
  # Vg needs to be defined differently... depending only on covariates
  Vm <- crossprod(Y.fit, Y) - N * tcrossprod(bar.Y)
  
  # Shrinkage estimate of Vr
  if (shrink) {
    # This computation needs to be checked...
    # I don't think it corresponds to what is in Ibrahim et al.
    Vr <- Vr/(N-2)
    mu.E <- sum(diag(Vr))/p
    # Can this apply be replaced by something vectorized...?
    rho <- sum(apply(res, 1, function(a, b) {
      sum(diag(crossprod(tcrossprod(a)-b)))
    }, b=Vr))
    rho <- rho/sum(diag(crossprod(Vr - diag(mu.E))))/N^2
    rho <- min(1, rho)
    Vr <- diag(rho * mu.E) + (1 - rho)*Vr
    Vr <- Vr*(N - 2)
  }
  
  # Computing PCEV
  temp <- eigen(Vr, symmetric=TRUE)
  Ur <- temp$vectors
  diagD <- temp$values
  value <- 1/sqrt(diagD)
  root.Vr <- Ur %*% diag(value) %*% t(Ur)
  mainMatrix <- root.Vr %*% Vm %*% root.Vr
  temp1 <- eigen(mainMatrix, symmetric=TRUE)
  weights <- root.Vr %*% temp1$vectors
  d <- temp1$values
  
  return(list("residual" = Vr,
              "model" = Vm,
              "weights" = weights[,1],
              "rootVr" = root.Vr,
              "largestRoot" = d[1],
              "rho" = rho))
}

#' @describeIn estimatePcev
estimatePcev.PcevBlock <- function(pcevObj, shrink, index) {
  p <- ncol(pcevObj$Y)
  N <- nrow(pcevObj$Y)
  
  if (is.null(index) || p != length(index)) {
    stop("index should have length equal to number of response variables")
  }
  
  d <- length(unique(index))
  Ypcev <- matrix(NA, nrow = N, ncol = d)
  weights <- rep_len(0, p)
  rootVr <- list("first" = vector("list", d), 
                 "second" = NA)
  
  for (i in 1:d) {
    pcevObj_red <- pcevObj 
    pcevObj_red$Y <- pcevObj$Y[, index == i, drop = FALSE]
    class(pcevObj_red) <- "PcevClassical"
    result <- estimatePcev(pcevObj_red, shrink)
    weights[index==i] <- result$weights
    Ypcev[,i] <- Y.red %*% weights[index==i]
    rootVr$first[[i]] <- result$rootVr
  }
  pcevObj_total <- pcevObj
  pcevObj_total$Y <- Ypcev
  class(pcevObj_total) <- "PcevClassical"
  result <- estimatePcev(pcevObj_total, shrink)
  weight_step2 <- result$weights
  for (i in 1:d) {
    weights[index==i] <- weights[index==i]*weight_step2[i]
  }
  
  rootVr$second <- result$rootVr
  
  return(list("weights" = weights,
              "rootVr" = rootVr))
}