#' Wilks Lambda test.
#'
#' \code{WilksLambda} computes a p-value for the classical PCEV.
#'
#' The Wilks Lambda essentially tests whether the proportion of variance being
#' explained by the covariates is zero. This function computes the variance
#' decomposition which is used to compute the first PCEV. Note that this method
#' is currently only implemented for the case where there is only one covariate.
#'
#' @seealso \code{\link{computePCEV}}
#' @aliases wilkslambda wilksLambda Wilkslambda Wilks wilks
#' @param Y Matrix of values for a multivariate response.
#' @param x Vector of values for a covariate.
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#'   Default value is FALSE.
#' @return A list containing the variance components, the first PCEV, the
#'   p-value, etc.
WilksLambda <- function(Y, x, shrink = FALSE) {
  #initializing parameters
  rho <- NULL
  N <- nrow(Y)
  p <- ncol(Y)
  bar.Y <- colMeans(Y)

  #Variance decomposition
  fit <- lm.fit(cbind(rep_len(1, N), x), Y)
  Y.fit <- fit$fitted.values
  res <- Y - Y.fit
  Vr <- crossprod(res, Y)
  Vg <- crossprod(Y.fit, Y) - N * tcrossprod(bar.Y)
  if (shrink) {
    Vr <- Vr/(N-2)
    mu.E <- sum(diag(Vr))/p
    rho <- sum(apply(res, 1, function(a, b) {
                    sum(diag((tcrossprod(a)-b) %*% (tcrossprod(a)-b)))
                  }, b=Vr))
    rho <- rho/sum(diag((Vr-mu.E*diag(p)) %*% (Vr-mu.E*diag(p))))/N^2
    rho <- min(1,rho)
    Vr <- rho*mu.E*diag(p)+(1-rho)*Vr
    Vr <- Vr*(N-2)
  }
  #Computing PCEV
  temp <- eigen(Vr, symmetric=TRUE)
  Ur <- temp$vectors
  diagD <- temp$values
  value <- 1/sqrt(diagD)
  root.Vr <- Ur %*% diag(value) %*% t(Ur)
  m <- root.Vr %*% Vg %*% root.Vr
  temp1 <- eigen(m, symmetric=TRUE)
  PCEV <- root.Vr %*% temp1$vectors
  #Computing Wilks Lambda test statistic
  d <- temp1$values
  wilks.lambda <- (N-p-1)/p * d[1]
  #wilks.lambda = d[1]
  df1 <- p
  df2 <- N-p-1
  p.value <- pf(wilks.lambda, df1, df2, lower.tail = FALSE)

  return(list("environment" = Vr,
              "genetic" = Vg,
              "PCEV"=PCEV,
              "rootVr"=root.Vr,
              "values"=d,
              "pvalue"=p.value,
              "rho"=rho))
}

#' Principal Components of Explained Variance.
#'
#' \code{computePCEV} compute the first PCEV and tests its significance.
#'
#' This is the main function. It computes the PCEV using either the classical
#' method or the block-method. A p-value is also computed, testing the
#' significance of the PCEV. Note that the classical method is currently only
#' implemented for use with a single covariate.
#' @seealso \code{\link{WilksLambda}}
#' @param Y Matrix of values for a multivariate response.
#' @param x Vector of values for a covariate.
#' @param method Character string specifying which method to use: "all" or
#'   "block". Default value is "all".
#' @param index If type = "block", index is a vector describing the block to
#'   which the columns of Y correspond.
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#'   Default value is FALSE.
#' @return A list containing the first PCEV, the p-value, the estimate of the
#'   shrinkage factor, etc.


computePCEV <- function(Y, x, method = "all", index=NULL, shrink = FALSE) {
  #Check input
  if (!method %in% c("all", "block")){
    stop("Method should be \"all\" or \"block\"")
  }
  if (!is.numeric(index)) index <- NULL
  if(!is.numeric(Y) || !is.numeric(x)) {
    stop("Y and x should be numeric")
  }
  if (!is.logical(shrink)) shrink <- FALSE

  #Method by block
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

  #Classical method
  if (method == "all") {
    result <- WilksLambda(Y, x, shrink)
    weights <- result$PCEV[,1]
    Y.PCH <- Y %*% weights
    pvalue <- result$pvalue
    rho <- result$rho
  }

  #return results
  return(list("PCHvalues" = Y.PCH,
              "weights" = weights,
              "Pval"= pvalue,
              "rho"= rho))
}
