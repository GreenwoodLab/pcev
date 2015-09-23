# permutation methods----

#' Permutation p-value
#' 
#' sdf
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance? 
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, index is a vector
#'   describing the block to which individual response variables correspond.
#' @param nperm The number of permutations to perform.
permutePval <- function(pcevObj) UseMethod("permutePval")

#' @describeIn permutePval
permutePval.default <- function(pcevObj) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevBlock"))
}

#' @describeIn permutePval
permutePval.PcevClassical <- function(pcevObj, shrink, index, nperm) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- nrow(pcevObj$X) - 1
  df2 <- N - df1 + 1
  initFstat <- (sum((mean(PCEV) - fit$fitted.values)^2)/df1)/(sum(fit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  permutationPvalues <- replicate(nperm, expr = {
    tmp <- pcevObj
    tmp$Y <- tmp$Y[sample(N), ]
    
    tmpRes <- estimatePcev(tmpRes, shrink)
    tmpPCEV <- tmp$Y %*% tmpRes$weights
    
    tmpFit <- lm.fit(tmp$X, tmpPCEV)
    tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
    return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
  })
  
  pvalue <- mean(permutationPvalues < initPval)
  
  results$pvalue <- pvalue
  
  return(results)
}

#' @describeIn permutePval
permutePval.PcevBlock <- function(pcevObj, shrink, index, nperm) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- nrow(pcevObj$X) - 1
  df2 <- N - df1 + 1
  initFstat <- (sum((mean(PCEV) - fit$fitted.values)^2)/df1)/(sum(fit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  permutationPvalues <- replicate(nperm, expr = {
    tmp <- pcevObj
    tmp$Y <- tmp$Y[sample(N), ]
    
    tmpRes <- estimatePcev(tmpRes, shrink)
    tmpPCEV <- tmp$Y %*% tmpRes$weights
    
    tmpFit <- lm.fit(tmp$X, tmpPCEV)
    tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
    return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
  })
  
  pvalue <- mean(permutationPvalues < initPval)
  
  results$pvalue <- pvalue
  
  return(results)
}

# Wilks methods----

#' Wilks' lambda exact test
#' 
#' sdf
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance? 
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, index is a vector
#'   describing the block to which individual response variables correspond.
wilksPval <- function(pcevObj) UseMethod("wilksPval")

#' @describeIn wilksPval
wilksPval.default <- function(pcevObj) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevBlock"))
}

#' @describeIn wilksPval
wilksPval.PcevClassical <- function(pcevObj, shrink, index) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  d <- results$largestRoot
  wilksLambda <- (N-p-1)/p * d
  
  df1 <- p
  df2 <- N-p-1
  pvalue <- pf(wilksLambda, df1, df2, lower.tail = FALSE)
  results$pvalue <- pvalue
  
  return(results)
}

#' @describeIn wilksPval
wilksPval.PcevBlock <- function(pcevObj, shrink, index) {
  results <- estimatePcev(pcevObj, shrink, index)
  N <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  PCEV <- pcevObj$Y %*% results$weights
  
  fit <- lm.fit(pcevObj$X, PCEV)
  beta <- fit$coefficients[2]
  h2Hat <- beta^2/(1 + beta^2)
  
  df1 <- p
  df2 <- N-p-1
  pvalue <- pf((N-p-1) * N * h2Hat/p, df1, df2, lower.tail = FALSE)
  
  results$pvalue <- pvalue
  
  return(results)
}

# Roy's largest root methods----

#' Roy's largest root exact test
#' 
#' sdf
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance? 
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, index is a vector
#'   describing the block to which individual response variables correspond.
roysPval <- function(pcevObj) UseMethod("roysPval")

#' @describeIn roysPval
roysPval.default <- function(pcevObj) {
  stop(strwrap("Pcev is currently not implemented for 
                   multiple covariates and an exact inference method"))
}

#' @describeIn roysPval
roysPval.PcevClassical <- function(PcevObj, shrink, index) {
  stop(strwrap("Pcev is currently not implemented for 
                   multiple covariates and an exact inference method"))
}

#' @describeIn roysPval
roysPval.PcevBlock <- function(PcevObj, shrink, index) {
  stop(strwrap("Pcev is currently not implemented for 
                   multiple covariates and an exact inference method"))
}