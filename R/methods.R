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
#' @param ... Extra parameters.
#' @export
permutePval <- function(pcevObj, ...) UseMethod("permutePval")

#' @describeIn permutePval
permutePval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevBlock"),
       call. = FALSE)
}

#' @describeIn permutePval
permutePval.PcevClassical <- function(pcevObj, shrink, index, nperm, ...) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- nrow(pcevObj$X) - 1
  df2 <- N - df1 + 1
  initFstat <- (sum((mean(PCEV) - initFit$fitted.values)^2)/df1)/(sum(initFit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  permutationPvalues <- replicate(nperm, expr = {
    tmp <- pcevObj
    tmp$Y <- tmp$Y[sample(N), ]
    
    tmpRes <- try(estimatePcev(tmp, shrink, index), silent=TRUE)
    if(inherits(tmpRes, "try-error")) {
      return(NA)
    } else {
      tmpPCEV <- tmp$Y %*% tmpRes$weights
      
      tmpFit <- lm.fit(tmp$X, tmpPCEV)
      tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
      return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
    }
  })
  
  pvalue <- mean(permutationPvalues < initPval)
  
  results$pvalue <- pvalue
  
  return(results)
}

#' @describeIn permutePval
permutePval.PcevBlock <- function(pcevObj, shrink, index, nperm, ...) {
  results <- estimatePcev(pcevObj, shrink, index)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- nrow(pcevObj$X) - 1
  df2 <- N - df1 + 1
  initFstat <- (sum((mean(PCEV) - initFit$fitted.values)^2)/df1)/(sum(initFit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  permutationPvalues <- replicate(nperm, expr = {
    tmp <- pcevObj
    tmp$Y <- tmp$Y[sample(N), ]
    
    tmpRes <- try(estimatePcev(tmp, shrink, index), silent=TRUE)
    if(inherits(tmpRes, "try-error")) {
      return(NA)
    } else {
      tmpPCEV <- tmp$Y %*% tmpRes$weights
      
      tmpFit <- lm.fit(tmp$X, tmpPCEV)
      tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
      return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
    }
  })
  
  pvalue <- mean(permutationPvalues < initPval, na.rm = TRUE)
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
#' @param ... Extra parameters.
#' @export
wilksPval <- function(pcevObj, ...) UseMethod("wilksPval")

#' @describeIn wilksPval
wilksPval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevBlock"),
       call. = FALSE)
}

#' @describeIn wilksPval
wilksPval.PcevClassical <- function(pcevObj, shrink, index, ...) {
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
wilksPval.PcevBlock <- function(pcevObj, shrink, index, ...) {
  N <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  if (N - p <= 1) stop("To perform an exact test, we need N - p > 1", 
                       call. = FALSE)
  results <- estimatePcev(pcevObj, shrink, index)
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
#' @param ... Extra parameters.
#' @export 
roysPval <- function(pcevObj, ...) UseMethod("roysPval")

#' @describeIn roysPval
roysPval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical"),
       call. = FALSE)
}

#' @describeIn roysPval
roysPval.PcevClassical <- function(pcevObj, shrink, index, ...) {
  
  if (!requireNamespace("RMTstat", quietly = TRUE)) {
    stop("RMTstat needs to be installed in order to use the exact method with multiple covariates.", 
         call. = FALSE)
  }
  
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  q <- ncol(pcevObj$X) - 1
  d <- results$largestRoot
  theta <- d / (1 + d)
  
  nuH <- q
  nuE <- N - q - 1
  s <- min(p, nuH)
  m <- 0.5 * (abs(p - nuH) - 1)
  n <- 0.5 * (nuE - p - 1)
  one_third <- 1/3
  
  N <- 2 * (s + m + n) + 1 
  gamma <- 2 * asin( sqrt( (s - 0.5)/N ) )
  phi <- 2*asin( sqrt( (s + 2*m + 0.5)/N ) )
  
  mu <- 2 * log(tan( 0.5*(phi + gamma)))
  sigma <- (16/N^2)^one_third * ( sin(phi+gamma)^2*sin(phi)*sin(gamma)) ^(-1*one_third)
  
  TW <- (log(theta/(1-theta)) - mu)/sigma
  
  pvalue <- RMTstat::ptw(TW, beta=1, lower.tail = FALSE, log.p = FALSE)
  results$pvalue <- pvalue
  
  return(results)
  
}

#' @describeIn roysPval
roysPval.PcevBlock <- function(pcevObj, shrink, index, ...) {
  stop(strwrap("Pcev is currently not implemented for 
                   multiple covariates, estimation with blocks and an exact inference method"),
       call. = FALSE)
}

# Print method----

print.Pcev <- function(pcevRes, ...) {
  # Provide a summary of the results
  N <- nrow(pcevRes$pcevObj$Y)
  p <- ncol(pcevRes$pcevObj$Y)
  q <- ncol(pcevRes$pcevObj$X)
  if(q == 2) {
    exact <- "Roy's largest root test)"
  } else {
    exact <- "Wilks' lambda test)"
  }
  
  cat("\nPrincipal component of explained variance\n")
  cat("\nEstimation method:", pcevRes$methods[1])
  cat("\nInference method:", pcevRes$methods[2])
  if(pcevRes$methods[2] == "exact") {
    cat("\n(performed using", exact)
  }
  pvalue <- pcevRes$pvalue
  if(pvalue == 0 && pcevRes$methods[2] == "permutation") {
    pvalue <- paste0("< ", 1/pcevRes$nperm)
  }
  cat("\nP-value obtained:", pvalue, "\n")
  
  cat("\nVariable importance factors")
  if(p > 10) cat(" (truncated)\n") else cat("\n")
  cat(format(head(pcevRes$VIMP, n = 10), digits = 3), 
      "\n\n")
  
}