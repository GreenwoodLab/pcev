# Permutation methods----

#' Permutation p-value
#' 
#' Computes a p-value using a permutation procedure.
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or \code{PcevSingular} 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, \code{index} is a
#'   vector describing the block to which individual response variables
#'   correspond.
#' @param nperm The number of permutations to perform.
#' @param ... Extra parameters.
#' @export
permutePval <- function(pcevObj, ...) UseMethod("permutePval")

#' @rdname permutePval
permutePval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical, PcevSingular or PcevBlock"),
       call. = FALSE)
}

#' @rdname permutePval
permutePval.PcevClassical <- function(pcevObj, shrink, index, nperm, ...) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- ncol(pcevObj$X) - 1
  df2 <- N - ncol(pcevObj$X)
  initFstat <- (sum((mean(PCEV) - initFit$fitted.values)^2)/df1)/(sum(initFit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  if (nperm) {
    permutationPvalues <- replicate(nperm, expr = {
      tmp <- pcevObj
      tmp$Y <- tmp$Y[sample(N), ]
      
      tmpRes <- try(estimatePcev(tmp, shrink, index), silent = TRUE)
      if (inherits(tmpRes, "try-error")) {
        return(NA)
      } else {
        tmpPCEV <- tmp$Y %*% tmpRes$weights
        
        tmpFit <- lm.fit(tmp$X, tmpPCEV)
        tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
        return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
      }
    })
    
    pvalue <- mean(permutationPvalues < initPval)
  } else {
    pvalue <- NA
  }
  
  results$pvalue <- pvalue
  
  return(results)
}

#' @rdname permutePval
permutePval.PcevBlock <- function(pcevObj, shrink, index, nperm, ...) {
  results <- estimatePcev(pcevObj, shrink, index)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- ncol(pcevObj$X) - 1
  df2 <- N - ncol(pcevObj$X)
  initFstat <- (sum((mean(PCEV) - initFit$fitted.values)^2)/df1)/(sum(initFit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  if (nperm) {
    permutationPvalues <- replicate(nperm, expr = {
      tmp <- pcevObj
      tmp$Y <- tmp$Y[sample(N), ]
      
      tmpRes <- try(estimatePcev(tmp, shrink, index), silent = TRUE)
      if (inherits(tmpRes, "try-error")) {
        return(NA)
      } else {
        tmpPCEV <- tmp$Y %*% tmpRes$weights
        
        tmpFit <- lm.fit(tmp$X, tmpPCEV)
        tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
        return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
      }
    })
    
    pvalue <- mean(permutationPvalues < initPval)
  } else {
    pvalue <- NA
  }
  
  results$pvalue <- pvalue
  return(results)
}

#' @rdname permutePval
permutePval.PcevSingular <- function(pcevObj, shrink, index, nperm, ...) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- ncol(pcevObj$X) - 1
  df2 <- N - ncol(pcevObj$X)
  initFstat <- (sum((mean(PCEV) - initFit$fitted.values)^2)/df1)/(sum(initFit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  if (nperm) {
    permutationPvalues <- replicate(nperm, expr = {
      tmp <- pcevObj
      tmp$Y <- tmp$Y[sample(N), ]
      
      tmpRes <- try(estimatePcev(tmp, shrink, index), silent = TRUE)
      if (inherits(tmpRes, "try-error")) {
        return(NA)
      } else {
        tmpPCEV <- tmp$Y %*% tmpRes$weights
        
        tmpFit <- lm.fit(tmp$X, tmpPCEV)
        tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
        return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
      }
    })
    
    pvalue <- mean(permutationPvalues < initPval)
  } else {
    pvalue <- NA
  }
  
  results$pvalue <- pvalue
  
  return(results)
}

###################
# Wilks methods----

#' Wilks' lambda exact test
#' 
#' Computes a p-value using Wilks' Lambda.
#' 
#' The null distribution of this test statistic is only known in the case of a
#' single covariate, and therefore this is the only case implemented.
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, \code{index} is
#'   a vector describing the block to which individual response variables 
#'   correspond.
#' @param ... Extra parameters.
#' @export
wilksPval <- function(pcevObj, ...) UseMethod("wilksPval")

#' @rdname wilksPval
wilksPval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical, PcevBlock or PcevSingular"),
       call. = FALSE)
}

#' @rdname wilksPval
wilksPval.PcevClassical <- function(pcevObj, shrink, index, ...) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  d <- results$largestRoot
  wilksLambda <- (N - p - 1)/p * d
  
  df1 <- p
  df2 <- N - p - 1
  pvalue <- pf(wilksLambda, df1, df2, lower.tail = FALSE)
  results$pvalue <- pvalue
  
  return(results)
}

#' @rdname wilksPval
wilksPval.PcevSingular <- function(pcevObj, shrink, index, ...) {
  stop(strwrap("Wilks test is only available for estimation = all."),
       call. = FALSE)
}

#' @rdname wilksPval
wilksPval.PcevBlock <- function(pcevObj, shrink, index, ...) {
  stop("Only inference = \"permutation\" is available for the block approach.",
       call. = FALSE)
}

################################
# Roy's largest root methods----

#' Roy's largest root exact test
#' 
#' In the classical domain of PCEV applicability this function uses Johnstone's
#' approximation to the null distribution of ' Roy's Largest Root statistic.
#' It uses a location-scale variant of the Tracy-Wildom distribution of order 1.
#' 
#'For singular PCEV, where number of variables is higher than the number of observations,
#'one can choose between two aprroximations. First one, called 'Wishart', is based
#'on the first term of expansion of the joint distribution of roots for singular beta
#'ensemble. Second version of test, \code{estimation = "TW"}, assumes that largest root statistics follows
#'Tracy-Wildom as in classical case. See vignette for the domains of the applicability
#'of these tests.
#'    
#'   
#' 
#' Note that if \code{shrink} is set to \code{TRUE}, the location-scale 
#' parameters are estimated using a small number of permutations.
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, \code{index} is
#'   a vector describing the block to which individual response variables 
#'   correspond
#' @param nperm Number of permutations for Tracy-Widom empirical estimate.
#' @param ... Extra parameters.
#' @export
roysPval <- function(pcevObj, ...) UseMethod("roysPval")

#' @rdname roysPval
roysPval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevSingular"),
       call. = FALSE)
}

#' @rdname roysPval
roysPval.PcevClassical <- function(pcevObj, shrink, index, ...) {
  
  results <- estimatePcev(pcevObj, shrink)
  n <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  q <- ncol(pcevObj$X) - 1
  d <- results$largestRoot
  params <- JohnstoneParam(p, n - q, q)
  
  if (shrink) {
    # Estimate the null distribution using 
    # permutations and MLE
    null_dist <- replicate(25, expr = {
      tmp <- pcevObj
      tmp$Y <- tmp$Y[sample(n), ]
      
      tmpRes <- try(estimatePcev(tmp, shrink, index), silent = TRUE)
      if (inherits(tmpRes, "try-error")) {
        return(NA)
      } else {
        return(tmpRes$largestRoot)
      }
    })
    
    # Fit a location-scale version of TW distribution
    # Note: likelihood may throw warnings at some evaluations, 
    # which is OK
    oldw <- getOption("warn")
    options(warn = -1)
    res <- optim(c(params[1], params[2]), function(param) logLik(param, log(null_dist)), 
                 control = list(fnscale = -1))
    options(warn = oldw)
    
    mu1 <- res$par[1]
    sigma1 <- res$par[2]
    TW <- (log(d) - mu1)/sigma1
    
    pvalue <- RMTstat::ptw(TW, beta = 1, lower.tail = FALSE, log.p = FALSE)
  } else {
    # TW <- (log(theta/(1-theta)) - mu)/sigma
    TW <- (log(d) - params[1])/params[2]
    
    pvalue <- RMTstat::ptw(TW, beta = 1, lower.tail = FALSE, log.p = FALSE)
  }
  
  results$pvalue <- pvalue
  
  return(results)
  
}

#' @rdname roysPval
roysPval.PcevSingular <- function(pcevObj, shrink, index, nperm, ...) {
  results <- estimatePcev(pcevObj, shrink)
  n <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  q <- ncol(pcevObj$X) 
  d <- results$largestRoot
  if (missing(nperm)) nperm <- 50

  null_dist <- replicate(nperm, expr = {
    tmp <- pcevObj
    tmp$Y <- tmp$Y[sample(n), ]
    
    tmpRes <- try(estimatePcev(tmp, shrink, index), silent = TRUE)
    if (inherits(tmpRes, "try-error")) {
      return(NA)
    } else {
      return(tmpRes$largestRoot)
    }
  })
  
  # Use method of moments
  muTW <- -1.2065335745820
  sigmaTW <- sqrt(1.607781034581)
  
  muS <- mean(log(null_dist))
  sigmaS <- stats::sd(log(null_dist))
  
  sigma1 <- sigmaS/sigmaTW
  mu1 <- muS - sigma1 * muTW
  
  TW <- (log(d) - mu1)/sigma1
  pvalue <- RMTstat::ptw(TW, beta = 1, lower.tail = FALSE, log.p = FALSE)

  
  results$pvalue <- pvalue
  
  return(results)
}

#' @rdname roysPval
roysPval.PcevBlock <- function(pcevObj, shrink, index, ...) {
  stop("Only inference = \"permutation\" is available for the block approach.",
       call. = FALSE)
}

##################
# Print method----

#' @export
print.Pcev <- function(x, ...) {
  # Provide a summary of the results
  N <- nrow(x$pcevObj$Y)
  p <- ncol(x$pcevObj$Y)
  q <- ncol(x$pcevObj$X)
  if (x$Wilks) {
    exact <- "Wilks' lambda test)"
  } else {
    exact <- "Roy's largest root test)"
  }
  
  cat("\nPrincipal component of explained variance\n")
  cat("\n", N, "observations,", p, "response variables\n")
  cat("\nEstimation method:", x$methods[1])
  cat("\nInference method:", x$methods[2])
  if (x$methods[2] == "exact") {
    cat("\n(performed using", exact)
  }
  pvalue <- x$pvalue
  if (!is.na(pvalue)) {
    if (pvalue == 0) {
      if (x$methods[2] == "permutation") {
        pvalue <- paste0("< ", 1/x$nperm)
      }
      if (x$methods[1] == "exact") {
        pvalue <- "~ 0"
      }
    }
  }
  cat("\nP-value obtained:", pvalue, "\n")
  if (x$shrink) cat("\nShrinkage parameter rho was estimated at", x$rho, "\n")
  cat("\nVariable importance factors")
  if (p > 10) cat(" (truncated)\n") else cat("\n")
  cat(format(sort(x$VIMP, decreasing = TRUE)[1:10], digits = 3), 
      "\n\n")
  
}

########################################################
# Utility functions for null distribution estimation----
dtw_ls <- function(x, mu, sigma, beta=1, log=FALSE) {
  x1 <- (x - mu)/sigma
  # values for dtw are only available for the 
  # interval -10 to 6
  x1 <- as.numeric(x1 >= -10)*x1 - 10*as.numeric(x1 < -10)
  x1 <- as.numeric(x1 <= 6)*x1 + 6*as.numeric(x1 > 6)
  return(RMTstat::dtw(x1, beta, log)/sigma)
}

logLik <- function(param, data) {
  mu <- param[1]
  sigma <- param[2]
  data <- data[!is.na(data)]
  
  lL <- sum(log(dtw_ls(data, mu, sigma, log = FALSE)))
  
  return(lL)
}

blockMatrixDiagonal <- function(matrix_list) {  
  
  dimensions <- sapply(matrix_list, nrow)
  finalDimension <- sum(dimensions)
  finalMatrix <- matrix(0, nrow = finalDimension, ncol = finalDimension)
  index <- 1
  for (k in 1:length(dimensions)) {
    finalMatrix[index:(index + dimensions[k] - 1), 
                index:(index + dimensions[k] - 1)] <- matrix_list[[k]]
    index <- index + dimensions[k]
  }
  
  return(finalMatrix)
}

JohnstoneParam <- function(p, m, n) {
  s_serif <- min(n, p)
  m_serif <- 0.5*(abs(n - p) - 1)
  n_serif <- 0.5*(abs(m - p) - 1)
  N_serif <- 2*(s_serif + m_serif + n_serif) + 1
  one_third <- 1/3
  
  gamma <- 2 * asin( sqrt( (s_serif - 0.5)/N_serif ) )
  phi <- 2 * asin( sqrt( (s_serif + 2*m_serif + 0.5)/N_serif ) )
  
  mu <- 2 * log(tan( 0.5*(phi + gamma)))
  sigma <- (16/N_serif^2)^one_third * (sin(phi + gamma)^2 * sin(phi) * sin(gamma))^(-1*one_third)
  
  return(c(mu, sigma))
}
