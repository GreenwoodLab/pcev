roysPval.PcevSingular <- function(pcevObj, shrink, distrib = c("Wishart", "TW"), index, ...) {
  results <- estimatePcev(pcevObj, shrink)
  n <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  q <- ncol(pcevObj$X) 
  d <- results$largestRoot
  distrib <- match.arg(distrib)
  
  if (distrib == 'Wishart'){
    if(!requireNamespace("rootWishart")) {
      stop("This requires the rootWishart package. You should use distrib = TW instead.")
    }
    resid <- results$residual
    s <- q - 1
    r <- n - s
    trRessq <- sum(resid ^ 2) 
    trRes <- sum(diag(resid))
    b <- (trRes/r/p)^2 / ((trRessq-trRes^2/r)/(r-1)/(r+2)/p)
    pvalue <- 1-rootWishart::singleWishart(d*p*b, r, s, mprec = TRUE) 
  } else {
    null_dist <- replicate(25, expr = {
      tmp <- pcevObj
      tmp$Y <- tmp$Y[sample(n), ]
      
      tmpRes <- try(estimatePcev(tmp, shrink, index), silent=TRUE)
      if(inherits(tmpRes, "try-error")) {
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
    pvalue <- RMTstat::ptw(TW, beta=1, lower.tail = FALSE, log.p = FALSE)
  }
  
  results$pvalue <- pvalue
  
  return(results)
  
}

#Function to find the level in the factor variable X
#having the strongest effect on Y.

#Anova helper-function
#' @importFrom stats anova as.formula fitted.values lm
linregfn <- function(y, covars){
  Dis <- covars
  y <- unlist(y)
  modelstring <- "y ~ covars"  
  fit <- lm(as.formula(modelstring))
  resid <- fitted.values(fit)
  sfit <- summary(fit)
  anovF <- anova(fit)
  nc1 <- nrow(sfit$coefficients)
  nc2 <- nrow(anovF) - 1
  pvals <- c(sfit$coefficients[2:nc1, 1], sfit$coefficients[2:nc1,4], anovF[1:nc2,5])
  names(pvals) <- c(paste0('c_', rownames(sfit$coefficients)[2:nc1]), rownames(sfit$coefficients)[2:nc1], rownames(anovF)[1:nc2])
  return(list(pvals = pvals, 
              resid = resid))
}

#function to find the most deviant subpopulation of samples
# @export
fimp <- function(data, covars){
  pvals <- list()
  l <- unique(covars)
  for (i in l) {
    cnew <- ifelse(covars == i, 1, 0)
    resFn <- linregfn(data, as.factor(cnew))
    pvals <- c(pvals, resFn$pvals['covars'])
  }
  allp <- linregfn(data, as.factor(covars))$pvals['covars']
  minp <- which(unlist(pvals) == min(unlist(pvals))); 
  return(list(minp = minp, pval = pvals[minp], 
              lev = l[minp], allp = allp))
}