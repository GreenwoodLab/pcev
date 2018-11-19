############################################
# For helper functions that are not exported
############################################
computeVIMP <- function(pcevObj, list, signed=FALSE) {
  
  VIMP <- cor(pcevObj$Y, list$PCEV)[,1]
  if (!signed) {
    VIMP <- abs(VIMP)
  }
  
  return(VIMP)
}

shrink_est <- function(Vr, res){
  # port of matlab code from http://www.econ.uzh.ch/faculty/wolf/publications.html#9
  # Ledoit, O. and Wolf, M. (2004).
  # Honey, I shrunk the sample covariance matrix.
  # Journal of Portfolio Management 30, Volume 4, 110-119.
  p <- ncol(res); n <- nrow(res)
  
  # Compute sample covariance matrix using the de-meaned returns
  sample <- Vr/n
  
  # Compute prior
  var <- matrix(diag(sample), ncol = 1)
  sqrtvar <- sqrt(var)
  tmpMat <- matrix(rep(sqrtvar, p), nrow = p)
  rBar <- (sum(sum(sample / (tmpMat * t(tmpMat)))) - p) / (p * (p - 1))
  prior <- rBar * tmpMat * t(tmpMat)
  diag(prior) <- var
  
  # What is called pi-hat
  y <- res^2
  phiMat <- crossprod(y) / n - 2 * crossprod(res) * sample / n + sample^2
  phi <- sum(phiMat)
  
  # What is called rho-hat
  term1 <- crossprod(res^3, res) / n
  help <- crossprod(res)/n
  helpDiag <- matrix(diag(help), ncol = 1)
  term2 <- matrix(rep(helpDiag, p), ncol = p, byrow = FALSE) * sample
  term3 <- help * matrix(rep(var, p), ncol = p, byrow = FALSE)
  term4 <- matrix(rep(var, p), ncol = p, byrow = FALSE) * sample
  thetaMat <- term1 - term2 - term3 + term4
  diag(thetaMat) <- 0
  rho <- sum(diag(phiMat)) + rBar * sum(sum(tcrossprod(1 / sqrtvar, sqrtvar) * thetaMat))
  
  # What is called gamma-hat
  gamma <- norm(sample - prior, "F")^2
  
  # Compute shrinkage constant
  kappa <- (phi - rho) / gamma
  shrinkage <- max(0, min(1, kappa / n))
  
  # Compute the estimator
  sigma <- shrinkage * prior + (1 - shrinkage) * sample
  sigma <- n * sigma
  out <- list(cov = sigma, 
              rho = shrinkage)
  return(out)
}

########################################################
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

getFitted <- function(x, y) {
  if (FALSE) {
    return(lm.fit(x, y)$fitted.values)
  } else {
    coef <- matrix(NA, nrow = ncol(y),
                   ncol = ncol(x))
    for (i in seq_len(ncol(y))) {
      # X_des <- model.matrix(y[,i] ~ x)
      model <- lm.fit(x, na.omit(y[,i]))
      coef[i,] <- model$coefficients
    }
    return(tcrossprod(x, coef))
  }
}
