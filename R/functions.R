# Works for x in 1 dimension
WilksLambda <- function(Y, x, shrink = FALSE) {
  rho <- NULL
  N <- nrow(Y)
  p <- ncol(Y)
  bar.Y <- colMeans(Y)
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

  temp <- eigen(Vr, symmetric=TRUE)
  Ur <- temp$vectors
  diagD <- temp$values
  value <- 1/sqrt(diagD)
  root.Vr <- Ur %*% diag(value) %*% t(Ur)
  m <- root.Vr %*% Vg %*% root.Vr
  temp1 <- eigen(m, symmetric=TRUE)
  PCEV <- root.Vr %*% temp1$vectors
  d <- temp1$values
  wilks.lambda <- (N-p-1)/p * d[1]
  #wilks.lambda = d[1]
  df1 <- p
  df2 <- N-p-1
  p.value <- pf(wilks.lambda, df1, df2, lower.tail = FALSE)

  return(list("environment" = Vr,
              "genetic" = Vg,
              "PCEV"=PCEV,
              "root.Vr"=root.Vr,
              "values"=d,
              "p.value"=p.value,
              "rho"=rho))
}


# The following function computes the PCEV in two stages.
# The index parameter indicates to which cluster each responses belong.
# Type = "all" or "block"
# Works only for one covariate x
# Shrinkage = TRUE if we want to use shrinkage for the environmental variance matrix

PCEV <- function(Y, x, index, type = "all", shrink = FALSE) {
  if (type == "block") {
    if (ncol(Y) != length(index)) stop("index should have length equal to ncol(Y)")

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
    pvalue <- result$p.value
    rho <- result$rho
  }
  if (type == "all") {
    result <- WilksLambda(Y, x, shrink)
    weights <- result$PCEV[,1]
    Y.PCH <- Y %*% weights
    pvalue <- result$p.value
    rho <- result$rho
  }
  return(list("PCH.values" = Y.PCH,
              "weights" = weights,
              "Pval"= pvalue,
              "rho"= rho))
}
