# Influence measure for observations----
# From Pendergast 2011

#' Influence measure on the PCEV
#' 
#' This function computes a measure of the influence of each observation on the PCEV.
#' 
#' This measure is computed following equation 14 in Pendergast (2011).
#' 
#' @param model an object of class \code{Pcev}. 
#' @param ... further arguments passed to or from other methods.
#' @return Influence measures
#' @export
#' @importFrom stats influence
#' @rdname influence
influence.Pcev <- function(model, ...) {
  if(ncol(model$pcevObj$X) < 3) {
    stop("Need at least two covariates to be able to compute influence measures.")
  }
  if(inherits(model$pcevObj, "PcevBlock")) {
    stop("Only available for classical and singular estimation")
  }
  out <- SCIA(vectors = model$genEigSol$sampleComp, 
              values = model$genEigSol$sampleVal)
  names(out) <- rownames(model$pcevObj$Y)
  
  return(out)
}

SCIA <- function(vectors, values) {
  vectors_S <- vectors[,1,drop=TRUE]
  vectors_notS <- vectors[,-1,drop=FALSE]
  values_S <- values[1]
  values_notS <- values[-1]
  
  multiplier <- values_notS/(values_S * (values_S - values_notS)^2)
  multiplier <- matrix(multiplier, nrow=length(vectors_S),
                       ncol=length(multiplier), byrow = TRUE)
  numerator <- vectors_notS^2 * vectors_S^2
  
  influence <- rowSums(numerator * multiplier)  
  
  return(influence)
}