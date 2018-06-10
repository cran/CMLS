rsinvsqrt <- function(X, method = c("chol", "eigen")){
  # ridge stabilized inverse square root
  X <- as.matrix(X)
  pts <- nrow(X)
  if(ncol(X) != pts) stop("Input 'X' must be a symmetric matrix.")
  eps <- .Machine$double.eps * pts
  if(method[1] == "eigen"){
    Xeig <- eigen(X, symmetric = TRUE)
    nvals <- sum(Xeig$values > (eps * Xeig$values[1]))
    if(nvals < pts) Xeig$values <- Xeig$values + (eps * Xeig$values[1] - Xeig$values[pts])
    for(k in 1:pts) Xeig$vectors[,k] <- Xeig$vectors[,k] / sqrt(Xeig$values[k])
    return(Xeig$vectors)
  } else {
    Xeig <- eigen(X, symmetric = TRUE, only.values = TRUE)
    nvals <- sum(Xeig$values > (eps * Xeig$values[1]))
    if(nvals < pts) X <- X + diag(rep(eps * Xeig$values[1] - Xeig$values[pts], pts))
    return(solve(chol(X)))
  }
  
} # end Rinv.R