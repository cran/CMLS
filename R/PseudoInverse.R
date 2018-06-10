PseudoInverse <-
  function(x, symmetric, tol = NULL){
    
    x <- as.matrix(x)
    if(missing(symmetric)){
      symmetric <- isSymmetric(x)
    } else {
      symmetric <- as.logical(symmetric[1])
      if(symmetric && nrow(x) != ncol(x)) {
        warning("Input 'x' is not a square matrix, but input 'symmetric' is TRUE.\n  Setting symmetric = FALSE.")
        symmetric <- FALSE
      }
    }
    
    if(is.null(tol)){
      tol <- .Machine$double.eps * max(dim(x))
    }
    
    if(symmetric){
      xeig <- eigen(x, symmetric = TRUE)
      nze <- sum(xeig$values > (tol * xeig$values[1]))
      if(nze > 1L){
        xinv <- xeig$vectors[,1:nze] %*% diag(1 / xeig$values[1:nze]) %*% t(xeig$vectors[,1:nze])
      } else {
        xinv <- tcrossprod(xeig$vectors[,1]) / xeig$values[1]
      }
    } else {
      xsvd <- svd(x)
      nze <- sum(xsvd$d > (tol * xsvd$d[1]))
      if(nze > 1L){
        xinv <- xsvd$v[,1:nze] %*% diag(1 / xsvd$d[1:nze]) %*% t(xsvd$u[,1:nze])
      } else {
        xinv <- tcrossprod(xsvd$v[,1], xsvd$u[,1]) / xsvd$d[1]
      }
    }
    xinv
    
  }