IsplineBasis <-
  function(x, df = NULL, knots = NULL, degree = 2, intercept = FALSE, 
           Boundary.knots = range(x)){
    
    # is intercept needed?
    if(intercept){
      if(!is.null(df)) df <- df - 1L
      X <- cbind(1, IsplineBasis(x, df = df, knots = knots, degree = degree,
                                 Boundary.knots = Boundary.knots))
      attr(X, "intercept") <- TRUE
      return(X)
    }
    
    # get needed M-spline basis
    if(!is.null(df)) df <- df + 1L
    degree <- degree + 1L
    Xmat <- MsplineBasis(x = x, df = df, knots = knots, degree = degree,
                         intercept = TRUE, Boundary.knots = Boundary.knots)
    
    # form knots and get info
    Aknots <- sort(c(rep(attr(Xmat, "Boundary.knots"), degree+1L), attr(Xmat ,"knots")))
    nx <- length(x)
    nk <- length(Aknots)
    df <- ncol(Xmat)
    
    # make design matrix
    X <- matrix(0, nrow=nx, ncol=df-1L)
    for(h in 1:nx){
      j <- sum(x[h] >= Aknots)
      for(i in 2:df){
        if(i > j){
          X[h,i-1L] <- 0
        } else if(i < (j - degree + 1)){
          X[h,i-1L] <- 1 
        } else {
          for(m in i:j){
            X[h,i-1L] <- X[h,i-1L] + (Aknots[m+degree+1] - Aknots[m]) * Xmat[h,m] / (degree + 1)
          }
        }
      } # end for(i in 2:df)
    } # end for(h in 1:nx)
    
    dimnames(X) <- list(1:nx, 1:ncol(X))
    a <- list(degree = degree, knots = attr(Xmat ,"knots"), 
              intercept = intercept, Boundary.knots = Boundary.knots)
    attributes(X) <- c(attributes(X), a)
    return(X)
    
  }