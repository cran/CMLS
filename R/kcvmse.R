kcvmse <- 
  function(x, Xmat, Ymat, nfolds, foldid, ...){
    cverr <- rep(0, nfolds)
    for(kk in 1:nfolds){
      indx <- which(foldid == kk)
      Bhat <- cmls(X = Xmat[-indx,,drop=FALSE], Y = Ymat[-indx,,drop=FALSE], 
                   const = x$const, df = x$df, 
                   degree = x$degree, intercept = x$intercept, ...)
      cverr[kk] <- mean((Ymat[indx,,drop=FALSE] - Xmat[indx,,drop=FALSE] %*% Bhat)^2)
    } # end for(kk in 1:nfolds)
    mean(cverr)
  } # end kcvmse.R