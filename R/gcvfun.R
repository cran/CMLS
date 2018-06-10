gcvfun <-
  function(x, Xmat, Ymat, nfolds, foldid, ...){
    Bhat <- cmls(X = Xmat, Y = Ymat, const = x$const, 
                 df = x$df, degree = x$degree, 
                 intercept = x$intercept, ...)
    edf <- sum(attr(Bhat, "df"))
    nm <- prod(dim(Ymat))
    top <- mean((Ymat - Xmat %*% Bhat)^2)
    bot <- (1 - edf / nm)^2
    top / bot
  } # end gcvfun.R