cv.cmls <- 
  function(X, Y, nfolds = 2, foldid = NULL, parameters = NULL,
           const = "uncons", df = 10, degree = 3, intercept = TRUE,
           mse = TRUE, parallel = FALSE, cl = NULL, verbose = TRUE, ...){
    # k-fold Cross-Validation for cmls.R
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: June 3, 2018
    
    # check inputs
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(Y)
    m <- ncol(Y)
    p <- ncol(X)
    if(nrow(X) != n) stop("Need nrow(Y) == nrow(X).")
    
    # get folds
    if(is.null(foldid)){
      nfolds <- as.integer(nfolds[1])
      if(nfolds < 1L) stop("Input 'nfolds' must be a positive integer.")
      foldid <- NULL
      if(nfolds > 1L) foldid <- sample(rep(1:nfolds, length.out = n))
    } else {
      if(length(foldid) != n) stop("Input 'foldid' must be a vector satisfying: nrow(X) == length(foldid)")
      foldid <- as.factor(foldid)
      nfolds <- nlevels(foldid)
      if(nfolds < 2L) stop("Input 'foldid' must contain more than 1 unique fold label.")
      if(nfolds > n) stop("Input 'foldid' must contain n or fewer levels.")
      foldid <- as.integer(foldid)
    }
    
    # check parallel and cl
    if(parallel && !any(class(cl) == "cluster")) {
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?cv_cmls")
    }
    
    # parameters for cv tuning
    if(is.null(parameters)){
      
      # get combinations of parameters
      parameters <- expand.grid(const = const,
                                df = df,
                                degree = degree,
                                intercept = intercept,
                                cvloss = NA)
      parameters$const <- as.character(parameters$const)
      
      # remove duplicates
      const.smooth <- c("smooth", "smonon", "smoper", "smpeno", 
                        "ortsmo", "orsmpe", "monsmo", "mosmno", "unismo", "unsmno", 
                        "unsmpe", "unsmpn")
      rmix <- !(parameters$const %in% const.smooth)
      parameters$df[rmix] <- parameters$degree[rmix] <- parameters$intercept[rmix] <- NA
      parameters <- unique(parameters)
      npar <- nrow(parameters)
      rownames(parameters) <- 1:npar
      
    } else {
      
      parameters <- as.data.frame(parameters)
      npar <- nrow(parameters)
      if(is.null(parameters$const)) stop("Input 'parameters' must contain a column named 'const'.")
      if(is.null(parameters$df)) stop("Input 'parameters' must contain a column named 'df'.")
      if(is.null(parameters$degree)) stop("Input 'parameters' must contain a column named 'degree'.")
      if(is.null(parameters$intercept)) stop("Input 'parameters' must contain a column named 'intercept'.")
      if(is.null(parameters$cvloss)) parameters$cvloss <- rep(NA, npar)
      
    } # end if(is.null(parameters))
    
    # tune model
    XtX <- NULL
    if(nfolds == 1L) XtX <- crossprod(X)
    if(parallel){
      
      parList <- split(parameters, f = seq(npar))
      cvloss <- parSapply(cl = cl, X = parList, 
                          FUN = ifelse(nfolds == 1L, "gcvfun", ifelse(mse, "kcvmse", "kcvmae")), 
                          Xmat = X, Ymat = Y, nfolds = nfolds, 
                          foldid = foldid, XtX = XtX, ...)
      parameters$cvloss <- cvloss
      
    } else {
      
      if(verbose) pbar <- txtProgressBar(min = 0, max = npar, style = 3)
      
      # loop through parameters
      if(nfolds == 1L){
        
        for(jj in 1:npar) {
          parameters$cvloss[jj] <- gcvfun(x = parameters[jj,], Xmat = X, Ymat = Y,
                                          nfolds = nfolds, foldid = foldid, XtX = XtX, ...)
          if(verbose) setTxtProgressBar(pbar, jj)
        }
        
      } else {
        
        if(mse){
          for(jj in 1:npar) {
            parameters$cvloss[jj] <- kcvmse(x = parameters[jj,], Xmat = X, Ymat = Y,
                                            nfolds = nfolds, foldid = foldid, XtX = XtX, ...)
            if(verbose) setTxtProgressBar(pbar, jj)
          }
        } else {
          for(jj in 1:npar) {
            parameters$cvloss[jj] <- kcvmae(x = parameters[jj,], Xmat = X, Ymat = Y,
                                            nfolds = nfolds, foldid = foldid, XtX = XtX, ...)
            if(verbose) setTxtProgressBar(pbar, jj)
          }
        } # end if(mse)
        
      } # end if(nfolds == 1L)
      
      if(verbose) close(pbar)
      
    } # end if(parallel)
    
    ix <- sort(parameters$cvloss, index.return = TRUE)$ix
    return(list(best.parameters = parameters[ix[1],],
                top5.parameters = parameters[ix[1:min(length(ix),5)],],
                full.parameters = parameters,
                type = ifelse(nfolds == 1L, "GCV",
                              ifelse(mse, "KCV-MSE", "KCV-MAE"))))
    
} # end cv.cmls.R