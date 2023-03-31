cmls <- 
  function(X, Y, const = "uncons", struc = NULL, 
           z = NULL, df = 10, degree = 3, intercept = TRUE, 
           backfit = FALSE, maxit = 1e3, eps = 1e-10, 
           del = 1e-6, XtX = NULL, mode.range = NULL){
    # Constrained Multivariate Least Squares
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: March 31, 2021
    
    # Finds the B that minimizes:  sum( ( Y - X %*% B )^2 )
    # subject to specified constraints on rows of B
    
    
    ######******######   CONSTRAINT OPTIONS   ######******######
    
    # __________________________________________________________
    # uncons = unconstrained
    # nonneg = non-negative
    # period = periodicity
    # pernon = periodicity and non-negativity
    # __________________________________________________________
    # smooth = smoothness
    # smonon = smoothness and non-negative
    # smoper = smoothness and periodic
    # smpeno = smoothness and periodic and non-negative
    # __________________________________________________________
    # orthog = orthogonal
    # ortnon = orthogonal and non-negative
    # ortsmo = orthogonal and smooth
    # orsmpe = orthogonal and smooth and periodic
    # __________________________________________________________
    # moninc = monotonic increasing
    # monnon = monotonic increasing and non-negative
    # monsmo = monotonic increasing and smooth
    # mosmno = monotonic increasing and smooth and non-negative
    # __________________________________________________________
    # unimod = unimodal
    # uninon = unimodal and non-negative
    # uniper = unimodal and periodic
    # unpeno = unimodal and periodic and non-negative
    # __________________________________________________________
    # unismo = unimodal and smooth
    # unsmno = unimodal and smooth and and non-negative
    # unsmpe = unimodal and smooth and periodic
    # unsmpn = unimodal and smooth and periodic and non-negative
    
    
    
    ######******######   INITITAL CHECKS   ######******######
    
    catchError <- TRUE
    
    # check inputs
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(Y)
    m <- ncol(Y)
    p <- ncol(X)
    if(nrow(X) != n) stop("Need nrow(Y) == nrow(X).")
    
    # check const
    if(is.na(const) | is.null(const) | missing(const)) const <- "uncons"
    const <- as.character(const[1])
    const.types <- c("uncons", "nonneg", "period", "pernon",
                     "smooth", "smonon", "smoper", "smpeno",
                     "orthog", "ortnon", "ortsmo", "orsmpe",
                     "moninc", "monnon", "monsmo", "mosmno",
                     "unimod", "uninon", "uniper", "unpeno", 
                     "unismo", "unsmno", "unsmpe", "unsmpn")
    if(!any(const == const.types)) stop("Invalid 'const' input.")
    
    # check struc
    if(!is.null(struc)){
      struc <- as.matrix(struc)
      if(nrow(struc) != p | ncol(struc) != m) stop("Input 'struc' must satisfy:  nrow(struc) == ncol(X)  &&  ncol(struc) == ncol(Y)")
      schck <- max(abs( ((struc == TRUE) + (struc == FALSE)) - 1 ))
      if(schck > 0) stop("Input 'struc' must be a matrix of logicals (TRUE/FALSE).")
      if(any(const == c("orthog", "ortnon", "ortsmo", "orsmpe"))){
        if(p > 1){
          StS <- tcrossprod(struc)
          maxLT <- max(abs(StS[lower.tri(StS)]))
          if(maxLT > 0) stop("When using orthogonality constraints, the input 'struc' must\n  contain orthogonal rows (i.e., one or less TRUE in each column).")
        }
      }
    }
    
    # check XtX
    if(is.null(XtX)){
      XtX <- crossprod(X)
    } else {
      if(nrow(XtX) != p | ncol(XtX) != p) stop("Inputs 'X' and 'XtX' are incompatible:  XtX = crossprod(X)")
    }
    
    # check inputs relevant to smoothness updates
    const.smooth <- c("smooth", "smonon", "smoper", "smpeno",
                      "ortsmo", "orsmpe", "monsmo", "mosmno", 
                      "unismo", "unsmno", "unsmpe", "unsmpn")
    if(any(const == const.smooth)){
      
      # check z
      if(is.null(z)){
        z <- seq(0, 1, length.out = m)
      } else {
        z <- as.numeric(z)
        if(length(z) != m) stop("Inputs 'Y' and 'z' are incompatible:\nNeed ncol(Y) == length(z)")
      }
      
      # check df
      df <- as.integer(df[1])
      
      # check degree
      degree <- as.integer(degree[1])
      if(degree < 1L) stop("Input 'degree' should be greater than 0.") 
      
      # check intercept
      intercept <- intercept[1]
      if(!is.logical(intercept)) stop("Input 'intercept' must be a logical (TRUE/FALSE).")
      
      # check df (relative to degree and intercept)
      mindf <- degree + 1 - !intercept
      if(any(const == c("monsmo", "mosmno")) && intercept) mindf <- mindf + 1L
      if(df < mindf){
        df <- mindf
        warning("Input 'df' is too small!\n  Resetting 'df' to minimum possible value = ", mindf)
      } 
      
    } # end if(any(const == const.smooth))
    
    # check inputs relevant to unimodality updates
    const.unimod <- c("unimod", "uninon", "uniper", "unpeno", 
                      "unismo", "unsmno", "unsmpe", "unsmpn")
    if(any(const == const.unimod) & !is.null(mode.range)){
      mode.range <- as.matrix(mode.range)
      if(nrow(mode.range) != 2L | ncol(mode.range) != p) stop("Input 'mode.range' must be a 2 x p matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for each row of B.")
      mode.range <- matrix(as.integer(mode.range), nrow = 2L, ncol = p)
      if(any(mode.range[1,] < 1L)) stop("First row of 'mode.range' must contain integers greater than or equal to one.")
      if(any(mode.range[2,] > m)) stop("Second row of 'mode.range' must contain integers less than or equal to m.")
      if(any((mode.range[2,] - mode.range[1,]) < 0)) stop("Input 'mode.range' must satisfy:  mode.range[1,j] <= mode.range[2,j]")
    } # end if(any(const == const.unimod))
    
    
    ######******######   UNCONSTRAINED   ######******######
    
    # unconstrained (uncons): unstructured or structured
    if(const == "uncons"){
      if(is.null(struc)){
        B <- PseudoInverse(XtX, symmetric = TRUE) %*% crossprod(X, Y)
        Bdf <- rep(ncol(B), nrow(B))
      } else {
        B <- matrix(0, nrow = p, ncol = m)
        XtY <- crossprod(X, Y)
        Rinv <- rsinvsqrt(XtX)
        for(kk in 1:m){
          ntrue <- sum(struc[,kk])
          if(ntrue > 0L){
            meq <- p - ntrue
            Amat <- matrix(0, nrow = p, ncol = meq)
            Amat[!struc[,kk],] <- diag(meq)
            B[,kk] <- solveQP(Dmat = Rinv, dvec = XtY[,kk], Amat = Amat, meq = meq,
                              factorized = TRUE, catchError = catchError)$solution
          }
        }
        Bdf <- rowSums(struc)
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(const == "uncons")
    
    
    
    ######******######   NON-NEGATIVITY   ######******######
    
    # non-negative (nonneg): unstructured or structured
    if(const == "nonneg"){
      B <- matrix(0, nrow = p, ncol = m)
      Bdf <- matrix(1, nrow = p, ncol = m)
      XtY <- crossprod(X, Y)
      Rinv <- rsinvsqrt(XtX)
      Amat <- diag(p)
      if(is.null(struc)){
        for(kk in 1:m) {
          Bfit <- solveQP(Dmat = Rinv, dvec = XtY[,kk], Amat = Amat, 
                          factorized = TRUE, catchError = catchError)
          Bdf[Bfit$iact,kk] <- 0
          B[,kk] <- Bfit$solution
        }
      } else {
        for(kk in 1:m) {
          sidx <- sort(struc[,kk], index.return = TRUE)$ix
          cidx <- (1:p)[sidx]
          Amati <- Amat[, sidx, drop = FALSE]
          meqi <- sum(!struc[,kk])
          Bfit <- solveQP(Dmat = Rinv, dvec = XtY[,kk], Amat = Amati, meq = meqi, 
                          factorized = TRUE, catchError = catchError)
          Bdf[cidx[Bfit$iact],kk] <- 0
          B[,kk] <- Bfit$solution
        }
      } # end if(is.null(struc))
      Bdf <- rowSums(Bdf)
      attr(B, "df") <- Bdf
      return(B)
    } # end if(const == "nonneg")
    
    
    
    ######******######   PERIODICITY   ######******######
    
    # period or pernon: unstructured or structured
    if(any(const == c("period", "pernon"))){
      nonneg <- ifelse(const == "pernon", TRUE, FALSE)
      A <- matrix(rep(c(-1, 0, 1), c(1, m-2, 1)), nrow = m, ncol = 1)
      if(is.null(struc)){
        if(nonneg) A <- cbind(A, diag(m))
        B <- t(mlsei(X = X, Y = Y, A = A, meq = 1, 
                     backfit = backfit, maxit = maxit, eps = eps, del = del, 
                     XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } else {
        Amat <- vector("list", p)
        if(nonneg) Aeye <- diag(m)
        meq <- rep(1, p)
        for(jj in 1:p){
          meq[jj] <- 1L + sum(!struc[jj,])
          if(nonneg){
            Amat[[jj]] <- cbind(A, Aeye[,sort(struc[jj,], index.return = TRUE)$ix])
          } else {
            if(meq[jj] > 1L){
              Anew <- matrix(0, nrow = m, ncol = meq[jj] - 1L)
              Anew[!struc[jj,],] <- diag(meq[jj] - 1L)
              Amat[[jj]] <- cbind(A, Anew)
            } else {
              Amat[[jj]] <- A
            }
          }
        } # end for(jj in 1:p)
        B <- t(mlsei(X = X, Y = Y, A = Amat, meq = meq, 
                     backfit = backfit, maxit = maxit, eps = eps, del = del, 
                     XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(any(const == c("period", "pernon")))
    
    
    
    ######******######   SMOOTHNESS   ######******######
    
    # smooth or smoper: unstructured or structured
    if(any(const == c("smooth", "smoper"))){
      periodic <- ifelse(const == "smoper", TRUE, FALSE)
      #zindx <- seq(0, 1, length.out = m)
      if(is.null(struc)){
        #Z <- MsplineBasis(x = zindx, df = df, degree = degree, intercept = intercept)
        Z <- MsplineBasis(x = z, df = df, degree = degree, intercept = intercept)
        Z <- scale(Z, center = FALSE)
        if(periodic){
          A <- t(Z[m,,drop=FALSE] - Z[1,,drop=FALSE])
          B <- mlsei(X = X, Y = Y, Z = Z, A = A, meq = 1,
                     backfit = backfit, maxit = maxit, eps = eps, del = del,
                     XtX = XtX, catchError = catchError)
          Bdf <- attr(B, "df")
          B <- t(Z %*% B)
        } else {
          B <- PseudoInverse(X) %*% tcrossprod(Y, PseudoInverse(Z))
          B <- tcrossprod(B, Z)
          Bdf <- rep(df, p)
        } # end if(periodic)
      } else {
        Amat <- Zmat <- vector("list", p)
        meq <- rep(1, p)
        for(jj in 1:p){
          #zjj <- zindx[struc[jj,]]
          zjj <- z[struc[jj,]]
          ndf <- max(round(df * length(zjj) / m), mindf)
          Z <- matrix(0, nrow = m, ncol = ndf)
          Z[struc[jj,],] <- MsplineBasis(x = zjj, df = ndf, degree = degree, intercept = intercept)
          Zmat[[jj]] <- scale(Z, center = FALSE)
          if(periodic) Amat[[jj]] <- t(Zmat[[jj]][m,,drop=FALSE] - Zmat[[jj]][1,,drop=FALSE])
        } # end for(jj in 1:p)
        if(periodic){
          BB <- mlsei(X = X, Y = Y, Z = Zmat, A = Amat, meq = meq, 
                      backfit = backfit, maxit = maxit, eps = eps, del = del, 
                      XtX = XtX, simplify = FALSE, catchError = catchError)
        } else {
          BB <- mlsei(X = X, Y = Y, Z = Zmat, backfit = backfit,
                      maxit = maxit, eps = eps, del = del, 
                      XtX = XtX, simplify = FALSE, catchError = catchError)
        }
        Bdf <- attr(BB, "df")
        B <- matrix(0, nrow = p, ncol = m)
        for(jj in 1:p) B[jj,] <- Zmat[[jj]] %*% BB[[jj]]
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(any(const == c("smooth", "smoper")))
    
    # smonon or smpeno: unstructured or structured
    if(any(const == c("smonon", "smpeno"))){
      periodic <- ifelse(const == "smpeno", TRUE, FALSE)
      #zindx <- seq(0, 1, length.out = m)
      if(is.null(struc)){
        #Z <- MsplineBasis(x = zindx, df = df, degree = degree, intercept = intercept)
        Z <- MsplineBasis(x = z, df = df, degree = degree, intercept = intercept)
        Z <- scale(Z, center = FALSE)
        A <- t(Z)
        meq <- as.integer(periodic)
        if(periodic) A <- cbind(Z[m,] - Z[1,], A)
        B <- mlsei(X = X, Y = Y, Z = Z, A = A, meq = meq,
                   backfit = backfit, maxit = maxit, eps = eps, del = del,
                   XtX = XtX, catchError = catchError)
        Bdf <- attr(B, "df")
        B <- t(Z %*% B)
      } else {
        Amat <- Zmat <- vector("list", p)
        meq <- rep(as.integer(periodic), p)
        for(jj in 1:p){
          #zjj <- zindx[struc[jj,]]
          zjj <- z[struc[jj,]]
          ndf <- max(round(df * length(zjj) / m), mindf)
          Z <- matrix(0, nrow = m, ncol = ndf)
          Z[struc[jj,],] <- MsplineBasis(x = zjj, df = ndf, degree = degree, intercept = intercept)
          Zmat[[jj]] <- scale(Z, center = FALSE)
          if(periodic){
            Amat[[jj]] <- cbind(Zmat[[jj]][m,] - Zmat[[jj]][1,], t(Zmat[[jj]]))
          } else {
            Amat[[jj]] <- t(Zmat[[jj]])
          }
        } # end for(jj in 1:p)
        BB <- mlsei(X = X, Y = Y, Z = Zmat, A = Amat, meq = meq, 
                    backfit = backfit, maxit = maxit, eps = eps, del = del, 
                    XtX = XtX, simplify = FALSE, catchError = catchError)
        Bdf <- attr(BB, "df")
        B <- matrix(0, nrow = p, ncol = m)
        for(jj in 1:p) B[jj,] <- Zmat[[jj]] %*% BB[[jj]]
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(any(const == c("smonon", "smpeno")))
    
    
    
    ######******######   ORTHOGONALITY   ######******######
    
    # orthog: unstructured and structured
    if(const == "orthog"){
      if(is.null(struc)){
        XtYsvd <- svd(crossprod(X, Y))
        B <- tcrossprod(XtYsvd$u, XtYsvd$v)
        Bdf <- rep(m - (p+1)/2, p)
      } else {
        B <- matrix(0, nrow = p, ncol = m)
        dXtX <- diag(XtX)
        for(kk in 1:m){
          trueID <- which(struc[,kk])
          if(length(trueID) > 0) B[trueID , kk] <- crossprod(X[,trueID], Y[,kk]) / dXtX[trueID]
        }
        Bdf <- rowSums(struc)
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(const == "orthog")
    
    # ortnon: unstructured and structured
    if(const == "ortnon"){
      if(is.null(struc)){
        
        B <- matrix(0, nrow = p, ncol = m)
        YtY <- colSums(Y^2)
        XtY <- crossprod(X, Y)
        dXtX <- diag(XtX)
        zeros <- rep(0, p)
        for(kk in 1:m){
          coefs <- pmax(XtY[,kk] / dXtX, zeros)
          if(any(coefs > 0)){
            ssevals <- rep(YtY[kk], p) - 2 * XtY[,kk] * coefs + dXtX * coefs^2
            minid <- which.min(ssevals)
            B[minid,kk] <- coefs[minid]
          }
        }
        Bdf <- rowSums(B > 0)
      } else {
        B <- matrix(0, nrow = p, ncol = m)
        YtY <- colSums(Y^2)
        XtY <- crossprod(X, Y)
        dXtX <- diag(XtX)
        zeros <- rep(0, p)
        for(kk in 1:m){
          ix <- which(struc[,kk])
          B[ix,kk] <- max(XtY[ix,kk] / dXtX[ix], 0)
        }
        Bdf <- rowSums(B > 0)
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(const == "ortnon")
    
    # ortsmo or orsmpe: unstructured and structured
    if(any(const == c("ortsmo","orsmpe"))){
      periodic <- ifelse(const == "orsmpe", TRUE, FALSE)
      #zindx <- seq(0, 1, length.out = m)
      if(is.null(struc)){
        # Z <- MsplineBasis(x = zindx, df = df, degree = degree, 
        #                   intercept = intercept, periodic = periodic)
        Z <- MsplineBasis(x = z, df = df, degree = degree, 
                          intercept = intercept, periodic = periodic)
        Z <- scale(Z, center = FALSE)
        Zsvd <- svd(Z)
        Tsvd <- svd(crossprod(Zsvd$u, crossprod(Y, X)))
        B <- Zsvd$v %*% diag(1/Zsvd$d) %*% tcrossprod(Tsvd$u, Tsvd$v)
        Bdf <- rep(nrow(B) - (ncol(B)+1)/2, ncol(B))
        B <- t(Z %*% B)
      } else {
        if(periodic){
          B <- cmls(X = X, Y = Y, const = "smoper", struc = struc,
                    backfit = backfit, maxit = maxit, eps = eps, del = del,
                    df = df, degree = degree, intercept = intercept, XtX = XtX)
          Bdf <- attr(B, "df")
        } else {
          B <- cmls(X = X, Y = Y, const = "smooth", struc = struc,
                    backfit = backfit, maxit = maxit, eps = eps, del = del,
                    df = df, degree = degree, intercept = intercept, XtX = XtX)
          Bdf <- attr(B, "df")
        }
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(any(const == c("ortsmo","orsmpe")))
    
    
    
    ######******######   MONOTONICITY   ######******######
    
    # moninc and monnon: unstructured or structured
    if(any(const == c("moninc", "monnon"))){
      A <- diag(m)
      A <- A[,2:m] - A[,1:(m-1)]
      if(const == "monnon") A <- cbind(rep(c(1, 0), c(1, m-1)), A)
      if(is.null(struc)){
        B <- t(mlsei(X = X, Y = Y, A = A, backfit = backfit,
                     maxit = maxit, eps = eps, del = del, XtX = XtX, 
                     catchError = catchError))
        Bdf <- attr(B, "df")
      } else {
        Amat <- vector("list", p)
        meq <- rep(0, p)
        for(jj in 1:p) {
          meq[jj] <- sum(!struc[jj,])
          if(meq[jj] > 0){
            Anew <- matrix(0, nrow = m, ncol = meq[jj])
            Anew[!struc[jj,],] <- diag(meq[jj])
            Amat[[jj]] <- cbind(Anew, A)
          } else {
            Amat[[jj]] <- A
          }
        } # end for(jj in 1:p)
        B <- t(mlsei(X = X, Y = Y, A = Amat, meq = meq, backfit = backfit,
                     maxit = maxit, eps = eps, del = del, XtX = XtX, 
                     catchError = catchError))
        Bdf <- attr(B, "df")
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(any(const == c("moninc", "monnon")))
    
    # monsmo and mosmno: unstructured or structured
    if(any(const == c("monsmo", "mosmno"))){
      nonneg <- ifelse(const == "mosmno", TRUE, FALSE)
      #zindx <- seq(0, 1, length.out = m)
      if(is.null(struc)){
        #Z <- IsplineBasis(x = zindx, df = df, degree = degree, intercept = intercept)
        Z <- IsplineBasis(x = z, df = df, degree = degree, intercept = intercept)
        A <- t(Z[2:m,] - Z[1:(m-1),])
        if(nonneg) A <- cbind(Z[1,], A)
        B <- mlsei(X = X, Y = Y, Z = Z, A = A, backfit = backfit,
                   maxit = maxit, eps = eps, del = del, XtX = XtX, 
                   catchError = catchError)
        Bdf <- attr(B, "df")
        B <- t(Z %*% B)
      } else {
        Amat <- Zmat <- vector("list", p)
        for(jj in 1:p){
          #zjj <- zindx[struc[jj,]]
          zjj <- z[struc[jj,]]
          ndf <- max(round(df * length(zjj) / m), mindf)
          Z <- matrix(0, nrow = m, ncol = ndf)
          Z[struc[jj,],] <- IsplineBasis(x = zjj, df = ndf, degree = degree, intercept = intercept)
          Zmat[[jj]] <- scale(Z, center = FALSE)
          if(nonneg){
            Amat[[jj]] <- cbind(Zmat[[jj]][1,], t(Zmat[[jj]][2:m,] - Zmat[[jj]][1:(m-1),]))
          } else {
            Amat[[jj]] <- t(Zmat[[jj]][2:m,] - Zmat[[jj]][1:(m-1),])
          }
        } # end for(jj in 1:p)
        BB <- mlsei(X = X, Y = Y, Z = Zmat, A = Amat, backfit = backfit,
                    maxit = maxit, eps = eps, del = del, XtX = XtX,
                    simplify = FALSE, catchError = catchError)
        Bdf <- attr(BB, "df")
        B <- matrix(0, nrow = p, ncol = m)
        for(jj in 1:p) B[jj,] <- Zmat[[jj]] %*% BB[[jj]]
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(any(const == c("monsmo", "mosmno")))
    
    
    
    ######******######   UNIMODALITY   ######******######
    
    # unimod: unstructured or structured
    if(const == "unimod"){
      if(is.null(struc)){
        B <- t(mlsun(X = X, Y = Y, mode.range = mode.range,
                     maxit = maxit, eps = eps, del = del, 
                     XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } else {
        Amat <- vector("list", p)
        meq <- rep(0, p)
        for(jj in 1:p) {
          meq[jj] <- sum(!struc[jj,])
          if(meq[jj] > 0){
            Amat[[jj]] <- matrix(0, nrow = m, ncol = meq[jj])
            Amat[[jj]][!struc[jj,],] <- diag(meq[jj])
          } else {
            Amat[[jj]] <- matrix(0, nrow = m, ncol = 1)
          }
        }
        if(is.null(mode.range)) mode.range <- apply(struc, 1, function(x) range(which(x)))
        B <- t(mlsun(X = X, Y = Y, A = Amat, meq = meq, 
                     mode.range = mode.range, maxit = maxit, eps = eps, 
                     del = del, XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(const == "unimod")
    
    # uninon: unstructured or structured
    if(const == "uninon"){
      if(is.null(struc)){
        A <- cbind(rep(c(1, 0), c(1, m - 1)), rep(c(0, 1), c(m - 1, 1)))
        B <- t(mlsun(X = X, Y = Y, A = A, mode.range = mode.range,
                     maxit = maxit, eps = eps, del = del, 
                     XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } else {
        Annc <- cbind(rep(c(1, 0), c(1, m - 1)), rep(c(0, 1), c(m - 1, 1)))
        Amat <- vector("list", p)
        meq <- rep(0, p)
        for(jj in 1:p) {
          meq[jj] <- sum(!struc[jj,])
          if(meq[jj] > 0){
            Ajj <- matrix(0, nrow = m, ncol = meq[jj])
            Ajj[!struc[jj,],] <- diag(meq[jj])
            Amat[[jj]] <- cbind(Ajj, Annc)
          } else {
            Amat[[jj]] <- Annc
          }
        }
        if(is.null(mode.range)) mode.range <- apply(struc, 1, function(x) range(which(x)))
        B <- t(mlsun(X = X, Y = Y, A = Amat, meq = meq, 
                     mode.range = mode.range, maxit = maxit, eps = eps, 
                     del = del, XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(const == "uninon")
    
    # uniper: unstructured or structured
    if(const == "uniper"){
      A <- matrix(rep(c(-1, 0, 1), c(1, m - 2, 1)), nrow = m, ncol = 1)
      if(is.null(struc)){
        B <- t(mlsun(X = X, Y = Y, A = A, meq = 1,
                     mode.range = mode.range, maxit = maxit, 
                     eps = eps, del = del, XtX = XtX, 
                     catchError = catchError))
        Bdf <- attr(B, "df")
      } else {
        Amat <- vector("list", p)
        meq <- rep(1, p)
        for(jj in 1:p) {
          meq[jj] <- sum(!struc[jj,]) + 1L
          if(meq[jj] > 1L){
            Ajj <- matrix(0, nrow = m, ncol = meq[jj] - 1L)
            Ajj[!struc[jj,],] <- diag(meq[jj] - 1L)
            Amat[[jj]] <- cbind(A, Ajj)
          } else {
            Amat[[jj]] <- A
          }
        }
        if(is.null(mode.range)) mode.range <- apply(struc, 1, function(x) range(which(x)))
        B <- t(mlsun(X = X, Y = Y, A = Amat, meq = meq,
                     mode.range = mode.range, maxit = maxit, eps = eps, 
                     del = del, XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(const == "uniper")
    
    # unpeno: unstructured or structured
    if(const == "unpeno"){
      Aeq <- matrix(rep(c(-1, 0, 1), c(1, m - 2, 1)), nrow = m, ncol = 1)
      Ann <- cbind(rep(c(1, 0), c(1, m - 1)), rep(c(0, 1), c(m - 1, 1)))
      if(is.null(struc)){
        B <- t(mlsun(X = X, Y = Y, A = cbind(Aeq, Ann), meq = 1,
                     mode.range = mode.range, maxit = maxit, eps = eps, 
                     del = del, XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } else {
        Amat <- vector("list", p)
        meq <- rep(1, p)
        for(jj in 1:p) {
          meq[jj] <- sum(!struc[jj,]) + 1L
          if(meq[jj] > 1L){
            Ajj <- matrix(0, nrow = m, ncol = meq[jj] - 1L)
            Ajj[!struc[jj,],] <- diag(meq[jj] - 1L)
            Amat[[jj]] <- cbind(Aeq, Ajj, Ann)
          } else {
            Amat[[jj]] <- cbind(Aeq, Ann)
          }
        }
        if(is.null(mode.range)) mode.range <- apply(struc, 1, function(x) range(which(x)))
        B <- t(mlsun(X = X, Y = Y, A = Amat, meq = meq,
                     mode.range = mode.range, maxit = maxit, eps = eps, 
                     del = del, XtX = XtX, catchError = catchError))
        Bdf <- attr(B, "df")
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(const == "unpeno")
    
    
    
    ######******######   UNIMODALITY AND SMOOTHNESS   ######******######
    
    # unismo and unsmpe: unstructured or structured
    if(any(const == c("unismo", "unsmpe"))){
      periodic <- ifelse(const == "unsmpe", TRUE, FALSE)
      #zindx <- seq(0, 1, length.out = m)
      if(is.null(struc)){
        #Z <- MsplineBasis(x = zindx, df = df, degree = degree, intercept = intercept)
        Z <- MsplineBasis(x = z, df = df, degree = degree, intercept = intercept)
        Z <- scale(Z, center = FALSE)
        if(periodic){
          A <- t(Z[m,,drop=FALSE] - Z[1,,drop=FALSE])
          B <- mlsun(X = X, Y = Y, Z = Z, A = A, meq = 1,
                     mode.range = mode.range, maxit = maxit, eps = eps,
                     del = del, XtX = XtX, catchError = catchError)
        } else {
          B <- mlsun(X = X, Y = Y, Z = Z, mode.range = mode.range,
                     maxit = maxit, eps = eps, del = del,
                     XtX = XtX, catchError = catchError)
        }
        Bdf <- attr(B, "df")
        B <- t(Z %*% B)
      } else {
        Amat <- Zmat <- vector("list", p)
        meq <- rep(1, p)
        for(jj in 1:p){
          #zjj <- zindx[struc[jj,]]
          zjj <- z[struc[jj,]]
          ndf <- max(round(df * length(zjj) / m), mindf)
          Z <- matrix(0, nrow = m, ncol = ndf)
          Z[struc[jj,],] <- MsplineBasis(x = zjj, df = ndf, degree = degree, intercept = intercept)
          Zmat[[jj]] <- scale(Z, center = FALSE)
          if(periodic) Amat[[jj]] <- t(Zmat[[jj]][m,,drop=FALSE] - Zmat[[jj]][1,,drop=FALSE])
        } # end for(jj in 1:p)
        if(is.null(mode.range)) mode.range <- apply(struc, 1, function(x) range(which(x)))
        if(periodic){
          BB <- mlsun(X = X, Y = Y, Z = Zmat, A = Amat, meq = meq,
                      mode.range = mode.range, maxit = maxit, eps = eps, 
                      del = del, XtX = XtX, simplify = FALSE, catchError = catchError)
        } else {
          BB <- mlsun(X = X, Y = Y, Z = Zmat, mode.range = mode.range,
                      maxit = maxit, eps = eps, del = del, XtX = XtX,
                      simplify = FALSE, catchError = catchError)
        }
        Bdf <- attr(BB, "df")
        B <- matrix(0, nrow = p, ncol = m)
        for(jj in 1:p) B[jj,] <- Zmat[[jj]] %*% BB[[jj]]
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(any(const == c("unismo", "unsmpe")))
    
    # unsmno and unsmpn: unstructured or structured
    if(any(const == c("unsmno","unsmpn"))){
      periodic <- ifelse(const == "unsmpn", TRUE, FALSE)
      #zindx <- seq(0, 1, length.out = m)
      if(is.null(struc)){
        #Z <- MsplineBasis(x = zindx, df = df, degree = degree, intercept = intercept)
        Z <- MsplineBasis(x = z, df = df, degree = degree, intercept = intercept)
        Z <- scale(Z, center = FALSE)
        A <- t(Z[c(1,m),])
        meq <- as.integer(periodic)
        if(periodic) A <- cbind(Z[m,] - Z[1,], A)
        B <- mlsun(X = X, Y = Y, Z = Z, A = A, meq = meq,
                   mode.range = mode.range, maxit = maxit, eps = eps, 
                   del = del, XtX = XtX, catchError = catchError)
        Bdf <- attr(B, "df")
        B <- t(Z %*% B)
      } else {
        Amat <- Zmat <- vector("list", p)
        meq <- rep(as.integer(periodic), p)
        for(jj in 1:p){
          #zjj <- zindx[struc[jj,]]
          zjj <- z[struc[jj,]]
          ndf <- max(round(df * length(zjj) / m), mindf)
          Z <- matrix(0, nrow = m, ncol = ndf)
          Z[struc[jj,],] <- MsplineBasis(x = zjj, df = ndf, degree = degree, intercept = intercept)
          Zmat[[jj]] <- scale(Z, center = FALSE)
          if(periodic){
            Amat[[jj]] <- cbind(Zmat[[jj]][m,] - Zmat[[jj]][1,], t(Zmat[[jj]][c(1,m),]))
          } else {
            Amat[[jj]] <- t(Zmat[[jj]][c(1,m),])
          }
        } # end for(jj in 1:p)
        if(is.null(mode.range)) mode.range <- apply(struc, 1, function(x) range(which(x)))
        BB <- mlsun(X = X, Y = Y, Z = Zmat, A = Amat, meq = meq,
                    mode.range = mode.range, maxit = maxit, eps = eps, 
                    del = del, XtX = XtX, simplify = FALSE, catchError = catchError)
        Bdf <- attr(BB, "df")
        B <- matrix(0, nrow = p, ncol = m)
        for(jj in 1:p) B[jj,] <- Zmat[[jj]] %*% BB[[jj]]
      } # end if(is.null(struc))
      attr(B, "df") <- Bdf
      return(B)
    } # end if(any(const == c("unsmno","unsmpn")))
    
    
  }