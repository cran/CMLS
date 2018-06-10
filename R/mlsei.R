mlsei <- 
  function(X, Y, Z, A, b, meq,
           backfit = FALSE, maxit = 1000, 
           eps = 1e-10, del = 1e-6,
           XtX = NULL, ZtZ = NULL, 
           simplify = TRUE, catchError = FALSE){
    # Solve a Multivariate Least Squares Problem
    # with Equality/Inequality Constraints
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: June 2, 2018
    
    # finds the coefficient matrix B (nu x p) that minimizes
    # sum( ( Y - X %*% t(Z %*% B) )^2 )
    # subject to t(A) %*% B[,j] >= b
    # or 
    # sum( ( Y - sum( X[,j] %*% t(Z[[j]] %*% B[,j]) ) )^2 )
    # subject to t(A[[j]]) %*% B[,j] >= b[[j]]
    
    # Inputs
    #  1) X is n x p predictor matrix (left side).
    #  2) Y is n x m response matrix.
    #  3) Z is m x nu predictor matrix (right side). 
    #     Defaults to m x m identity matrix. 
    #     Can also input a list of length p where  
    #     Z[[j]] is predictor matrix for B[,j].
    #  4) A is nu x r matrix where A[,r] defines the r-th constraint
    #     such that sum(A[,r] * B[,j]) >= b[r] for all j,r.
    #     Can also input a list of length p where  
    #     A[[j]] is contraint matrix for B[,j].
    #  5) b is r x 1 vector defining r-th constraint. 
    #     Defaults to vector of zeros.
    #     Can also input a list of length p where  
    #     b[[j]] is contraint vector for B[,j].
    #  6) meq is a scalar giving # of equality constraints. 
    #     First meq columns are equality (defaults to zero).
    #     Can also input a list of length p where  
    #     meq[[j]] is # of equality constraints for B[,j].
    #  7) backfit is a logical (TRUE = backfit, FALSE = vectorize)
    #     Defaults to scalable computation via backfitting.
    #  8) maxit is maximum # of iterations for backfitting.
    #  9-10) eps and del are epsilon and delta for backfitting.
    #     Convergence:  eps > mean((B - Bold)^2) / (mean(Bold^2) + del)
    # 11-12) XtX = crossprod(X) and ZtZ = crossprod(Z)
    # 13) simply = TRUE will return a matrix B if possible
    #     simply = FALSE will return a list B if Z is a list
    #     Input is ignored if Z is missing or a matrix.
    # 14) errors = TRUE will produce an error (from solve.QP) if constraints are inconsistent.
    #     errors = FALSE uses tryCatch to handle an error from solve.QP
    #     (will return 0 vector if feasible, otherwise NA vector)
    
    # if A is missing, no linear equality constraints are imposed
    
    
    ######***######   CHECK INPUTS   ######***######
    
    ### check X and Y
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(Y)
    if(nrow(X) != n) stop("Inputs 'X' and 'Y' must be matrices with the same number of rows.")
    m <- ncol(Y)
    p <- ncol(X)
    
    ### check Z
    Zmissing <- sameNU <- TRUE
    nu <- rep(m, p)
    if(!missing(Z)){
      Zmissing <- FALSE
      if(is.list(Z)){
        Zlist <- TRUE
        if(length(Z) != p) stop("When transpose = TRUE and 'Z' is a list, it must satisfy:  length(Z) = ncol(X).")
        nu <- rep(0, p)
        for(jj in 1:p){
          Z[[jj]] <- as.matrix(Z[[jj]])
          if(nrow(Z[[jj]]) != m) stop("When transpose = TRUE and 'Z' is a list, each element must satisfy:  nrow(Z[[j]]) = ncol(Y).")
          nu[jj] <- ncol(Z[[jj]])
        }
        sameNU <- ifelse( max(abs(nu[1] - nu)) < .Machine$double.eps , TRUE, FALSE)
      } else {
        Zlist <- FALSE
        Z <- as.matrix(Z)
        if(nrow(Z) != m) stop("Inputs 'Y' and 'Z' must satisfy:  ncol(Y) = nrow(Z)")
        nu <- rep(ncol(Z), p)
        sameNU <- TRUE
      } # end if(is.list(Z))
    } # end if(!missing(Z))
    
    ### check A
    Amissing <- sameNCON <- TRUE
    ncon <- rep(0L, p)
    if(!missing(A)){
      Amissing <- FALSE
      if(is.list(A)){
        Alist <- TRUE
        if(length(A) != p) stop("When 'A' is a list, it must satisfy:  length(A) = ncol(X).")
        ncon <- rep(0, p)
        for(jj in 1:p){
          A[[jj]] <- as.matrix(A[[jj]])
          if(nrow(A[[jj]]) != nu[jj]) stop(ifelse(!Zmissing && Zlist, 
                                                  "When 'A' is a list, each element must satisfy:  nrow(A[[j]]) = ncol(Z[[j]]).",
                                                  "When 'A' is a list, each element must satisfy:  nrow(A[[j]]) = ncol(Z)."))
          ncon[jj] <- ncol(A[[jj]])
        }
        sameNCON <- ifelse( max(abs(ncon[1] - ncon)) < .Machine$double.eps , TRUE, FALSE)
      } else {
        Alist <- FALSE
        if(!sameNU) stop("If 'Z' is a list with ncol(Z[[i]]) != ncol(Z[[j]]) for some i,j\n then 'A' must be a list where the elements satisfy:  nrow(A[[j]]) = ncol(Z[[j]])")
        A <- as.matrix(A)
        if(nrow(A) != nu[1]) stop(ifelse(Zlist,
                                         "When 'A' is a matrix, inputs 'A' and 'Z' must satisfy:  nrow(A) = ncol(Z[[j]])",
                                         "When 'A' is a matrix, inputs 'A' and 'Z' must satisfy:  nrow(A) = ncol(Z)"))
        ncon <- rep(ncol(A), p)
      } # end if(is.list(A))
    } # end if(!missing(A))
    
    ### check b
    if(!Amissing){
      if(missing(b)){
        if(Alist){
          b <- vector("list", p)
          for(jj in 1:p) b[[jj]] <- rep(0, ncon[jj])
        } else {
          b <- rep(0, ncon[1])
        } # end if(Alist)
      } else {
        blist <- ifelse(is.list(b), TRUE, FALSE)
        if(Alist != blist) stop("Inputs 'A' and 'b' must satisfy one of the following:\n 1) t(A) %*% B[,j] >= b, or\n 2) t(A[[j]]) %*% B[,j] >= b[[j]]")
        if(blist){
          if(length(b) != p) stop("When 'b' is a list, it must satisfy:  length(b[[j]]) = ncol(X).")
          for(jj in 1:p){
            b[[jj]] <- as.vector(b[[jj]])
            if(length(b[[jj]]) != ncon[jj]) stop(paste0("Element number ",jj," of the inputs 'A' and 'b' are incompatible:\n  ncol(A[[",jj,"]]) != length(b[[",jj,"]])"))
          }
        } else {
          b <- as.vector(b)
          if(length(b) != ncon[1]) stop("Inputs 'A' and 'b' are incompatible:  ncol(A) != length(b)")
        } # end if(Alist)
      } # end if(missing(b))
    } # end if(!Amissing)
    
    ### check meq
    if(!Amissing){
      if(missing(meq)){
        if(Alist){
          meq <- rep(0, p)
        } else {
          meq <- 0L
        } # end if(Alist)
      } else {
        meq <- unlist(meq)
        if(Alist){
          if(length(meq) != p) stop("When 'A' is a list, it must satisfy:  length(A) = length(meq).")
          for(jj in 1:p){
            meq[jj] <- as.integer(meq[jj])
            if(meq[jj] < 0) stop(paste0("Element number ",jj," of the input 'meq' must be a non-negative integer."))
            if(meq[jj] > ncon[jj]) stop(paste0("Element number ",jj," of the inputs 'A' and 'meq' are incompatible:\n  ncol(A[[",jj,"]]) < meq[",jj,"]"))
          }
        } else {
          meq <- as.integer(meq[1])
          if(meq < 0) stop("Input 'meq' must be a non-negative integer.")
          if(meq > ncon[1]) stop("Inputs 'A' and 'meq' are incompatible:  ncol(A) < meq")
        } # end if(Alist)
      } # end if(missing(b))
    } # end if(!Amissing)
    
    ### check backfitting information
    backfit <- backfit[1]
    if(!any(backfit == c(TRUE, FALSE))) stop("Input 'backfit' must be a logical (TRUE/FALSE).")
    maxit <- as.integer(maxit[1])
    if(maxit < 1L) stop("Input 'maxit' must be a positive integer.")
    eps <- as.numeric(eps[1])
    if(eps < .Machine$double.eps) stop("Input 'eps' must be a positive scalar (greater than or equal to machine epsilon).")
    if(del < 0) stop("Input 'del' must be a positive scalar.")
    
    ### check XtX
    XtXmissing <- TRUE
    if(!is.null(XtX)){
      XtXmissing <- FALSE
      XtX <- as.matrix(XtX)
      if(nrow(XtX) != p | ncol(XtX) != p) stop("Input 'XtX' is must be equal to:  XtX = crossprod(X)")
    }
    
    ### check ZtZ
    ZtZmissing <- TRUE
    if(!Zmissing && !is.null(ZtZ)){
      ZtZmissing <- FALSE
      if(Zlist){
        for(jj in 1:p){
          ZtZ[[jj]] <- as.matrix(ZtZ[[jj]])
          if(nrow(ZtZ[[jj]]) != nu[jj] | ncol(ZtZ[[jj]]) != nu[jj]) stop("Input 'ZtZ[[jj]]' is must be equal to:  ZtZ[[jj]] = crossprod(Z[[jj]])")
        }
      } else {
        ZtZ <- as.matrix(ZtZ)
        if(nrow(ZtZ) != nu[1] | ncol(ZtZ) != nu[1]) stop("Input 'ZtZ' is must be equal to:  ZtZ = crossprod(Z)")
      }
    }
    
    
    ######***######   MISSING A (UNCONSTRAINED)   ######***######
    
    if(Amissing){
      
      if(Zmissing){
        # unconstrained multivariate least squares (Z = identity)
        
        B <- crossprod(Y, PseudoInverse(t(X), symmetric = FALSE))
        Bdf <- rep(m, p)
        attr(B, "df") <- Bdf
        return(B)
        
      } else {
        # both sides least squares (Z != identity)
        
        # is Z a list?
        if(Zlist){
          
          # setup LS problem
          if(XtXmissing) XtX <- crossprod(X)
          dXtX <- diag(XtX)
          B <- Bold <- Rinv <- vector("list", p)
          for(jj in 1:p) {
            B[[jj]] <- Bold[[jj]] <- rep(0, nu[jj])
            if(ZtZmissing){
              Rinv[[jj]] <- tcrossprod(rsinvsqrt(crossprod(Z[[jj]])))
            } else {
              Rinv[[jj]] <- tcrossprod(rsinvsqrt(ZtZ[[jj]]))
            }
          }
          varid <- 1:p
          
          # get crossproducts
          YtX <- crossprod(Y, X)
          ZtYtX <- vector("list", p)
          for(jj in 1:p) ZtYtX[[jj]] <- crossprod(Z[[jj]], YtX)
          
          # solve LS problem (via backfitting)
          epschk <- eps + 1
          iter <- 0
          while(epschk > eps & iter < maxit){
            for(jj in 1:p){
              Hmat <- matrix(0, nrow = m, ncol = 1)
              if(p > 1) for(hh in varid[-jj]) Hmat <- Hmat + Z[[hh]] %*% B[[hh]] * XtX[hh,jj]
              dvec <- (ZtYtX[[jj]][,jj] - crossprod(Z[[jj]], Hmat)) / dXtX[jj]
              B[[jj]] <- Rinv[[jj]] %*% dvec
            } # end for(jj in 1:p)
            epschk <- mean((unlist(B) - unlist(Bold))^2) / (mean(unlist(Bold)^2) + del)
            Bold <- B
            iter <- iter + 1L
          } # end while(epschk > eps & iter < maxit)
          
          # return solution
          Bdf <- nu
          if(sameNU && simplify) B <- matrix(unlist(B), nrow = nu[1], ncol = p)
          attr(B, "df") <- Bdf
          attr(B, "iter") <- iter
          return(B)
          
        } else {
          
          B <- PseudoInverse(Z, symmetric = FALSE) %*% crossprod(Y, PseudoInverse(t(X), symmetric = FALSE))
          Bdf <- rep(nu[1], p)
          attr(B, "df") <- Bdf
          return(B)
          
        } # end if(Zlist)
        
      } # end if(Zmissing)
      
    } # end if(Amissing)
    
    
    ######***######   NON-MISSING A (E/I CONSTRAINTS)   ######***######
    
    if(Zmissing){
      # Z = identity
      
      if(backfit){
        # scalable computation via back-fitting
        
        # setup QP problem
        if(XtXmissing) XtX <- crossprod(X)
        dXtX <- diag(XtX)
        B <- Bold <- matrix(0, nrow = m, ncol = p)
        Bdf <- rep(0, p)
        Rinv <- diag(m)
        if(!Alist){
          Amat <- A
          bvec <- b
          meqc <- meq
        }
        
        # get crossproducts
        YtX <- crossprod(Y, X)
        
        # solve QP problem (via backfitting)
        epschk <- eps + 1
        iter <- 0
        while(epschk > eps & iter < maxit){
          for(jj in 1:p){
            if(p > 1){
              dvec <- (YtX[,jj] - B[,-jj,drop=FALSE] %*% XtX[-jj,jj,drop=FALSE]  ) / dXtX[jj]
            } else {
              dvec <- YtX[,jj] / dXtX[jj]
            }
            if(Alist){
              Amat <- A[[jj]]
              bvec <- b[[jj]]
              meqc <- meq[jj]
            }
            Bfit <- solveQP(Dmat = Rinv, dvec = dvec, Amat = Amat, bvec = bvec, meq = meqc, 
                            factorized = TRUE, catchError = catchError)
            B[,jj] <- Bfit$solution
            nact <- length(Bfit$iact)
            if(nact == 1L && Bfit$iact == 0) nact <- 0L
            Bdf[jj] <- max(m - nact, 0)
          } # end for(jj in 1:p)
          epschk <- mean((B - Bold)^2) / (mean(Bold^2) + del)
          Bold <- B
          iter <- iter + 1L
        } # end while(epschk > eps & iter < maxit)
        
        # return solution
        attr(B, "df") <- Bdf
        attr(B, "iter") <- iter
        return(B)
        
      } else {
        # simultaneous estimation via vectorization
        
        # setup QP problem
        betas <- paste0("beta", 1:p)
        conlab <- factor(rep(betas, ncon), levels = betas)
        if(Alist){
          Amat <- list2bdiag(A)
          bvec <- unlist(b)
          meqc <- sum(meq)
          nconcs <- c(0, cumsum(ncon))
          if(any(meq > 0)){
            meqid <- NULL
            for(jj in 1:p) if(meq[jj] > 0) meqid <- c(meqid, 1:meq[jj] + nconcs[jj])
            idvec <- 1:ncol(Amat)
            meqidx <- c(meqid, idvec[-meqid])
            Amat <- Amat[,meqidx,drop=FALSE]
            bvec <- bvec[meqidx]
            conlab <- conlab[meqidx]
          }
        } else {
          Amat <- kronecker(diag(p), A)
          bvec <- rep(b, p)
          meqc <- meq * p
          if(meq > 0){
            meqid <- NULL
            for(jj in 1:p) meqid <- c(meqid, 1:meq + (jj - 1) * ncon)
            idvec <- 1:ncol(Amat)
            meqidx <- c(meqid, idvec[-meqid])
            Amat <- Amat[,meqidx,drop=FALSE]
            bvec <- bvec[meqidx]
            conlab <- conlab[meqidx]
          }
        }
        if(XtXmissing) XtX <- crossprod(X)
        Rinv <- kronecker(rsinvsqrt(XtX), diag(m))
        dvec <- c(crossprod(Y, X))
        
        # solve QP problem (via vectorization)
        Bfit <- solveQP(Dmat = Rinv, dvec = dvec, Amat = Amat, bvec = bvec, meq = meqc, 
                        factorized = TRUE, catchError = catchError)
        B <- matrix(Bfit$solution, nrow = m, ncol = p)
        
        # return solution
        Bdf <- rep(m, p)
        nact <- length(Bfit$iact)
        if(nact == 1L && Bfit$iact == 0) nact <- 0L
        if(nact > 0L){
          for(jj in 1:p){
            nactj <- sum(conlab[Bfit$iact] == betas[jj])
            Bdf[jj] <- max(Bdf[jj] - nactj, 0)
          }
        }
        attr(B, "df") <- Bdf
        return(B)
        
      } # end if(backfit)
      
    } else {
      
      # is Z a list?
      if(Zlist){
        # back-fitting is the only option
        
        # setup QP problem
        if(XtXmissing) XtX <- crossprod(X)
        dXtX <- diag(XtX)
        B <- Bold <- Rinv <- vector("list", p)
        Bdf <- rep(0, p)
        for(jj in 1:p) {
          B[[jj]] <- Bold[[jj]] <- rep(0, nu[jj])
          if(ZtZmissing){
            Rinv[[jj]] <- rsinvsqrt(crossprod(Z[[jj]]))
          } else {
            Rinv[[jj]] <- rsinvsqrt(ZtZ[[jj]])
          }
        }
        varid <- 1:p
        if(!Alist){
          Amat <- A
          bvec <- b
          meqc <- meq
        }
        
        # get crossproducts
        YtX <- crossprod(Y, X)
        ZtYtX <- vector("list", p)
        for(jj in 1:p) ZtYtX[[jj]] <- crossprod(Z[[jj]], YtX)
        
        # solve QP problem (via backfitting)
        epschk <- eps + 1
        iter <- 0
        while(epschk > eps & iter < maxit){
          for(jj in 1:p){
            Hmat <- matrix(0, nrow = m, ncol = 1)
            if(p > 1) for(hh in varid[-jj]) Hmat <- Hmat + Z[[hh]] %*% B[[hh]] * XtX[hh,jj]
            dvec <- (ZtYtX[[jj]][,jj] - crossprod(Z[[jj]], Hmat)) / dXtX[jj]
            if(Alist){
              Amat <- A[[jj]]
              bvec <- b[[jj]]
              meqc <- meq[jj]
            }
            Bfit <- solveQP(Dmat = Rinv[[jj]], dvec = dvec, Amat = Amat, bvec = bvec, meq = meqc, 
                            factorized = TRUE, catchError = catchError)
            B[[jj]] <- Bfit$solution
            nact <- length(Bfit$iact)
            if(nact == 1L && Bfit$iact == 0) nact <- 0L
            Bdf[jj] <- max(nu[jj] - nact, 0)
          } # end for(jj in 1:p)
          epschk <- mean((unlist(B) - unlist(Bold))^2) / (mean(unlist(Bold)^2) + del)
          Bold <- B
          iter <- iter + 1L
        } # end while(epschk > eps & iter < maxit)
        
        # return solution
        if(sameNU && simplify) B <- matrix(unlist(B), nrow = nu[1], ncol = p)
        attr(B, "df") <- Bdf
        attr(B, "iter") <- iter
        return(B)
        
      } else {
        
        # Z is a matrix: backfit or vectorize?
        if(backfit){
          # scalable computation via back-fitting
          
          # setup QP problem
          if(XtXmissing) XtX <- crossprod(X)
          dXtX <- diag(XtX)
          B <- Bold <- matrix(0, nrow = nu[1], ncol = p)
          Bdf <- rep(0, p)
          if(ZtZmissing) ZtZ <- crossprod(Z)
          Rinv <- rsinvsqrt(ZtZ)
          if(!Alist){
            Amat <- A
            bvec <- b
            meqc <- meq
          }
          
          # get crossproducts
          ZtYtX <- crossprod(Z, crossprod(Y, X))
          
          # solve QP problem (via backfitting)
          epschk <- eps + 1
          iter <- 0
          while(epschk > eps & iter < maxit){
            for(jj in 1:p){
              if(p > 1){
                dvec <- (ZtYtX[,jj] - ZtZ %*% B[,-jj,drop=FALSE] %*% XtX[-jj,jj,drop=FALSE]  ) / dXtX[jj]
              } else {
                dvec <- ZtYtX[,jj] / dXtX[jj]
              }
              if(Alist){
                Amat <- A[[jj]]
                bvec <- b[[jj]]
                meqc <- meq[jj]
              }
              Bfit <- solveQP(Dmat = Rinv, dvec = dvec, Amat = Amat, bvec = bvec, meq = meqc, 
                              factorized = TRUE, catchError = catchError)
              B[,jj] <- Bfit$solution
              nact <- length(Bfit$iact)
              if(nact == 1L && Bfit$iact == 0) nact <- 0L
              Bdf[jj] <- max(nu[1] - nact, 0)
            } # end for(jj in 1:p)
            epschk <- mean((B - Bold)^2) / (mean(Bold^2) + del)
            Bold <- B
            iter <- iter + 1L
          } # end while(epschk > eps & iter < maxit)
          
          # return solution
          attr(B, "df") <- Bdf
          attr(B, "iter") <- iter
          return(B)
          
        } else {
          # simultaneous estimation via vectorization
          
          # setup QP problem
          betas <- paste0("beta", 1:p)
          conlab <- factor(rep(betas, ncon), levels = betas)
          if(Alist){
            Amat <- list2bdiag(A)
            bvec <- unlist(b)
            meqc <- sum(meq)
            nconcs <- c(0, cumsum(ncon))
            if(any(meq > 0)){
              meqid <- NULL
              for(jj in 1:p) if(meq[jj] > 0) meqid <- c(meqid, 1:meq[jj] + nconcs[jj])
              idvec <- 1:ncol(Amat)
              meqidx <- c(meqid, idvec[-meqid])
              Amat <- Amat[,meqidx,drop=FALSE]
              bvec <- bvec[meqidx]
              conlab <- conlab[meqidx]
            }
          } else {
            Amat <- kronecker(diag(p), A)
            bvec <- rep(b, p)
            meqc <- meq * p
            if(meq > 0){
              meqid <- NULL
              for(jj in 1:p) meqid <- c(meqid, 1:meq + (jj - 1) * ncon)
              idvec <- 1:ncol(Amat)
              meqidx <- c(meqid, idvec[-meqid])
              Amat <- Amat[,meqidx,drop=FALSE]
              bvec <- bvec[meqidx]
              conlab <- conlab[meqidx]
            }
          }
          if(XtXmissing) XtX <- crossprod(X)
          if(ZtZmissing) ZtZ <- crossprod(Z)
          Rinv <- kronecker(rsinvsqrt(XtX), rsinvsqrt(ZtZ))
          dvec <- c(crossprod(Z, crossprod(Y, X)))
          
          # solve QP problem (via vectorization)
          Bfit <- solveQP(Dmat = Rinv, dvec = dvec, Amat = Amat, bvec = bvec, meq = meqc, 
                          factorized = TRUE, catchError = catchError)
          B <- matrix(Bfit$solution, nrow = nu[1], ncol = p)
          
          # return solution
          Bdf <- nu
          nact <- length(Bfit$iact)
          if(nact == 1L && Bfit$iact == 0) nact <- 0L
          if(nact > 0L){
            for(jj in 1:p){
              nactj <- sum(conlab[Bfit$iact] == betas[jj])
              Bdf[jj] <- max(Bdf[jj] - nactj, 0)
            }
          }
          attr(B, "df") <- Bdf
          return(B)
          
        } # end if(backfit)
        
      } # end if(Zlist)
      
    } # end if(Zmissing)
    
  }