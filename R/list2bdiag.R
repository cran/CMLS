list2bdiag <- 
  function(x){
    
    if(!is.list(x)) stop("Input 'x' must be a list.")
    xrow <- xcol <- rep(0, length(x))
    square <- TRUE
    for(j in 1:length(x)){
      x[[j]] <- as.matrix(x[[j]])
      xrow[j] <- nrow(x[[j]])
      xcol[j] <- ncol(x[[j]])
      if(xrow[j] != xcol[j]) square <- FALSE
    }
    
    if(square){
      bdim <- sum(xrow)
      csdim <- c(0, cumsum(xrow))
      bdmat <- matrix(0, nrow = bdim, ncol = bdim)
      for(j in 1:length(x)){
        indx <- seq(csdim[j] + 1, csdim[j+1])
        bdmat[indx, indx] <- x[[j]]
      }
    } else {
      bdrow <- sum(xrow)
      csrow <- c(0, cumsum(xrow))
      bdcol <- sum(xcol)
      cscol <- c(0, cumsum(xcol))
      bdmat <- matrix(0, nrow = bdrow, ncol = bdcol)
      for(j in 1:length(x)){
        rowindx <- seq(csrow[j] + 1, csrow[j+1])
        colindx <- seq(cscol[j] + 1, cscol[j+1])
        bdmat[rowindx, colindx] <- x[[j]]
      }
      
    } # end if(square)
    
    bdmat
    
  }