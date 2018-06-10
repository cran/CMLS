solveQP <- 
  function(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE, catchError = TRUE){
    
    if(missing(bvec)) bvec <- rep(0, ncol(Amat))
    if(catchError){
      return(tryCatch(solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq, factorized = factorized),
                      error = function(e) solveQPerror(Amat, bvec, meq)))
    } else {
      return(solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq, factorized = factorized))
    } # end if(catchError)
    
  } # end solveQP