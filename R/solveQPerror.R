solveQPerror <- 
  function(Amat, bvec, meq){
    ncoef <- nrow(Amat)
    ncons <- ncol(Amat)
    Atz <- rep(0, ncons)
    if(meq > 0){
      if(meq == ncons){
        if(any(abs(Atz - bvec) > ncoef * .Machine$double.eps)){
          sol <- list(solution = rep(NA, ncoef), value = NA, iact = NA)
        } else {
          sol <- list(solution = rep(0, ncoef), value = 0, iact = 1:ncons)
        }
      } else {
        eqcon <- any(abs(Atz[1:meq] - bvec[1:meq]) > ncoef * .Machine$double.eps)
        incon <- any(Atz[(meq+1):ncons] < bvec[(meq+1):ncons])
        if(eqcon | incon){
          sol <- list(solution = rep(NA, ncoef), value = NA, iact = NA)
        } else {
          sol <- list(solution = rep(0, ncoef), value = 0, iact = 1:ncons)
        }
      }
    } else {
      if(any(Atz < bvec)){
        sol <- list(solution = rep(NA, ncoef), value = NA, iact = NA)
      } else {
        sol <- list(solution = rep(0, ncoef), value = 0, iact = 1:ncons)
      }
    }
    sol
  } # end solveQPerror