\name{mlsei}
\alias{mlsei}
\title{
Multivariate Least Squares with Equality/Inequality Constraints
}
\description{
Finds the \eqn{q} x \eqn{p} matrix \code{B} that minimizes the multivariate least squares problem
\tabular{c}{
\code{ sum(( Y - X \%*\% t(Z \%*\% B) )^2) }
} 
subject to \code{t(A) \%*\% B[,j] >= b} for all \code{j = 1:p}. Unique basis functions and constraints are allowed for each column of \code{B}.
}
\usage{
mlsei(X, Y, Z, A, b, meq,
      backfit = FALSE, maxit = 1000, 
      eps = 1e-10, del = 1e-6,
      XtX = NULL, ZtZ = NULL, 
      simplify = TRUE, catchError = FALSE)
}
\arguments{
  \item{X}{Matrix of dimension \eqn{n} x \eqn{p}.}
  \item{Y}{Matrix of dimension \eqn{n} x \eqn{m}.}
  \item{Z}{Matrix of dimension \eqn{m} x \eqn{q}. Can also input a list (see Note). If missing, then \code{Z = diag(m)} so that \eqn{q = m}.}
  \item{A}{Constraint matrix of dimension \eqn{q} x \eqn{r}. Can also input a list (see Note). If missing, no constraints are imposed.}
  \item{b}{Consraint vector of dimension \eqn{r} x 1. Can also input a list (see Note). If missing, then \code{b = rep(0, r)}.}
  \item{meq}{The first \code{meq} columns of \code{A} are equality constraints, and the remaining \code{r - meq} are inequality constraints. Can also input a vector (see Note). If missing, then \code{meq = 0}.}
  \item{backfit}{Estimate \code{B} via back-fitting (\code{TRUE}) or vectorization (\code{FALSE}). See Details.}
  \item{maxit}{Maximum number of iterations for back-fitting algorithm. Ignored if \code{backfit = FALSE}.}
  \item{eps}{Convergence tolerance for back-fitting algorithm. Ignored if \code{backfit = FALSE}.}
  \item{del}{Stability tolerance for back-fitting algorithm. Ignored if \code{backfit = FALSE}.}
  \item{XtX}{Crossproduct matrix:  \code{XtX = crossprod(X)}.}
  \item{ZtZ}{Crossproduct matrix:  \code{ZtZ = crossprod(Z)}.}
  \item{simplify}{If \code{Z} is a list, should \code{B} be returned as a matrix (if possible)? See Note.}
  \item{catchError}{If \code{catchError = FASLE}, an error induced by \code{\link{solve.QP}} will be returned. Otherwise \code{\link{tryCatch}} will be used in attempt to catch the error.}

}
\details{
If \code{backfit = FALSE} (default), a closed-form solution is used to estimate \code{B} whenever possible. Otherwise a back-fitting algorithm is used, where the columns of \code{B} are updated sequentially until convergence. The backfitting algorithm is determined to have converged when 

\code{mean((B.new - B.old)^2) < eps * (mean(B.old^2) + del)}, 

where \code{B.old} and \code{B.new} denote the parameter estimates at iterations \eqn{t} and \eqn{t+1} of the backfitting algorithm.
}
\value{
If \code{Z} is a list with \eqn{q_j = q} for all \eqn{j = 1,\ldots,p}, then...
\item{\code{B}}{is returned as a \eqn{q} x \eqn{p} matrix when \code{simplify = TRUE}}
\item{\code{B}}{is returned as a list of length \eqn{p} when \code{simplify = FALSE}}

If \code{Z} is a list with \eqn{q_j \neq q} for some \eqn{j}, then \code{B} is returned as a list of length \eqn{p}.

Otherwise \code{B} is returned as a \eqn{q} x \eqn{p} matrix.

}
\references{
Goldfarb, D., & Idnani, A. (1983). A numerically stable dual method for solving strictly convex quadratic programs. Mathematical Programming, 27, 1-33. \doi{10.1007/BF02591962}

Helwig, N. E. (in prep). Constrained multivariate least squares in R.

Ten Berge, J. M. F. (1993). Least Squares Optimization in Multivariate Analysis. Volume 25 of M & T Series. DSWO Press, Leiden University. ISBN: 9789066950832

Turlach, B. A., & Weingessel, A. (2019). quadprog: Functions to solve Quadratic Programming Problems. R package version 1.5-8. https://CRAN.R-project.org/package=quadprog
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The \code{Z} input can also be a list of length \eqn{p} where \code{Z[[j]]} contains a \eqn{m} x \eqn{q_j} matrix. If \eqn{q_j = q} for all \eqn{j = 1,\ldots,p} and \code{simplify = TRUE}, the output \code{B} will be a matrix. Otherwise \code{B} will be a list of length \eqn{p} where \code{B[[j]]} contains a \eqn{q_j} x 1 vector.

The \code{A} and \code{b} inputs can also be lists of length \eqn{p} where \code{t(A[[j]]) \%*\% B[,j] >= b[[j]]} for all \eqn{j = 1,\ldots,p}. If \code{A} and \code{b} are lists of length \eqn{p}, the \code{meq} input should be a vector of length \eqn{p} indicating the number of equality constraints for each element of \code{A}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cmls}} calls this function for several of the constraints.
}
\examples{
######***######   GENERATE DATA   ######***######

# make X
set.seed(2)
n <- 50
m <- 20
p <- 2
Xmat <- matrix(rnorm(n*p), nrow = n, ncol = p)

# make B (which satisfies all constraints except monotonicity)
x <- seq(0, 1, length.out = m)
Bmat <- rbind(sin(2*pi*x), sin(2*pi*x+pi)) / sqrt(4.75)
struc <- rbind(rep(c(TRUE, FALSE), each = m / 2),
               rep(c(FALSE, TRUE), each = m / 2))
Bmat <- Bmat * struc

# make noisy data
set.seed(1)
Ymat <- Xmat \%*\% Bmat + rnorm(n*m, sd = 0.25)


######***######   UNCONSTRAINED   ######***######

# unconstrained
Bhat.cmls <- cmls(X = Xmat, Y = Ymat, const = "uncons")
Bhat.mlsei <- t(mlsei(X = Xmat, Y = Ymat))
mean((Bhat.cmls - Bhat.mlsei)^2)

# unconstrained and structured (note: cmls is more efficient)
Bhat.cmls <- cmls(X = Xmat, Y = Ymat, const = "uncons", struc = struc)
Amat <- vector("list", p)
meq <- rep(0, p)
for(j in 1:p){
   meq[j] <- sum(!struc[j,])
   if(meq[j] > 0){
      A <- matrix(0, nrow = m, ncol = meq[j])
      A[!struc[j,],] <- diag(meq[j])
      Amat[[j]] <- A
   } else {
      Amat[[j]] <- matrix(0, nrow = m, ncol = 1)
   }
}
Bhat.mlsei <- t(mlsei(X = Xmat, Y = Ymat, A = Amat, meq = meq))
mean((Bhat.cmls - Bhat.mlsei)^2)


######***######   NON-NEGATIVITY   ######***######

# non-negative
Bhat.cmls <- cmls(X = Xmat, Y = Ymat, const = "nonneg")
Bhat.mlsei <- t(mlsei(X = Xmat, Y = Ymat, A = diag(m)))
mean((Bhat.cmls - Bhat.mlsei)^2)

# non-negative and structured (note: cmls is more efficient)
Bhat.cmls <- cmls(X = Xmat, Y = Ymat, const = "nonneg", struc = struc)
eye <- diag(m)
meq <- rep(0, p)
for(j in 1:p){
   meq[j] <- sum(!struc[j,])
   Amat[[j]] <- eye[,sort(struc[j,], index.return = TRUE)$ix]
}
Bhat.mlsei <- t(mlsei(X = Xmat, Y = Ymat, A = Amat, meq = meq))
mean((Bhat.cmls - Bhat.mlsei)^2)


# see internals of cmls.R for further examples

}

\keyword{ models }
\keyword{ multivariate }
\keyword{ regression }
\keyword{ optimize }
\keyword{ smooth }
