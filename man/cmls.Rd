\name{cmls}
\alias{cmls}
\title{
Solve a Constrained Multivariate Least Squares Problem
}
\description{
Finds the \eqn{p} x \eqn{m} matrix \code{B} that minimizes the multivariate least squares problem
\tabular{c}{
\code{ sum(( Y - X \%*\% B )^2) }
} 
subject to the specified constraints on the rows of \code{B}.
}
\usage{
cmls(X, Y, const = "uncons", struc = NULL, 
     df = 10, degree = 3, intercept = TRUE,
     backfit = FALSE, maxit = 1e3, eps = 1e-10, 
     del = 1e-6, XtX = NULL, mode.range = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Matrix of dimension \eqn{n} x \eqn{p}.}
  \item{Y}{Matrix of dimension \eqn{n} x \eqn{m}.}
  \item{const}{Constraint code. See \code{\link{const}} for the 24 available options.}
  \item{struc}{Structural constraints (defaults to unstructured). See Note.}
  \item{df}{Degrees of freedom for the spline basis (for smoothness constraints). See Note.}
  \item{degree}{Polynomial degree for the spline basis (for smoothness constraints). See Note.}
  \item{intercept}{Logical indicating whether the spline basis should contain an intercept (for smoothness constraints). See Note.}
  \item{backfit}{Estimate \code{B} via back-fitting (\code{TRUE}) or vectorization (\code{FALSE}). See Details.}
  \item{maxit}{Maximum number of iterations for back-fitting algorithm. Ignored if \code{backfit = FALSE}.}
  \item{eps}{Convergence tolerance for back-fitting algorithm. Ignored if \code{backfit = FALSE}.}
  \item{del}{Stability tolerance for back-fitting algorithm. Ignored if \code{backfit = FALSE}.}
  \item{XtX}{Crossproduct matrix:  \code{XtX = crossprod(X)}.}
  \item{mode.range}{Mode search ranges (for unimodal constraints). See Note.}
}
\details{
If \code{backfit = FALSE} (default), a closed-form solution is used to estimate \code{B} whenever possible. Otherwise a back-fitting algorithm is used, where the rows of \code{B} are updated sequentially until convergence. The backfitting algorithm is determined to have converged when 

\code{mean((B.new - B.old)^2) < eps * (mean(B.old^2) + del)}, 

where \code{B.old} and \code{B.new} denote the parameter estimates at iterations \eqn{t} and \eqn{t+1} of the backfitting algorithm.
}
\value{
Returns the estimated matrix \code{B} with attribute "df" (degrees of freedom), which gives the df for each row of \code{B}.
}
\references{
Goldfarb, D., & Idnani, A. (1983). A numerically stable dual method for solving strictly convex quadratic programs. Mathematical Programming, 27, 1-33.

Helwig, N. E. (in prep). Constrained multivariate least squares in R.

Ten Berge, J. M. F. (1993). Least Squares Optimization in Multivariate Analysis. Volume 25 of M & T Series. DSWO Press, Leiden University. ISBN: 9789066950832

Turlach, B. A., & Weingessel, A. (2013). quadprog: Functions to solve Quadratic Programming Problems. R package version 1.5-5. https://CRAN.R-project.org/package=quadprog
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Structure constraints (\code{struc}) should be specified with a \eqn{p} x \eqn{m} matrix of logicals (TRUE/FALSE), such that FALSE elements indicate a weight should be constrained to be zero. Default uses unstructured weights, i.e., a \eqn{p} x \eqn{m} matrix of all TRUE values.

Inputs \code{df}, \code{degree}, and \code{intercept} are only applicable when using one of the 12 constraints that involves a spline basis, i.e., "smooth", "smonon", "smoper", "smpeno", "ortsmo", "orsmpe", "monsmo", "mosmno", "unismo", "unsmno", "unsmpe", "unsmpn".

Input \code{mode.range} is only applicable when using one of the 8 constraints that enforces unimodality: "unimod", "uninon", "uniper", "unpeno", "unismo", "unsmno", "unsmpe", "unsmpn". Mode search ranges (\code{mode.range}) should be specified with a 2 x \eqn{p} matrix of integers such that 

\code{1 <= mode.range[1,j] <= mode.range[2,j] <= m} for all \code{j = 1:p}. 

Default is  \code{mode.range = matrix(c(1, m), 2, p)}.
}

\seealso{
\code{\link{const}} prints/returns the contraint options.

\code{\link{cv.cmls}} performs k-fold cross-validation to tune the constraint options.

\code{\link{mlsei}} and \code{\link{mlsun}} are used to implement several of the constraints.
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
Bhat <- cmls(X = Xmat, Y = Ymat, const = "uncons")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unconstrained and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "uncons", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")


######***######   NON-NEGATIVITY   ######***######

# non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "nonneg")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "nonneg", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")


######***######   PERIODICITY   ######***######

# periodic
Bhat <- cmls(X = Xmat, Y = Ymat, const = "period")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# periodic and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "period", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# periodic and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "pernon")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# periodic and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "pernon", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")


######***######   SMOOTHNESS   ######***######

# smooth
Bhat <- cmls(X = Xmat, Y = Ymat, const = "smooth")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# smooth and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "smooth", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# smooth and periodic
Bhat <- cmls(X = Xmat, Y = Ymat, const = "smoper")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# smooth and periodic and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "smoper", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# smooth and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "smonon")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# smooth and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "smonon", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# smooth and periodic and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "smpeno")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# smooth and periodic and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "smpeno", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")


######***######   ORTHOGONALITY   ######***######

# orthogonal
Bhat <- cmls(X = Xmat, Y = Ymat, const = "orthog")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# orthogonal and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "orthog", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# orthgonal and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "ortnon")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# orthgonal and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "ortnon", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# orthogonal and smooth
Bhat <- cmls(X = Xmat, Y = Ymat, const = "ortsmo")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# orthogonal and smooth and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "ortsmo", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# orthogonal and smooth and periodic
Bhat <- cmls(X = Xmat, Y = Ymat, const = "orsmpe")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# orthogonal and smooth and periodic and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "orsmpe", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")


######***######   UNIMODALITY   ######***######

# unimodal
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unimod")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unimod", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "uninon")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "uninon", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and periodic
Bhat <- cmls(X = Xmat, Y = Ymat, const = "uniper")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and periodic and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "uniper", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and periodic and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unpeno")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and periodic and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unpeno", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")


######***######   UNIMODALITY AND SMOOTHNESS   ######***######

# unimodal and smooth
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unismo")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and smooth and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unismo", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and smooth and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unsmno")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and smooth and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unsmno", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and smooth and periodic
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unsmpe")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and smooth and periodic and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unsmpe", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and smooth and periodic and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unsmpn")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# unimodal and smooth and periodic and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "unsmpn", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")


######***######   MONOTONICITY   ######***######

# make B
x <- 1:m
Bmat <- rbind(1 / (1 + exp(-(x - quantile(x, 0.5)))), 
              1 / (1 + exp(-(x - quantile(x, 0.8)))))
struc <- rbind(rep(c(FALSE, TRUE), c(1 * m, 3 * m) / 4),
               rep(c(FALSE, TRUE), c(m, m) / 2))
Bmat <- Bmat * struc               

# make noisy data
set.seed(1)
Ymat <- Xmat \%*\% Bmat + rnorm(m*n, sd = 0.25)

# monotonic increasing
Bhat <- cmls(X = Xmat, Y = Ymat, const = "moninc")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# monotonic increasing and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "moninc", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# monotonic increasing and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "monnon")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# monotonic increasing and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "monnon", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# monotonic increasing and smooth
Bhat <- cmls(X = Xmat, Y = Ymat, const = "monsmo")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# monotonic increasing and smooth and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "monsmo", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# monotonic increasing and smooth and non-negative
Bhat <- cmls(X = Xmat, Y = Ymat, const = "mosmno")
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

# monotonic increasing and smooth and non-negative and structured
Bhat <- cmls(X = Xmat, Y = Ymat, const = "mosmno", struc = struc)
mean((Bhat - Bmat)^2)
attr(Bhat, "df")

}

\keyword{ models }
\keyword{ multivariate }
\keyword{ regression }
\keyword{ optimize }
\keyword{ smooth }
