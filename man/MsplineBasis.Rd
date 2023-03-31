\name{MsplineBasis}
\alias{MsplineBasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M-Spline Basis for Polynomial Splines
}
\description{
Generate the M-spline basis matrix for a polynomial spline.
}
\usage{
MsplineBasis(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
             Boundary.knots = range(x), periodic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the predictor variable. Missing values are \bold{not} allowed.
}
  \item{df}{
degrees of freedom; if specified the number of \code{knots} is defined as \code{df - degree - ifelse(intercept, 1, 0)}; the \code{knots} are placed at the quantiles of \code{x}
}
  \item{knots}{
the internal breakpoints that define the spline (typically the quantiles of \code{x}) 
}
  \item{degree}{
degree of the piecewise polynomial---default is 3 for cubic splines
}
  \item{intercept}{
if \code{TRUE}, the basis includes an intercept column  
}
  \item{Boundary.knots}{
boundary points for M-spline basis; defaults to min and max of \code{x}  
}
  \item{periodic}{
if \code{TRUE}, the M-spline basis is constrained to be periodic  
}
}
\details{
Syntax is adapted from the \code{bs} function in the \bold{splines} package (R Core Team, 2021).

Used for implementing various types of smoothness constraints in the \code{\link{cmls}} fucntion.
}
\value{
A matrix of dimension \code{c(length(x), df)} where either \code{df} was supplied or \code{df = length(knots) + degree + ifelse(intercept, 1, 0)}
}
\references{
R Core Team (2023). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Ramsay, J. O. (1988). Monotone regression splines in action. \emph{Statistical Science, 3}, 425-441. \doi{10.1214/ss/1177012761}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{IsplineBasis}}
}
\examples{
x <- seq(0, 1, length.out = 101)
M <- MsplineBasis(x, df = 8, intercept = TRUE)
M <- scale(M, center = FALSE)
plot(x, M[,1], ylim = range(M), t = "l")
for(j in 2:8) lines(x, M[,j], col = j)
}

\keyword{ regression }
\keyword{ smooth }
