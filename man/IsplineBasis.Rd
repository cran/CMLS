\name{IsplineBasis}
\alias{IsplineBasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
I-Spline Basis for Monotonic Polynomial Splines
}
\description{
Generate the I-spline basis matrix for a monotonic polynomial spline.
}
\usage{
IsplineBasis(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
             Boundary.knots = range(x))
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
degree of the M-spline basis---default is 3 for cubic splines
}
  \item{intercept}{
if \code{TRUE}, the basis includes an intercept column  
}
  \item{Boundary.knots}{
boundary points for M-spline basis; defaults to min and max of \code{x}  
}
}
\details{
Syntax is adapted from the \code{bs} function in the \bold{splines} package (R Core Team, 2021).

Used for implementing monotonic smoothness constraints in the \code{\link{cmls}} fucntion.
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
\note{
I-spline basis functions are created by integrating M-spline basis functions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MsplineBasis}}
}
\examples{
x <- seq(0, 1, length.out = 101)
I <- IsplineBasis(x, df = 8, intercept = TRUE)
plot(x, I[,1], ylim = c(0, 1), t = "l")
for(j in 2:8) lines(x, I[,j], col = j)
}

\keyword{ regression }
\keyword{ smooth }
