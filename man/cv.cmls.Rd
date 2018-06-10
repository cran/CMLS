\name{cv.cmls}
\alias{cv.cmls}
\title{
Cross-Validation for cmls
}
\description{
Does k-fold or generalized cross-validation to tune the constraint options for \code{\link{cmls}}. Tunes the model with respect to any combination of the arguments \code{const}, \code{df}, \code{degree}, and/or \code{intercept}.
}
\usage{
cv.cmls(X, Y, nfolds = 2, foldid = NULL, parameters = NULL,
        const = "uncons", df = 10, degree = 3, intercept = TRUE,
        mse = TRUE, parallel = FALSE, cl = NULL, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Matrix of dimension \eqn{n} x \eqn{p}.}
  \item{Y}{Matrix of dimension \eqn{n} x \eqn{m}.}
  \item{nfolds}{Number of folds for k-fold cross-validation. Ignored if \code{foldid} argument is provided. Set \code{nfolds=1} for generalized cross-validation (GCV).}
  \item{foldid}{Factor or integer vector of length \eqn{n} giving the fold identification for each observation.}
  \item{parameters}{Parameters for tuning. Data frame with columns \code{const}, \code{df}, \code{degree}, and \code{intercept}. See Details.}
  \item{const}{Parameters for tuning. Character vector specifying constraints for tuning. See Details.}
  \item{df}{Parameters for tuning. Integer vector specifying degrees of freedom for tuning. See Details.}
  \item{degree}{Parameters for tuning. Integer vector specifying polynomial degrees for tuning. See Details.}
  \item{intercept}{Parameters for tuning. Logical vector specifying intercepts for tuning. See Details.}
  \item{mse}{If \code{TRUE} (default), the mean squared error is used as the CV loss function. Otherwise the mean absolute error is used.}
  \item{parallel}{Logical indicating if \code{\link{parSapply}} should be used. See Examples.}
  \item{cl}{Cluster created by \code{\link{makeCluster}}. Only used when \code{parallel = TRUE}.}
  \item{verbose}{If \code{TRUE}, tuning progress is printed via \code{\link{txtProgressBar}}. Ignored if \code{parallel = TRUE}.}
  \item{...}{Additional arguments to the \code{\link{cmls}} function, e.g., \code{struc}, \code{backfit}, etc.}
}
  
\details{
The parameters for tuning can be supplied via one of two options: 

(A) Using the \code{parameters} argument. In this case, the argument \code{parameters} must be a data frame with columns \code{const}, \code{df}, \code{degree}, and \code{intercept}, where each row gives a combination of parameters for the CV tuning.

(B) Using the \code{const}, \code{df}, \code{degree}, and \code{intercept} arguments. In this case, the \code{\link{expand.grid}} function is used to create the \code{parameters} data frame, which contains all combinations of the arguments \code{const}, \code{df}, \code{degree}, and \code{intercept}. Duplicates are removed before the CV tuning.
}
\value{
\item{best.parameters}{Best combination of parameters, i.e., the combination that minimizes the \code{cvloss}.}
\item{top5.parameters}{Top five combinations of parameters, i.e., the combinations that give the five smallest values of the \code{cvloss}.}
\item{full.parameters}{Full set of parameters. Data frame with \code{cvloss} (GCV, MSE, or MAE) for each combination of \code{parameters}.}
}
\references{
Helwig, N. E. (in prep). Constrained multivariate least squares in R.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\seealso{
See the \code{\link{cmls}} and \code{\link{const}} functions for further details on the available constraint options.
}

\examples{
# make X
set.seed(1)
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
Ymat <- Xmat \%*\% Bmat + rnorm(n*m, sd = 0.5)


# 5-fold CV:  tune df (5,...,15) for const = "smooth"
kcv <- cv.cmls(X = Xmat, Y = Ymat, nfolds = 5,
               const = "smooth", df = 5:15)
kcv$best.parameters
kcv$top5.parameters
plot(kcv$full.parameters$df, kcv$full.parameters$cvloss, t = "b")


\donttest{

# sample foldid for 5-fold CV
set.seed(2)
foldid <- sample(rep(1:5, length.out = n))


# 5-fold CV:  tune df (5,...,15) w/ all 20 relevant constraints (no struc)
#             using sequential computation (default)
myconst <- as.character(const(print = FALSE)$label[-c(13:16)])
system.time({
  kcv <- cv.cmls(X = Xmat, Y = Ymat, foldid = foldid,
                 const = myconst, df = 5:15)
})
kcv$best.parameters
kcv$top5.parameters


# 5-fold CV:  tune df (5,...,15) w/ all 20 relevant constraints (no struc)
#             using parallel package for parallel computations
myconst <- as.character(const(print = FALSE)$label[-c(13:16)])
system.time({
   cl <- makeCluster(detectCores())
   kcv <- cv.cmls(X = Xmat, Y = Ymat, foldid = foldid,
                  const = myconst, df = 5:15,
                  parallel = TRUE, cl = cl)
   stopCluster(cl)                  
})
kcv$best.parameters
kcv$top5.parameters


# 5-fold CV:  tune df (5,...,15) w/ all 20 relevant constraints (w/ struc)
#             using sequential computation (default)
myconst <- as.character(const(print = FALSE)$label[-c(13:16)])
system.time({
  kcv <- cv.cmls(X = Xmat, Y = Ymat, foldid = foldid,
                 const = myconst, df = 5:15, struc = struc)
})
kcv$best.parameters
kcv$top5.parameters


# 5-fold CV:  tune df (5,...,15) w/ all 20 relevant constraints (w/ struc)
#             using parallel package for parallel computations
myconst <- as.character(const(print = FALSE)$label[-c(13:16)])
system.time({
  cl <- makeCluster(detectCores())
  kcv <- cv.cmls(X = Xmat, Y = Ymat, foldid = foldid,
                 const = myconst, df = 5:15, struc = struc,
                 parallel = TRUE, cl = cl)
  stopCluster(cl)
})
kcv$best.parameters
kcv$top5.parameters

} % end donttest

}

\keyword{ models }
\keyword{ multivariate }
\keyword{ regression }
\keyword{ optimize }
\keyword{ smooth }
