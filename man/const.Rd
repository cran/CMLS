\name{const}
\alias{const}
\title{
Print or Return Constraint Options for cmls
}
\description{
Prints or returns six letter constraint codes for \code{\link{cmls}}, along with corresponding descriptions.
}
\usage{
const(x, print = TRUE)
}
\arguments{
  \item{x}{
Vector of six letter constraint codes. If missing, prints/returns all 24 options.
}
  \item{print}{
Should constraint information be printed (\code{print = TRUE}) or returned as a data frame (\code{print = FALSE}).
}
}
\value{
Prints (or returns) constraint codes and descriptions.
}
\references{
Helwig, N. E. (in prep). Constrained multivariate least squares in R.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\seealso{
Constraints are used in the \code{\link{cmls}} function.
}
\examples{
# print some constraints
const(c("uncons", "smpeno"))

# return some constraints
const(c("uncons", "smpeno"), print = FALSE)

# print all constraints
const()

# return all constraints
const(print = FALSE)

}
