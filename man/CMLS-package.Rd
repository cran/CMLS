\name{CMLS-package}
\alias{CMLS-package}
\alias{CMLS}
\docType{package}
\title{
\packageTitle{CMLS}
}
\description{
\packageDescription{CMLS}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{CMLS}
\packageIndices{CMLS}
The \code{\link{cmls}} function provides a user-friendly interface for solving the MLS problem with 24 common constraint options (the \code{\link{const}} function prints or returns the different contraint options). The \code{\link{cv.cmls}} function does k-fold or generalized cross-validation to tune the constraint options of the \code{\link{cmls}} function. The \code{\link{mlsei}} function solves the MLS problem subject to user-specified equality and/or inequality (E/I) constraints on the coefficients. The \code{\link{mlsun}} function solves the MLS problem subject to unimodality constraints and user-specified E/I constraints on the coefficients.
}
\author{
\packageAuthor{CMLS}

Maintainer: \packageMaintainer{CMLS}
}
\references{
Goldfarb, D., & Idnani, A. (1983). A numerically stable dual method for solving strictly convex quadratic programs. Mathematical Programming, 27, 1-33.

Helwig, N. E. (in prep). Constrained multivariate least squares in R.

Ten Berge, J. M. F. (1993). Least Squares Optimization in Multivariate Analysis. Volume 25 of M & T Series. DSWO Press, Leiden University. ISBN: 9789066950832

Turlach, B. A., & Weingessel, A. (2013). quadprog: Functions to solve Quadratic Programming Problems. R package version 1.5-5. https://CRAN.R-project.org/package=quadprog
}
\keyword{ package }
\examples{
# See examples for cmls, cv.cmls, mlsei, and mlsun
}
