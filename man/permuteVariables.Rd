% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qspray.R
\name{permuteVariables}
\alias{permuteVariables}
\title{Permute variables}
\usage{
permuteVariables(qspray, permutation)
}
\arguments{
\item{qspray}{a \code{qspray} object}

\item{permutation}{a permutation}
}
\value{
A \code{qspray} object.
}
\description{
Permute the variables of a \code{qspray} polynomial.
}
\examples{
library(qspray)
f <- function(x, y, z) {
  x^2 + 5*y + z - 1
}
x <- qlone(1)
y <- qlone(2)
z <- qlone(3)
P <- f(x, y, z)
permutation <- c(3, 1, 2)
Q <- permuteVariables(P, permutation)
Q == f(z, x, y) # should be TRUE
}
