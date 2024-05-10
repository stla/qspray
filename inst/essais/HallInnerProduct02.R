library(qspray)
library(gmp)

spray1 <- jack::SchurPol(3, c(3, 1))
spray2 <- jack::SchurPol(4, c(3, 1))
HallInnerProduct(spray1, spray2)

alpha <- gmp::as.bigq(2L)
spray1 <- jack::JackPol(3, c(2,1), alpha, which = "P")
spray2 <- jack::JackPol(3, c(1,1,1), alpha, which = "P")
HallInnerProduct(spray1, spray1, alpha = alpha)
t <- alpha
3*t^3/(t^2 + 3/2*t + 1/2)
(2*t^3 + t^2)/(t + 2) 
t^3/6 + t^2/2 + t/3

poly <- PSFpoly(4, c(4))
HallInnerProduct(poly, poly, alpha = t)
4*t
poly <- PSFpoly(4, c(3,1))
HallInnerProduct(poly, poly, alpha = t)
3*t^2
poly <- PSFpoly(4, c(2,2))
HallInnerProduct(poly, poly, alpha = t)
8*t^2
poly <- PSFpoly(4, c(2,1,1))
HallInnerProduct(poly, poly, alpha = t)
4*t^3
poly <- PSFpoly(4, c(1,1,1,1))
HallInnerProduct(poly, poly, alpha = t)
24 * t^4

HallInnerProduct(
  jack::JackPol(3, c(2,1), alpha, which = "J"), 
  jack::JackPol(3, c(2,1), alpha, which = "P"), 
  alpha
)
t + 2 # no

# 
library(jack)

#' Skew Jack polynomial
#' @description Computes a skew Jack polynomial.
#' 
#' @param n positive integer, the number of variables
#' @param lambda outer integer partition of the skew partition
#' @param mu inner integer partition of the skew partition; it must be a 
#'   subpartition of \code{lambda}
#' @param alpha the Jack parameter, an integer or a \code{bigq} number, positive
#' @param which which Jack polynomial, \code{"J"}, \code{"P"} or \code{"Q"}
#'
#' @return A \code{qspray} polynomial.
#' @export
#' @importFrom gmp is.bigq as.bigq
#' @importFrom partitions parts
#' @importFrom qspray PSPexpression HallInnerProduct
#' 
#' @details
#' Gröbner bases are used in the algorithm, and this is slow.
#' 
#'
#' @examples
#' SkewJackPol(3, c(3,1), c(2), 2)
SkewJackPol <- function(n, lambda, mu, alpha, which = "J") {
  stopifnot(jack:::isPositiveInteger(n))
  stopifnot(jack:::isPartition(lambda), jack:::isPartition(mu))
  mu <- c(mu, rep(0L, length(lambda) - length(mu)))
  if(any(lambda - mu < 0L)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  stopifnot(jack:::isPositiveInteger(alpha) || is.bigq(alpha))
  alpha <- as.bigq(alpha)
  stopifnot(alpha >= 0L)
  Jlambda <- PSPexpression(JackPolCPP(n, lambda, alpha, which))
  Jmu     <- JackPolCPP(n, mu, alpha, which)
  nus <- LRskew(lambda, mu, "list")$nu
  #nus <- partitions::parts(sum(lambda) - sum(mu))
  terms <- lapply(nus, function(nu) {
    Jnu <- JackPolCPP(n, nu, alpha, which)
    coeff <- HallInnerProduct(Jlambda, Jmu * Jnu, alpha) /
      HallInnerProduct(Jnu, Jnu, alpha) 
    coeff * Jnu
  })
  Reduce(`+`, terms)
}

SkewJackPol(3, c(3,1), c(1), 1L)
SkewSchurPol(3, c(3,1), c(1))


