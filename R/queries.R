#' @title Get coefficient
#' @description Get the coefficient corresponding to the given exponents.
#' 
#' @param qspray a \code{qspray} object
#' @param exponents a vector of exponents
#'
#' @return The coefficient as a \code{bigq} number.
#' @export
#' @importFrom gmp as.bigq
#'
#' @examples
#' library(qspray)
#' x <- qlone(1)
#' y <- qlone(2)
#' p <- 2*x^2 + 3*y - 5
#' getCoefficient(p, 2)
#' getCoefficient(p, c(2, 0))
#' getCoefficient(p, c(0, 1))
#' getCoefficient(p, 0) # the constant term
#' getCoefficient(p, 3)
getCoefficient <- function(qspray, exponents) {
  stopifnot(isExponents(exponents))
  exponents <- removeTrailingZeros(exponents)
  n <- arity(qspray)
  if(length(exponents) > n) {
    return(as.bigq(0L))
  }
  powers <- vapply(qspray@powers, function(pows) {
    toString(grow(pows, n))
  }, character(1L))
  i <- match(toString(grow(exponents, n)), powers)
  if(is.na(i)) {
    as.bigq(0L)
  } else {
    as.bigq(qspray@coeffs[i])
  }
}

