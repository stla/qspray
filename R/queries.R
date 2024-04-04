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
#' getCoefficient(p, c(2, 0)) # same as getCoefficient(p, 2)
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

#' @title Get the constant term of a qspray polynomial.
#'
#' @param qspray a \code{qspray} object
#'
#' @return A \code{bigq} number.
#' @export
getConstantTerm <- function(qspray) {
  getCoefficient(qspray, integer(0L))
}

#' @title Whether a qspray polynomial is constant
#' @description Checks whether a \code{qspray} object defines a constant 
#'   polynomial.
#'
#' @param qspray a \code{qspray} object
#'
#' @return A Boolean value.
#' @export
isConstantQspray <- function(qspray) {
  arity(qspray) == 0L
}

#' @title Whether a qspray polynomial is null
#' @description Checks whether a \code{qspray} object defines the zero 
#'   polynomial.
#'
#' @param qspray a \code{qspray} object
#'
#' @return A Boolean value.
#' @export
isQzero <- function(qspray) {
  isConstantQspray(qspray) && (getConstantTerm(qspray) == 0L)
}

#' @title Whether a qspray polynomial is the unit polynomial
#' @description Checks whether a \code{qspray} object defines the unit 
#'   polynomial.
#'
#' @param qspray a \code{qspray} object
#'
#' @return A Boolean value.
#' @export
isQone <- function(qspray) {
  isConstantQspray(qspray) && (getConstantTerm(qspray) == 1L)
}

#' @title Whether two qsprays are collinear
#' @description Checks whether two qsprays are collinear, that is, whether 
#'   they are equal up to a scalar factor.
#'
#' @param qspray1,qspray2 two \code{qspray} objects 
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' library(qspray)
#' qspray1 <- qsprayMaker(string = "1/2 x^(1, 1) + 4 x^(0, 2) + 5")
#' qspray2 <- "4/7" * qspray1
#' collinearQsprays(qspray1, qspray2)
collinearQsprays <- function(qspray1, qspray2) {
  if(qspray1 == qzero() && qspray2 == qzero()) {
    return(TRUE)
  }
  if(qspray1 == qzero() && qspray2 != qzero()) {
    return(FALSE)
  }
  if(qspray1 != qzero() && qspray2 == qzero()) {
    return(FALSE)
  }
  M1 <- powersMatrix(qspray1)
  M2 <- powersMatrix(qspray2)
  if(nrow(M1) != nrow(M2) || ncol(M1) != ncol(M2)) {
    return(FALSE)
  }
  ordr1 <- lexorder(M1)
  ordr2 <- lexorder(M2)
  powers1 <- M1[ordr1, , drop = FALSE]
  coeffs1 <- as.bigq(qspray1@coeffs[ordr1])
  powers2 <- M2[ordr2, , drop = FALSE]
  coeffs2 <- as.bigq(qspray2@coeffs[ordr2])
  if(any(powers1 != powers2)) {
    return(FALSE)
  }
  r <- coeffs2[1L] / coeffs1[1L]
  all(r * coeffs1 == coeffs2)
}

