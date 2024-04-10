setGeneric(
  "numberOfVariables", function(x) {
    NULL
  }
)
setGeneric(
  "numberOfTerms", function(qspray) {
    NULL
  }
)
setGeneric(
  "getCoefficient", function(qspray, exponents) {
    NULL
  }
)
setGeneric(
  "getConstantTerm", function(qspray) {
    NULL
  }
)
setGeneric(
  "isConstant", function(x) {
    NULL
  }
)
setGeneric(
  "isUnivariate", function(x) {
    NULL
  }
)
setGeneric(
  "isQzero", function(qspray) {
    NULL
  }
)
setGeneric(
  "isQone", function(qspray) {
    NULL
  }
)

#' @name numberOfVariables
#' @aliases numberOfVariables,qspray-method 
#' @docType methods
#' @title Number of variables in a 'qspray' polynomial
#' @description Number of variables involved in a \code{qspray} object.
#'
#' @param x a \code{qspray} object
#'
#' @return An integer.
#' @export
#' @note The number of variables in the \code{qspray} object \code{qlone(d)} 
#'   is \code{d}, not \code{1}.
setMethod(
  "numberOfVariables", "qspray", 
  function(x) {
    max(0L, arity(x))
  }
)

#' @name numberOfTerms
#' @aliases numberOfTerms,qspray-method 
#' @docType methods
#' @title Number of terms in a 'qspray' polynomial
#' @description Number of terms in the polynomial defined by a 
#'   \code{qspray} object.
#'
#' @param qspray a \code{qspray} object
#'
#' @return An integer.
#' @export
setMethod(
  "numberOfTerms", "qspray", 
  function(qspray) {
    length(qspray@powers)
  }
)

#' @name getCoefficient
#' @aliases getCoefficient,qspray-method 
#' @docType methods
#' @title Get a coefficient in a 'qspray' polynomial
#' @description Get the coefficient corresponding to the given sequence of 
#'   exponents.
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
setMethod(
  "getCoefficient", c("qspray", "numeric"), 
  function(qspray, exponents) {
    stopifnot(isExponents(exponents))
    exponents <- removeTrailingZeros(exponents)
    n <- arity(qspray)
    if(length(exponents) > n) {
      coeff <- 0L
    } else {
      powers <- vapply(qspray@powers, function(pows) {
        toString(grow(pows, n))
      }, character(1L))
      i <- match(toString(grow(exponents, n)), powers)
      if(is.na(i)) {
        coeff <- 0L
      } else {
        coeff <- qspray@coeffs[[i]]
      }
    }
    as.bigq(coeff)
  }
)

#' @name getConstantTerm
#' @aliases getConstantTerm,qspray-method 
#' @docType methods
#' @title Get the constant term of a 'qspray' polynomial
#' @description Get the constant term of a \code{qspray} polynomial.
#'
#' @param qspray a \code{qspray} object
#'
#' @return A \code{bigq} number.
#' @export
setMethod(
  "getConstantTerm", "qspray", 
  function(qspray) {
    getCoefficient(qspray, integer(0L))
  }
)

#' @name isConstant
#' @aliases isConstant,qspray-method 
#' @docType methods
#' @title Whether a 'qspray' polynomial is constant
#' @description Checks whether a \code{qspray} object defines a constant 
#'   polynomial.
#'
#' @param x a \code{qspray} object
#'
#' @return A Boolean value.
#' @export
setMethod(
  "isConstant", "qspray", 
  function(x) {
    numberOfVariables(x) == 0L
  }
)

#' @name isUnivariate
#' @aliases isUnivariate,qspray-method
#' @docType methods
#' @title Whether a 'qspray' is univariate
#' @description Checks whether a \code{qspray} object defines a
#'   univariate polynomial.
#'
#' @param x a \code{qspray} object
#'
#' @return A Boolean value.
#' @export
#' @note It is considered that a constant \code{qspray} is univariate.
setMethod(
  "isUnivariate", "qspray",
  function(x) {
    numberOfVariables(x) %in% c(0L, 1L)
  }
)

#' @name isQzero
#' @aliases isQzero,qspray-method 
#' @docType methods
#' @title Whether a qspray polynomial is null
#' @description Checks whether a \code{qspray} object defines the zero 
#'   polynomial.
#'
#' @param qspray a \code{qspray} object
#'
#' @return A Boolean value.
#' @export
setMethod(
  "isQzero", "qspray", 
  function(qspray) {
    arity(qspray) == -Inf
  }
)

#' @name isQone
#' @aliases isQone,qspray-method 
#' @docType methods
#' @title Whether a qspray polynomial is the unit polynomial
#' @description Checks whether a \code{qspray} object defines the unit 
#'   polynomial.
#'
#' @param qspray a \code{qspray} object
#'
#' @return A Boolean value.
#' @export
setMethod(
  "isQone", "qspray", 
  function(qspray) {
    isConstant(qspray) && getConstantTerm(qspray) == 1L
  }
)

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
