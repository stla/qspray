#' @include qspray.R
NULL

#' @title Ordered 'qspray'
#' @description Reorders the terms of a \code{qspray} object according to the 
#'   lexicographic order of the powers. 
#'
#' @param qspray a \code{qspray} object
#'
#' @return A \code{qspray} object. It defines the same polynomial as the 
#'   input \code{qspray} object but it is ordered.
#' @export
#'
#' @examples
#' qspray <- rQspray()
#' qspray == orderedQspray(qspray) # should be TRUE
orderedQspray <- function(qspray) {
  M <- powersMatrix(qspray)
  if(ncol(M) > 0L) {
    lex <- lexorder(M)
    qspray@powers <- qspray@powers[lex]
    qspray@coeffs <- qspray@coeffs[lex]
  }
  qspray
}

#' @title Partial derivative
#' @description Partial derivative of a qspray polynomial.
#'
#' @param qspray object of class \code{qspray}
#' @param i integer, the dimension to differentiate with respect to
#' @param derivative integer, how many times to differentiate
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' library(qspray)
#' x <- qlone(1)
#' y <- qlone(2)
#' qspray <- 2*x  + 3*x*y
#' derivQspray(qspray, 1)
derivQspray <- function(qspray, i, derivative = 1) {
  stopifnot(inherits(qspray, "qspray"))
  stopifnot(isNonnegativeInteger(i))
  stopifnot(isPositiveInteger(derivative))
  if(i > arity(qspray)) {
    dqspray <- qzero()
  } else {
    n    <- integer(length = i)
    n[i] <- as.integer(derivative)
    drv  <- qspray_deriv(qspray@powers, qspray@coeffs, n)
    dqspray <- qspray_from_list(drv)
  }
  passShowAttributes(qspray, dqspray)
}

#' @title Partial differentiation
#' @description Partial differentiation of a qspray polynomial.
#'
#' @param qspray object of class \code{qspray}
#' @param orders integer vector, the orders of the differentiation
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' library(qspray)
#' x <- qlone(1)
#' y <- qlone(2)
#' qspray <- x + 2*y  + 3*x*y
#' dQspray(qspray, c(1, 1))
#' derivQspray(derivQspray(qspray, 1), 2)
dQspray <- function(qspray, orders) {
  stopifnot(inherits(qspray, "qspray"))
  for(i in seq_along(orders)) {
    stopifnot(isPositiveInteger(orders[i]))
  }
  orders <- removeTrailingZeros(orders)
  if(length(orders) > arity(qspray)) {
    dqspray <- qzero()
  } else {
    n    <- as.integer(orders)
    drv  <- qspray_deriv(qspray@powers, qspray@coeffs, n)
    dqspray <- qspray_from_list(drv)
  }
  passShowAttributes(qspray, dqspray)
}

setGeneric(
  "permuteVariables", function(x, permutation) {
    NULL
  }
)

#' @name permuteVariables
#' @aliases permuteVariables,qspray,numeric-method 
#' @docType methods
#' @title Permute variables
#' @description Permute the variables of a \code{qspray} polynomial.
#'
#' @param x a \code{qspray} object
#' @param permutation a permutation
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' library(qspray)
#' f <- function(x, y, z) {
#'   x^2 + 5*y + z - 1
#' }
#' x <- qlone(1)
#' y <- qlone(2)
#' z <- qlone(3)
#' P <- f(x, y, z)
#' permutation <- c(3, 1, 2)
#' Q <- permuteVariables(P, permutation)
#' Q == f(z, x, y) # should be TRUE
setMethod(
  "permuteVariables", c("qspray", "numeric"), 
  function(x, permutation) {
    stopifnot(isPermutation(permutation))
    m <- arity(x)
    n <- length(permutation)
    if(m > n) {
      stop("Invalid permutation.")
    }
    permutation[permutation] <- seq_along(permutation)
    M <- powersMatrix(x)
    for(. in seq_len(n - m)) {
      M <- cbind(M, 0L)
    }
    M <- M[, permutation]
    powers <- apply(M, 1L, function(row) {
      removeTrailingZeros(row)
    }, simplify = FALSE)
    out <- new("qspray", powers = powers, coeffs = x@coeffs)
    passShowAttributes(x, out)
  }
)

setGeneric(
  "swapVariables", function(x, i, j) {
    NULL
  }
)

#' @name swapVariables
#' @aliases swapVariables,qspray,numeric,numeric-method 
#' @docType methods
#' @title Swap variables
#' @description Swap two variables in a \code{qspray} polynomial.
#'
#' @param x a \code{qspray} object
#' @param i,j indices of the variables to be swapped
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' library(qspray)
#' f <- function(x, y, z) {
#'   x^2 + 5*y + z - 1
#' }
#' x <- qlone(1)
#' y <- qlone(2)
#' z <- qlone(3)
#' P <- f(x, y, z)
#' Q <- swapVariables(P, 2, 3)
#' Q == f(x, z, y) # should be TRUE
setMethod(
  "swapVariables", c("qspray", "numeric", "numeric"), 
  function(x, i, j) {
    stopifnot(isNonnegativeInteger(i), isNonnegativeInteger(j))
    m <- arity(x)
    n <- max(m, i, j)
    permutation <- seq_len(n)
    permutation[i] <- j
    permutation[j] <- i
    M <- powersMatrix(x)
    for(. in seq_len(n - m)) {
      M <- cbind(M, 0L)
    }
    M <- M[, permutation]
    powers <- apply(M, 1L, function(row) {
      removeTrailingZeros(row)
    }, simplify = FALSE)
    out <- new("qspray", powers = powers, coeffs = x@coeffs)
    passShowAttributes(x, out)
  }
)