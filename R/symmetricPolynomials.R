#' @title Power sum function
#' @description Returns a power sum function as a polynomial.
#'
#' @param m integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#'   positive integers
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' library(qspray)
#' PSFpoly(3, c(3, 1))
PSFpoly <- function(m, lambda) {
  stopifnot(isNonnegativeInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0L]
  if(length(lambda) == 0L) {
    return(as.qspray(m))
  }
  if(any(lambda > m)) return(as.qspray(0L))
  out <- 1L
  for(k in lambda) {
    powers <- lapply(1L:m, function(i) {
      c(rep(0L, i-1L), k)
    })
    pk <- qsprayMaker(powers = powers, coeffs = rep("1", m))
    out <- out * pk
  }
  out
}

#' @title Monomial symmetric function
#' @description Returns a monomial symmetric function as a polynomial.
#'
#' @param m integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#'   positive integers
#'
#' @return A \code{qspray} object.
#' @importFrom DescTools Permn
#' @export
#'
#' @examples
#' library(qspray)
#' MSFpoly(3, c(3, 1))
MSFpoly <- function(m, lambda) {
  stopifnot(isNonnegativeInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(as.qspray(0L))
  kappa <- numeric(m)
  kappa[seq_along(lambda)] <- lambda
  perms <- Permn(kappa)
  n <- nrow(perms)
  powers <- lapply(1L:n, function(i) perms[i, ])
  coeffs <- rep("1", n)
  qsprayMaker(powers, coeffs)
}

#' @title Elementary symmetric polynomial
#' @description Returns an elementary symmetric function as a polynomial.
#'
#' @param m integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#'   positive integers
#'
#' @return A \code{qspray} object.
#' @importFrom DescTools Permn
#' @export
#'
#' @examples
#' library(qspray)
#' ESFpoly(3, c(3, 1))
ESFpoly <- function(m, lambda) {
  stopifnot(isNonnegativeInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0]
  if(any(lambda > m)) return(as.qspray(0L))
  out <- 1L
  for(k in seq_along(lambda)) {
    kappa <- integer(m)
    kappa[seq_len(lambda[k])] <- rep(1L, lambda[k])
    perms <- Permn(kappa)
    n <- nrow(perms)
    powers <- lapply(1L:n, function(i) perms[i, ])
    ek <- qsprayMaker(powers = powers, coeffs = rep("1", n))
    out <- out * ek
  }
  out
}

#' @title Check symmetry of a polynomial
#' @description Check whether a \code{qspray} polynomial is symmetric.
#'
#' @param qspray a \code{qspray} polynomial
#'
#' @return A Boolean value indicating whether the polynomial is symmetric. 
#' In addition, if it is \code{TRUE}, a \code{qspray} polynomial is attached to 
#' the output as an attribute named \code{"poly"}. This polynomial, say \eqn{P},
#' is such that \eqn{P(e_1, ..., e_n)} equals the input symmetric polynomial, 
#' where \eqn{e_i} is the i-th elementary symmetric polynomial 
#' (\code{ESFpoly(n, i)}).
#' @export
#'
#' @examples
#' e1 <- ESFpoly(3, 1)
#' e2 <- ESFpoly(3, 2)
#' e3 <- ESFpoly(3, 3)
#' f <- e1 + 2*e2 + 3*e3 + 4*e1*e3
#' isSymmetricPolynomial(f)
isSymmetricPolynomial <- function(qspray) {
  n <- arity(qspray)
  i_ <- seq_len(n)
  E <- lapply(i_, function(i) ESFpoly(n, i))
  Y <- lapply(i_, function(i) qlone(n + i))
  G <- lapply(i_, function(i) E[[i]] - Y[[i]])
  B <- groebner(G, TRUE, FALSE)
  constantTerm <- getCoefficient(qspray, integer(0L))
  g <- qdivision(qspray - constantTerm, B)
  check <- all(vapply(g@powers, function(pwr) {
    length(pwr) > n && all(pwr[1L:n] == 0L)
  }, logical(1L)))
  if(!check) {
    return(FALSE)
  }
  powers <- lapply(g@powers, function(pwr) {
    pwr[-(1L:n)]
  })
  P <- qsprayMaker(powers, g@coeffs) + constantTerm
  out <- TRUE
  attr(out, "poly") <- P
  out
}


#### ~ Hall inner product ~ ####

# qspray as a polynomial in the power sum polynomials
PSPpolyExpr <- function(qspray) {
  n <- arity(qspray)
  i_ <- seq_len(n)
  P <- lapply(i_, function(i) PSFpoly(n, i))
  Y <- lapply(i_, function(i) qlone(n + i))
  G <- lapply(i_, function(i) P[[i]] - Y[[i]])
  B <- groebner(G, TRUE, FALSE)
  constantTerm <- getCoefficient(qspray, integer(0L))
  g <- qdivision(qspray - constantTerm, B)
  check <- all(vapply(g@powers, function(pwr) {
    length(pwr) > n && all(pwr[1L:n] == 0L)
  }, logical(1L)))
  if(!check) {
    stop("PSPpolyExpr: the polynomial is not symmetric.")
  }
  powers <- lapply(g@powers, function(pwr) {
    pwr[-(1L:n)]
  })
  qsprayMaker(powers, g@coeffs) + constantTerm
}

# helper function for the Hall inner product
zlambda <- function(lambda, alpha) {
  # lambda is clean: no zero
  parts <- as.integer(unique(lambda))
  mjs <- vapply(parts, function(j) {
    sum(lambda == j)
  }, integer(1L))
  prod(factorial(mjs) * parts^mjs) * alpha^length(lambda)
}

#' @title Hall inner product
#' @description Hall inner product of two symmetric polynomials.
#'
#' @param qspray1,qspray2 two symmetric \code{qspray} polynomials
#' @param alpha parameter equal to \code{1} for the usual Hall inner product, 
#'   otherwise this is the "Jack parameter", an integer or a \code{bigq} number
#'
#' @return A \code{bigq} number.
#' @export
#' @importFrom gmp as.bigq is.bigq is.bigz
HallInnerProduct <- function(qspray1, qspray2, alpha = 1) {
  stopifnot(isInteger(alpha) || is.bigq(alpha) || is.bigz(alpha))
  alpha <- as.bigq(alpha)
  PSspray1 <- PSPpolyExpr(qspray1)
  PSspray2 <- PSPpolyExpr(qspray2)
  powers1 <- PSspray1@powers
  coeffs1 <- PSspray1@coeffs
  out <- as.bigq(0L)
  for(k in seq_along(powers1)) {
    pows <- powers1[[k]]
    coeff2 <- getCoefficient(PSspray2, pows)
    if(coeff2 != 0L) {
      lambda <- 
        unlist(lapply(rev(seq_along(pows)), function(i) rep(i, pows[i])))
      out <- out + as.bigq(coeffs1[k]) * coeff2 * zlambda(lambda, alpha)
    }
  }
  out
}
