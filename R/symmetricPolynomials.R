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
  lambda <- lambda[lambda != 0L]
  if(length(lambda) > m) return(as.qspray(0L))
  kappa                    <- numeric(m)
  kappa[seq_along(lambda)] <- lambda
  perms <- Permn(kappa)
  n     <- nrow(perms)
  powers <- lapply(1L:n, function(i) perms[i, ])
  coeffs <- rep("1", n)
  qsprayMaker(powers, coeffs)
}

#' @title Symmetric polynomial in terms of the monomial symmetric polynomials
#' @description Decomposition of a symmetric polynomial in the basis formed by 
#'   the monomial symmetric polynomials.
#'
#' @param qspray a \code{qspray} object defining a symmetric polynomial 
#' @param check Boolean, whether to check the symmetry
#'
#' @return A list defining the decomposition. Each element of this list is a 
#'   list with two elements: \code{coeff}, a \code{bigq} number, and 
#'   \code{lambda}, an integer partition; then this list corresponds to the 
#'   term \code{coeff * MSFpoly(n, lambda)}, where \code{n} is the number of 
#'   variables in the symmetric polynomial.
#' @export
#'
#' @examples
#' qspray <- PSFpoly(4, c(3, 1)) + ESFpoly(4, c(2, 2)) + 4L
#' MSPdecomposition(qspray)
MSPdecomposition <- function(qspray, check = TRUE) {
  constantTerm <- getCoefficient(qspray, integer(0L))
  M <- powersMatrix(qspray - constantTerm)
  M <- M[lexorder(M), , drop = FALSE]
  lambdas <- unique(apply(M, 1L, function(expnts) { 
    toString(sort(expnts[expnts != 0L], decreasing = TRUE))
  }))
  out <- lapply(lambdas, function(lambda) {
    lambda <- fromString(lambda)
    list(
      "coeff"  = getCoefficient(qspray, lambda),
      "lambda" = lambda 
    )
  })
  if(constantTerm != 0L) {
    out <- c(
      out, list(list("coeff" = constantTerm, "lambda" = integer(0L)))
    )
  }
  if(check) {
    n <- arity(qspray)
    check <- qzero()
    for(t in out) {
      coeff <- t[["coeff"]]
      lambda <- t[["lambda"]]
      check <- check + coeff * MSFpoly(n, lambda)
    }
    if(check != qspray) {
      stop("The polynomial is not symmetric.")
    }
  }
  out
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

#' @importFrom partitions compositions
#' @noRd
E_lambda_mu <- function(lambda, mu) {
  ell_lambda <- length(lambda)
  ell_mu     <- length(mu)
  if(ell_mu > ell_lambda) {
    return(0L)
  }
  # chaque composition donne les longueurs des nu_i 
  compos <- compositions(ell_lambda, ell_mu, include.zero = FALSE)
  compos <- Columns(compos)
  out <- Reduce(`+`, sapply(compos, function(compo) {
    decoupage(lambda, mu, compo)
  }, simplify = FALSE))
  if((ell_lambda - ell_mu) %% 2L == 0L) {
    out
  } else {
    -out
  }
}

decoupage <- function(lambda, mu, compo) {
  starts <- cumsum(c(0L, head(compo, -1L))) + 1L
  ends   <- cumsum(c(0L, head(compo, -1L))) + compo
  nus <- lapply(seq_along(compo), function(i) {
    lambda[(starts[i]):(ends[i])]
  })
  weights <- vapply(nus, function(nu) {
    as.integer(sum(nu))
  }, integer(1L))
  if(all(weights == mu)) {
    E_lambda_mu_term(mu, nus)
  } else {
    0L
  }
}

#' @importFrom gmp factorialZ
#' @noRd
E_lambda_mu_term <- function(mu, nus) {
  toMultiply <- sapply(seq_along(nus), function(i) {
    nu  <- nus[[i]]
    mjs <- vapply(as.integer(unique(nu)), function(j) {
      sum(nu == j)
    }, integer(1L))
    mu[i] * factorialZ(length(nu)-1L) / prod(factorialZ(mjs))
  }, simplify = FALSE)
  Reduce(`*`, toMultiply)
}

#' @importFrom gmp as.bigz
#' @importFrom partitions parts
#' @noRd
MSPinPSbasis <- function(mu) {
  mu <- as.integer(mu)
  lambdas <- Columns(parts(sum(mu)))
  out <- sapply(lambdas, function(lambda) {
    lambda <- lambda[lambda != 0L]
    E <- E_lambda_mu(mu, lambda)
    if(E != 0L) {
      list(
        "coeff"  = E / as.bigz(zlambda(lambda, alpha = 1L)),
        "lambda" = lambda
      )
    }
  }, simplify = FALSE)
  Filter(Negate(is.null), out)
}

# also used in the Hall inner product
zlambda <- function(lambda, alpha) {
  parts <- as.integer(unique(lambda[lambda != 0L]))
  mjs   <- vapply(parts, function(j) {
    sum(lambda == j)
  }, integer(1L))
  out <- prod(factorial(mjs) * parts^mjs) 
  if(alpha != 1L) {
    out <- out * alpha^length(lambda)
  }
  out
}


#' @title Symmetric polynomial in terms of the power sum polynomials.
#' @description Expression of a symmetric \code{qspray} polynomial as a 
#'   polynomial in the power sum polynomials.
#'
#' @param qspray a symmetric \code{qspray} polynomial (an error is returned if 
#'   it is not symmetric)
#'
#' @return A \code{qspray} polynomial, say \eqn{P}, such that 
#'   \eqn{P(p_1, ..., p_n)} equals the input symmetric polynomial, 
#'   where \eqn{p_i} is the i-th power sum polynomial (\code{PSFpoly(n, i)}).
#' @export
#' 
#' @note
#' This function has been exported because it is used in the \strong{jack} 
#'   package.
PSPexpression <- function(qspray) {
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
  out <- qsprayMaker(powers, g@coeffs) + constantTerm
  attr(out, "PSPexpression") <- TRUE
  out
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
  if(isTRUE(attr(qspray1, "PSPexpression"))) {
    PSspray1 <- qspray1
    PSspray2 <- PSPexpression(qspray2)
  } else {
    PSspray1 <- PSPexpression(qspray1)
    if(qspray2 == qspray1) {
      PSspray2 <- PSspray1
    } else {
      PSspray2 <- PSPexpression(qspray2)
    }
  }
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
