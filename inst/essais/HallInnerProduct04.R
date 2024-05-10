library(qspray)
library(gmp)
library(partitions)

#' @importFrom partitions compositions
#' @importFrom gmp as.bigq
#' @importFrom DescTools Permn
#' @noRd
E_lambda_mu <- function(lambda, mu) {
  ell_lambda <- length(lambda)
  ell_mu     <- length(mu)
  if(ell_mu > ell_lambda) {
    return(0L)
  }
  # chaque composition donne les longueurs des nu_i 
  compos <- compositions(ell_lambda, ell_mu, include.zero = FALSE)
  compos <- qspray:::Columns(compos)
  lambdas <- DescTools::Permn(lambda)
  L <- do.call(c, lapply(qspray:::Rows(lambdas), function(lambdaPerm) {
    Filter(Negate(is.null), lapply(compos, function(compo) {
      partitionSequences(lambdaPerm, mu, compo)
    }))
  }))
  if(length(L) == 0L) {
    return(as.bigq(0L))
  }
  out <- Reduce(`+`, lapply(L, function(nus) {
      E_lambda_mu_term(mu, nus)
  }))
  if((ell_lambda - ell_mu) %% 2L == 1L) {
    out <- -out
  } 
  return(out)
}

partitionSequences <- function(lambda, mu, compo) {
  starts <- cumsum(c(0L, head(compo, -1L))) + 1L
  ends   <- cumsum(c(0L, head(compo, -1L))) + compo
  nus <- lapply(seq_along(compo), function(i) {
    lambda[(starts[i]):(ends[i])]
  })
  weights <- vapply(nus, function(nu) {
    as.integer(sum(nu))
  }, integer(1L))
  test <- all(mu == weights) && all(vapply(nus, function(nu) {
    all(diff(nu) <= 0L)
  }, logical(1L)))
  if(test) {
    nus
  } else {
    NULL
  }
}

#' @importFrom gmp factorialZ
#' @noRd
E_lambda_mu_term <- function(mu, nus) {
  toMultiply <- lapply(seq_along(nus), function(i) {
    nu  <- nus[[i]]
    mjs <- vapply(as.integer(unique(nu)), function(j) {
      sum(nu == j)
    }, integer(1L))
    mu[i] * factorialZ(length(nu)-1L) / prod(factorialZ(mjs))
  })
  Reduce(`*`, toMultiply)
}

#' @importFrom gmp as.bigq as.bigz c_bigq
#' @importFrom partitions parts
#' @noRd
MSPinPSbasis <- function(mu) {
  mu <- as.integer(mu)
  partitions <- lapply(qspray:::Columns(parts(sum(mu))), qspray:::removeTrailingZeros)
  coeffs <- vector("list", length(partitions))
  k <- 1L
  for(lambda in partitions) {
    coeffs[[k]] <- E_lambda_mu(mu, lambda)
    k <- k + 1L    
  }
  qspray <- qsprayMaker(powers = partitions, coeffs = c_bigq(coeffs))
  lambdas <- qspray@powers
  weights <- as.bigq(qspray@coeffs)
  lapply(seq_along(weights), function(i) {
    lambda <- lambdas[[i]]
    list(
      "coeff"  = weights[i] / as.bigz(zlambda(lambda, alpha = 1L)),
      "lambda" = lambda
    )
  })
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



mu <- c(3L, 3L, 2L, 2L, 1L, 1L)
mu <- c(3L, 2L, 2L)
x <- MSPinPSbasis(mu)
n <- sum(mu)

library(qspray)
check <- qzero()
for(t in x) {
  coeff <- t[["coeff"]]
  if(coeff != 0L) {
    lambda <- t[["lambda"]]
    check <- check + coeff * PSFpoly(n, lambda)
  }
}

check == MSFpoly(n, mu)

# number of terms
choose(n, length(mu)) * nrow(DescTools::Permn(mu))
choose(n, length(mu)) * factorial(length(mu)) / prod(factorial(table(mu)))

