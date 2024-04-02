library(gmp)

#' @importFrom partitions compositions
#' @noRd
E_lambda_mu <- function(lambda, mu) {
  ell_lambda <- length(lambda)
  ell_mu     <- length(mu)
  if(ell_mu > ell_lambda) {
    return(0L)
  }
  # chaque composition donne les longueurs des nu_i 
  compos <- partitions::compositions(ell_lambda, ell_mu, include.zero = FALSE)
  compos <- lapply(seq_len(ncol(compos)), function(j) {
    compos[, j]
  })
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
    nu <- nus[[i]]
    mjs <- vapply(as.integer(unique(nu)), function(j) {
      sum(nu == j)
    }, integer(1L))
    mu[i] * factorialZ(length(nu)-1L) / prod(factorialZ(mjs))
  })
  Reduce(`*`, toMultiply)
}

#' @importFrom gmp as.bigz
#' @importFrom partitions parts
#' @noRd
MSPinPSbasis <- function(mu) {
  mu <- as.integer(mu)
  lambdas <- partitions::parts(sum(mu))
  z <- function(lambda) {
    parts <- unique(lambda[lambda != 0L])
    mjs <- vapply(parts, function(j) {
      sum(lambda == j)
    }, integer(1L))
    prod(factorial(mjs) * parts^mjs)
  }
  lambdas <- lapply(seq_len(ncol(lambdas)), function(j) {
    lambdas[, j]
  })
  out <- sapply(lambdas, function(lambda) {
    lambda <- lambda[lambda != 0L]
    E <- E_lambda_mu(mu, lambda)
    if(E != 0L) {
      list(
        "coeff"  = E / as.bigz(z(lambda)),
        "lambda" = lambda
      )
    }
  }, simplify = FALSE)
  Filter(Negate(is.null), out)
}


mu <- c(2L, 2L)
x <- MSPinPSbasis(mu)

library(qspray)
check <- qzero()
for(t in x) {
  coeff <- t[["coeff"]]
  if(coeff != 0L) {
    lambda <- t[["lambda"]]
    check <- check + coeff * PSFpoly(4, lambda)
  }
}

check == MSFpoly(4, mu)

################################################################################
# symmetric polynomial in MSP basis
qspray <- PSFpoly(4, c(3, 1)) + ESFpoly(4, c(2, 2)) + 4L
constantTerm <- getCoefficient(qspray, integer(0L))
M    <- qspray:::powersMatrix(qspray - constantTerm)
M    <- M[qspray:::lexorder(M), , drop = FALSE]
lambdas <- unique(apply(M, 1L, function(expnts) { 
  toString(sort(expnts[expnts != 0L], decreasing = TRUE))
}))
fromString <- function(string) {
  as.integer(strsplit(string, ",", fixed = TRUE)[[1L]])
}
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
out






