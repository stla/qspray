mu <- c(2, 1)
ell <- length(mu)
lambdas <- partitions::parts(sum(mu))

lambda1 <- c(2, 1)
nu1 <- c(2)
nu2 <- c(1)
E_lambda1_mu <- 2
#  1 * (2 * factorial(1) / 1) * (1 * factorial(0) / 1)

lambda2 <- c(1, 1, 1)
nu1 <- c(1, 1)
nu2 <- c(1)

E_lambda2_mu <- -1
#  (-1) * (2 * factorial(1) / factorial(2)) * (1 * factorial(0) / 1)

2 * qspray::PSFpoly(3, lambda1) / qspray:::zlambda(lambda1, 1) + 
  (-1) * qspray::PSFpoly(3, lambda2) / qspray:::zlambda(lambda2, 1)

qspray::MSFpoly(3, mu)

# il fallait calculer E_mu_lambda !
E_mu_lambda1 <- 2
E_mu_lambda2 <- 0  # c'est ce qu'il faudrait - ça m'a l'air bon
E_mu_lambda3 <- -3 # lambda3 = c(3)

############################################################
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
  }))
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

library(gmp)
E_lambda_mu_term <- function(mu, nus) {
  toMultiply <- sapply(seq_along(nus), function(i) {
    nu <- nus[[i]]
    mjs <- vapply(as.integer(unique(nu)), function(j) {
      sum(nu == j)
    }, integer(1L))
    mu[i] * factorialZ(length(nu)-1L) /
      prod(factorialZ(mjs))
  })
  Reduce(`*`, toMultiply)
}

coeffs <- function(mu) {
  lambdas <- partitions::parts(sum(mu))
  z <- function(lambda) {
    parts <- as.integer(unique(lambda[lambda != 0L]))
    mjs <- vapply(parts, function(j) {
      sum(lambda == j)
    }, integer(1L))
    prod(factorial(mjs) * parts^mjs)
  }
  lambdas <- lapply(seq_len(ncol(lambdas)), function(j) {
    lambdas[, j]
  })
  sapply(lambdas, function(lambda) {
    lambda <- lambda[lambda != 0L]
    E_lambda_mu(mu, lambda) / as.bigz(z(lambda))
  })
}


coeffs(c(2L, 1L))





