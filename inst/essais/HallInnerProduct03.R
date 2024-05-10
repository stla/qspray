library(qspray)
library(gmp)
library(partitions)

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
  compos <- qspray:::Columns(compos)
  # compos <- lapply(compos, qspray:::removeTrailingZeros)
  # compos <- Filter(function(compo) length(compo) >= ell_mu, compos)
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

# !!!!!! DOES NOT WORK !!!!!!!!!
#' @importFrom gmp as.bigz
#' @importFrom partitions parts
#' @noRd
MSPinPSbasis <- function(mu) {
  mu <- as.integer(mu)
  
  x <- lapply(mu, parts)
  Grid <- as.matrix(expand.grid(lapply(x, function(partitions) seq_len(ncol(partitions)))))
  powers <- vector("list", nrow(Grid))
  coeffs <- character(nrow(Grid))
  k <- 1L
  for(combo in qspray:::Rows(Grid)) {
    nus <- lapply(seq_along(mu), function(i) {
      qspray:::removeTrailingZeros(x[[i]][, combo[i]])
    })
    lambda <- do.call(c, nus)
    e <- E_lambda_mu(lambda, mu)
    print(e)
    powers[[k]] <- sort(lambda, decreasing = TRUE)
    coeffs[k] <- as.character(e)
    k <- k + 1L
  }
  qspray <- qsprayMaker(powers = powers, coeffs = coeffs)
  print(qspray)
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

MSPinPSbasis <- function(mu) {
  mu <- as.integer(mu)
  lambdas <- qspray:::Columns(parts(sum(mu)))
  coeffs <- vector("list", length(lambdas))
  k <- 1L
  x <- lapply(mu, parts)
  Grid <- as.matrix(expand.grid(lapply(x, function(partitions) seq_len(ncol(partitions)))))
  for(lambda in lambdas) {
    lambda <- qspray:::removeTrailingZeros(lambda)
    coeff <- as.bigq(0L)
    for(combo in qspray:::Rows(Grid)) {
      nus <- lapply(seq_along(mu), function(i) {
        qspray:::removeTrailingZeros(x[[i]][, combo[i]])
      })
      nu <- do.call(c, nus)
      coeff <- coeff + E_lambda_mu(nu, lambda)
    }
    coeffs[[k]] <- coeff
    print(coeff)
    k <- k + 1L    
  }
  qspray <- qsprayMaker(powers = lambdas, coeffs = c_bigq(coeffs))
  print(qspray)
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

MSPinPSbasis <- function(mu) {
  mu <- as.integer(mu)
  x <- lapply(mu, parts)
  Grid <- as.matrix(expand.grid(lapply(x, function(partitions) seq_len(ncol(partitions)))))
  partoches <- lapply(qspray:::Rows(Grid), function(combo) {
    nus <- lapply(seq_along(mu), function(i) {
      qspray:::removeTrailingZeros(x[[i]][, combo[i]])
    })
    do.call(c, nus)
  })
  partoches <- union(partoches, lapply(qspray:::Columns(parts(sum(mu))), qspray:::removeTrailingZeros))
  print(partoches)
  powers <- vector("list", length(partoches))
  coeffs <- character(length(partoches))
  k <- 1L
  for(lambda in partoches) {
    e <- E_lambda_mu(lambda, mu)
    print(e)
    powers[[k]] <- sort(lambda, decreasing = TRUE)
    coeffs[k] <- as.character(e)
    k <- k + 1L
  }
  qspray <- qsprayMaker(powers = powers, coeffs = coeffs)
  print(qspray)
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

MSPinPSbasis <- function(mu) {
  mu <- as.integer(mu)
  k <- 1L
  x <- lapply(mu, parts)
  Grid <- as.matrix(expand.grid(lapply(x, function(partitions) seq_len(ncol(partitions)))))
  partoches <- lapply(qspray:::Rows(Grid), function(combo) {
    nus <- lapply(seq_along(mu), function(i) {
      qspray:::removeTrailingZeros(x[[i]][, combo[i]])
    })
    do.call(c, nus)
  })
  partoches <- union(partoches, lapply(qspray:::Columns(parts(sum(mu))), qspray:::removeTrailingZeros))
  coeffs <- vector("list", length(partoches))
  powers <- vector("list", length(partoches))
  k <- 1L
  for(lambda in partoches) {
    coeff <- E_lambda_mu(mu, lambda)
    coeffs[[k]] <- coeff
    print(coeff)
    powers[[k]] <- sort(lambda, decreasing = TRUE)
    k <- k + 1L    
  }
  qspray <- qsprayMaker(powers = powers, coeffs = c_bigq(coeffs))
  print(qspray)
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



mu <- c(3L, 2L)
x <- MSPinPSbasis(mu)

library(qspray)
check <- qzero()
for(t in x) {
  coeff <- t[["coeff"]]
  if(coeff != 0L) {
    lambda <- t[["lambda"]]
    check <- check + coeff * PSFpoly(5, lambda)
  }
}

check == MSFpoly(5, mu)

pows <- check@powers
mm <- qsprayMaker(pows, coeffs = rep("1L", length(pows)))
