library(qspray)

PSFpoly <- function(m, lambda) {
  #stopifnot(isNonnegativeInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0]
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

PSPpolyExpr <- function(qspray) {
  n <- qspray:::arity(qspray)
  i_ <- seq_len(n)
  P <- lapply(i_, function(i) PSFpoly(n, i))
  Y <- lapply(i_, function(i) qlone(n + i))
  G <- lapply(i_, function(i) P[[i]] - Y[[i]])
  B <- groebner(G, TRUE, FALSE)
  g <- qdivision(qspray, B)
  check <- all(vapply(g@powers, function(pwr) {
    length(pwr) > n && all(pwr[1L:n] == 0L)
  }, logical(1L)))
  if(!check) {
    return(FALSE)
  }
  powers <- lapply(g@powers, function(pwr) {
    pwr[-(1L:n)]
  })
  qsprayMaker(powers, g@coeffs)
}

p <- PSFpoly(3, c(2)) * PSFpoly(3, c(3))^2 + PSFpoly(3, c(1))
PSPpolyExpr(p)
p <- PSFpoly(3, c(3, 3, 2))
q <- PSPpolyExpr(p)
powers <- q@powers
pows1 <- powers[[1L]]
unlist(lapply(rev(seq_along(pows1)), function(i) rep(i, pows1[i])))


#
zlambda <- function(lambda, alpha = 1L) {
  parts <- as.integer(unique(lambda))
  mjs <- vapply(parts, function(j) {
    sum(lambda == j)
  }, integer(1L))
  prod(factorial(mjs) * parts^mjs) * alpha^sum(lambda>0)
}

zlambda(c(3, 2, 2, 1))

#
library(gmp)
HallInnerProduct <- function(spray1, spray2, alpha = 1L) {
  PSspray1 <- PSPpolyExpr(spray1 - getCoefficient(spray1, integer(0L)))
  PSspray2 <- PSPpolyExpr(spray2 - getCoefficient(spray2, integer(0L)))
  n2 <- qspray:::arity(PSspray2)
  powers1 <- PSspray1@powers
  coeffs1 <- PSspray1@coeffs
  out <- as.bigq(0L)
  for(k in seq_along(powers1)) {
    pows <- powers1[[k]]
    coeff2 <- if(length(pows) <= n2) getCoefficient(PSspray2, pows) else 0L
    if(coeff2 != 0L) {
      lambda <- 
        unlist(lapply(rev(seq_along(pows)), function(i) rep(i, pows[i])))
      print(lambda)
      z <- zlambda(lambda, alpha)
      out <- out + as.bigq(coeffs1[k]) * coeff2 * z
    }
  }
  out
}

spray1 <- jack::SchurPol(3, c(3, 1))
spray2 <- jack::SchurPol(4, c(3, 1))
HallInnerProduct(spray1, spray2)

alpha <- as.bigq(2L)
t <- as.integer(alpha)
spray1 <- jack::JackPol(3, c(2, 1), alpha)
spray2 <- jack::JackPol(3, c(2, 1), alpha)
HallInnerProduct(spray1, spray2, alpha = t)
spray1 <- jack::ZonalPol(3, c(2,1))
spray2 <- jack::ZonalPol(3, c(2,1))
HallInnerProduct(spray1, spray2, alpha = t)
t <- alpha
3*t^3/(t^2 + 3/2*t + 1/2)
(2*t^3 + t^2)/(t + 2) 
t^3/6 + t^2/2 + t/3

HallInnerProduct(PSFpoly(4, c(1,1,1,1)), PSFpoly(4, c(1,1,1,1)), alpha = t)
24 * t^4

HallInnerProduct(
  jack::JackPol(3, c(2,1), alpha), PSFpoly(3, c(2,1)), alpha = t
)
t + 2 # no
