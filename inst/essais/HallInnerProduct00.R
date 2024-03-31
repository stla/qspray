library(qspray)

PSFpoly <- function(m, lambda) {
  #stopifnot(isNonnegativeInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0]
  if(length(lambda) == 0) {
    return(as.qspray(m))
  }
  if(any(lambda > m)) return(as.qspray(0))
  zeros <- rep(0L, m)
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
zlambda <- function(lambda) {
  parts <- as.integer(unique(lambda))
  mjs <- vapply(parts, function(j) {
    sum(lambda == j)
  }, integer(1L))
  prod(factorial(mjs) * parts^mjs)
}

zlambda(c(3, 2, 2, 1))

#
HallInnerProduct <- function(spray1, spray2) {
  PSspray1 <- PSPpolyExpr(spray1 - getCoefficient(spray1, 0L))
  PSspray2 <- PSPpolyExpr(spray2 - getCoefficient(spray2, 0L))
  
}


