library(qspray)

MSFpoly <- function(m, lambda){
  #stopifnot(m > 0L, isPositiveInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(as.qspray(0))
  kappa <- numeric(m)
  kappa[seq_along(lambda)] <- lambda
  perms <- DescTools::Permn(kappa)
  n <- nrow(perms)
  powers <- lapply(1L:n, function(i) perms[i, ])
  coeffs <- rep("1", n)
  qsprayMaker(powers, coeffs)
}


MSFpoly(3, c(3,1))
