library(qspray)

f <- function(x, y, z) {
  x^2 + 5*y + z - 1
}

x <- qlone(1); y <- qlone(2); z <- qlone(3)

p <- f(x, y, z)


permutation <- c(3, 2, 1)
permutation[permutation] <- seq_along(permutation)

qspray <- p
M <- do.call(rbind, lapply(qspray@powers, qspray:::grow, n = qspray:::arity(qspray)))
M <- M[, permutation]
powers <- apply(M, 1L, identity, simplify = FALSE)

q <- qsprayMaker(powers, p@coeffs)
f(z, x, y)





library(qspray)
f <- function(x, y, z) {
  x^2 + 5*y + z - 1
}
x <- qlone(1)
y <- qlone(2)
z <- qlone(3)
P <- f(x, y, z)
permutation <- c(3, 1, 2)
Q <- permuteVariables(P, permutation)
Q == f(z, x, y) # should be TRUE



