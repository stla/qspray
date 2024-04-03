library(qspray)
library(gmp)
#' @importFrom gmp c_bigq

qspray <- PSFpoly(4, c(3, 1)) + ESFpoly(4, c(2, 2)) # + 3L
qspray <- 4*MSFpoly(3, c(2)) #+ 3*PSFpoly(3, c(3))# + 3L
qspray <- PSFpoly(4, c(2,1,1))
mspdecomposition <- MSPdecomposition(qspray)

pspexpression <- qzero()
for(t in mspdecomposition) {
  x <- qspray:::MSPinPSbasis(t[["lambda"]])
  coeffs <- t[["coeff"]] * c_bigq(sapply(x, `[[`, "coeff", simplify = FALSE))
  powers <- lapply(x, function(xx) {
    lambda <- xx[["lambda"]]
    parts <- as.integer(unique(lambda[lambda != 0L]))
    vapply(1L:14L, function(j) {
      sum(lambda == j)
    }, integer(1L))
  })
  p <- qsprayMaker(powers, coeffs)
  print(p)
  pspexpression <- pspexpression + p
}

pspexpression # oui
#PSPexpression(qspray)
composeQspray(pspexpression, list(PSFpoly(4, 1), PSFpoly(4, 2)))#, PSFpoly(3, 3)))

x <- qlone(1)
y <- qlone(2)
composeQspray(x^2*y, list(PSFpoly(4, 1), PSFpoly(4, 2)))

