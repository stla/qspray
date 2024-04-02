library(qspray)
library(gmp)
#' @importFrom gmp c_bigq

qspray <- PSFpoly(4, c(3, 1)) + ESFpoly(4, c(2, 2)) # + 3L
qspray <- MSFpoly(3, c(2)) + PSFpoly(3, c(3))# + 3L
mspdecomposition <- MSPdecomposition(qspray)

pspexpression <- qzero()
for(t in mspdecomposition) {
  x <- qspray:::MSPinPSbasis(t[["lambda"]])
  coeffs <- c_bigq(sapply(x, `[[`, "coeff", simplify = FALSE))
  powers <- lapply(x, `[[`, "lambda")
  p <- qsprayMaker(powers, coeffs)
  pspexpression <- pspexpression + p
}

pspexpression # non ...
PSPexpression(qspray)
composeQspray(pspexpression, list(PSFpoly(3, 1)))

#
pspexpression <- qzero()
for(t in mspdecomposition) {
  x <- qspray:::MSPinPSbasis(t[["lambda"]])
  coeffs <- c_bigq(sapply(x, `[[`, "coeff", simplify = FALSE))
  powers <- lapply(x, function(xx) {
    lambda <- xx[["lambda"]]
    parts <- as.integer(unique(lambda[lambda != 0L]))
    vapply(1L:3L, function(j) {
      sum(lambda == j)
    }, integer(1L))
  })
  p <- qsprayMaker(powers, coeffs)
  pspexpression <- pspexpression + p
}

pspexpression # oui
PSPexpression(qspray)
composeQspray(pspexpression, list(PSFpoly(3, 1)))



