library(qspray)
library(gmp)
#' @importFrom gmp c_bigq

qspray <- PSFpoly(4, c(3, 1)) + ESFpoly(4, c(2, 2)) # + 3L
mspdecomposition <- MSPdecomposition(qspray)


pspexpression <- qzero()
for(t in mspdecomposition) {
  x <- qspray:::MSPinPSbasis(t[["lambda"]])
  coeffs <- c_bigq(sapply(x, `[[`, "coeff", simplify = FALSE))
  powers <- lapply(x, `[[`, "lambda")
  p <- qsprayMaker(powers, coeffs)
  pspexpression <- pspexpression + p
}







