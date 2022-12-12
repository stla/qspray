library(qspray)

Powers1 <- list(c(1L, 1L))
coeffs1 <- "3/1"
Powers2 <- list(c(1L, 1L))
coeffs2 <- "-3/1"

qspray_add(Powers1, coeffs1, Powers2, coeffs2)

qspray_mult(Powers1, coeffs1, Powers2, coeffs2)

qspray_power(Powers1, coeffs1, 0)


isExponents <- function(x) {
  is.numeric(x) && !anyNA(x) && all(x >= 0) && all(floor(x) == x)
}

isCoeffs <- function(x) {
  is.bigq(x) && !anyNA(x)
}

qspray <- function(powers, coeffs) {
  stopifnot(is.list(powers))
  check_powers <- all(vapply(powers, isExponents, FUN.VALUE = logical(1L)))
  if(!check_powers) {
    stop("Invalid `powers` list.")
  }
  powers <- lapply(powers, as.integer)
  if(isCoeffs(coeffs)) {
    stop("Invalid `coeffs` vector.")
  }
  if(length(powers) != length(coeffs)) {
    stop("`powers` and `coeffs` must have the same length.")
  }
  out <- list("powers" = powers, "coeffs" = as.character(coeffs))
  class(out) <- "qspray"
  out
}