isExponents <- function(x) {
  is.numeric(x) && !anyNA(x) && all(floor(x) == x)
}

isCoeffs <- function(x) {
  is.bigq(x) && !anyNA(x)
}

isPositiveInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x
}
