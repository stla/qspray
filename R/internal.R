isExponents <- function(x) {
  is.numeric(x) && !anyNA(x) && all(floor(x) == x)
}

isCoeffs <- function(x) {
  is.bigq(x) && !anyNA(x)
}

isInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && as.integer(x) == x
}

isPositiveInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x
}

isNonnegativeInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x && x != 0
}

isFraction <- function(x) {
  if(!is.character(x) || length(x) != 1L) {
    return(FALSE)
  }
  nd <- strsplit(x, "/")[[1L]]
  if(nd != 2L) {
    FALSE
  } else {
    n <- nd[1L]
    if(!grepl("^\\-*\\d+$")) {
      FALSE
    } else {
      d <- nd[2L]
      if(!grepl("^\\d+$") || grepl("^0+$")) {
        FALSE
      } else {
        TRUE
      }
    }
  }
}

isExponents <- function(x) {
  is.numeric(x) && !anyNA(x) && all(floor(x) == x)
}

isCoeffs <- function(x) {
  all(vapply(x, isFraction, FUN.VALUE = logical(1L)))
}
