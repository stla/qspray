isString <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
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
  x <- trimws(x)
  if(grepl("^\\-*\\d+$", x)) {
    return(TRUE)
  }
  nd <- trimws(strsplit(x, "/")[[1L]])
  if(length(nd) != 2L) {
    FALSE
  } else {
    n <- nd[1L]
    if(!grepl("^\\-*\\d+$", n)) {
      FALSE
    } else {
      d <- nd[2L]
      if(!grepl("^\\d+$", d) || grepl("^0+$", d)) {
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
  all(vapply(as.character(x), isFraction, FUN.VALUE = logical(1L)))
}

isPartition <- function(lambda){
  length(lambda) == 0L || 
    all(vapply(lambda, isPositiveInteger, FUN.VALUE = logical(1L))) && 
    all(diff(lambda) <= 0)
}

arity <- function(qspray) {
  max(lengths(qspray@powers))
}

removeTrailingZeros <- function(x) {
  n <- length(x)
  while(x[n] == 0 && n > 0L) {
    n <- n - 1L
  }
  head(x, n)
}
