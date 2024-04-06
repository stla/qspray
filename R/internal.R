isString <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

isStringVector <- function(x) {
  is.character(x) && !anyNA(x)
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

isFractionOrNA <- function(x) {
  is.na(x) || isFraction(x)
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

lexorder <- function(M) {
  do.call(
    order, 
    c(lapply(seq_len(ncol(M)), function(i) M[, i]), decreasing = TRUE)
  )
}

arity <- function(qspray) {
  suppressWarnings(max(lengths(qspray@powers)))
}

#' @importFrom utils head
#' @noRd
removeTrailingZeros <- function(x) {
  n <- length(x)
  while(x[n] == 0 && n > 0L) {
    n <- n - 1L
  }
  head(x, n)
}

isPermutation <- function(x) {
  setequal(x, seq_along(x))
}

fromString <- function(string) {
  as.integer(strsplit(string, ",", fixed = TRUE)[[1L]])
}

Columns <- function(M) {
  lapply(seq_len(ncol(M)), function(j) {
    M[, j]
  })
}

Rows <- function(M) {
  lapply(seq_len(nrow(M)), function(i) {
    M[i, ]
  })
}

grow <- function(powers, n) {
  c(powers, integer(n - length(powers)))
}

lexorder <- function(M){
  do.call(
    function(...) order(..., decreasing = TRUE), 
    lapply(seq_len(ncol(M)), function(i) M[, i])
  )
}
