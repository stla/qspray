library(qspray)

makeTheMatrix <- function(A) {
  d <- nrow(A)
  M <- list()
  for(i in seq_len(d)) {
    for(j in seq_len(d)) {
      indices <- toString(c(i, j))
      if(i == j) {
        entry <- A[i, j][] - qlone(1)
      } else {
        entry <- A[i, j][]
      }
      M[[indices]] <- entry
    }
  }
  attr(M, "nrow") <- d
  attr(M, "ncol") <- d
  M
}

dropRow <- function(M, i) {
  nrow <- attr(M, "nrow")
  ncol <- attr(M, "ncol")
  Mnew <- list()
  for(k in seq_len(nrow-1L)) {
    for(j in seq_len(ncol)) {
      ii <- if(k < i) k else k+1L
      Mnew[[toString(c(k, j))]] <- M[[toString(c(ii, j))]]
    }
  }
  attr(Mnew, "nrow") <- nrow - 1L
  attr(Mnew, "ncol") <- ncol
  Mnew
}

dropCol <- function(M, j) {
  nrow <- attr(M, "nrow")
  ncol <- attr(M, "ncol")
  Mnew <- list()
  for(i in seq_len(nrow)) {
    for(k in seq_len(ncol - 1L)) {
      jj <- if(k < j) k else k+1L
      Mnew[[toString(c(i, k))]] <- M[[toString(c(i, jj))]]
    }
  }
  attr(Mnew, "nrow") <- nrow
  attr(Mnew, "ncol") <- ncol - 1L
  Mnew
}

minorMatrix <- function(M, i, j) {
  dropCol(dropRow(M, i), j)
}

detLaplace <- function(M) {
  d <- attr(M, "nrow")
  if(d == 1L) {
    M[[toString(c(1L, 1L))]]
  } else {
    result <- 0L
    for(i in seq_len(d)) {
      entry <- M[[toString(c(i, 1L))]]
      if(i %% 2L == 0L) {
        result <- result + entry * detLaplace(minorMatrix(M, i, 1L))
      } else {
        result <- result - entry * detLaplace(minorMatrix(M, i, 1L))
      }
    }
    result
  }
}

#' @importFrom gmp is.matrixZQ
characteristicPolynomial <- function(A) {
  if(!is.matrix(A) && !is.matrixZQ(A)) {
    stop("A must be a matrix.")
  }
  stopifnot(nrow(A) == ncol(A))
  A <- as.bigq(A)
  if(anyNA(A)) {
    stop("Invalid matrix.")
  }
  detLaplace(makeTheMatrix(A))
}

A <- matrix(rpois(9L, 10), nrow = 3L, ncol = 3L)
p <- characteristicPolynomial(A)





