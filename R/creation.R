#' @title The null qspray polynomial
#' @description Returns the qspray polynomial identically equal to 0.
#' @return A \code{qspray} object.
#' @export
qzero <- function() {
  as.qspray(0L)
}

#' @title The unit qspray polynomial
#' @description Returns the qspray polynomial identically equal to 1.
#' @return A \code{qspray} object.
#' @export
qone <- function() {
  as.qspray(1L)
}

#' @title Polynomial variable
#' @description Create a polynomial variable.
#'
#' @param n nonnegative integer, the index of the variable
#'
#' @return A \code{qspray} object.
#' @export
#' @examples 
#' qlone(2)
qlone <- function(n) {
  stopifnot(isNonnegativeInteger(n))
  powers <- integer(n)
  powers[n] <- 1L
  new("qspray", powers = list(powers), coeffs = "1")
}

#' @title Make a 'qspray' object
#' @description Make a \code{qspray} object from a list of exponents and a 
#'   vector of coefficients.
#'
#' @param powers list of positive integer vectors
#' @param coeffs a vector such that each element of \code{as.character(coeffs)} 
#'   is a quoted integer or a quoted fraction; it must have the same length 
#'   as the \code{powers} list
#' @param string if not \code{NULL}, this argument takes precedence over 
#'   \code{powers} and \code{coeffs}; it must be a string representing a 
#'   multivariate polynomial; see the example
#'
#' @return A \code{qspray} object.
#' @export
#' @examples 
#' powers <- list(c(1, 1), c(0, 2))
#' coeffs <- c("1/2", "4")
#' qsprayMaker(powers, coeffs)
#' qsprayMaker(string = "1/2 x^(1, 1) + 4 x^(0, 2)")
qsprayMaker <- function(powers, coeffs, string = NULL) {
  if(!is.null(string)) {
    List <- stringToQspray(string)
    powers <- List[["powers"]]
    coeffs <- List[["coeffs"]]
  } 
  stopifnot(is.list(powers))
  check_powers <- all(vapply(powers, isExponents, FUN.VALUE = logical(1L)))
  if(!check_powers) {
    stop("Invalid `powers` list.")
  }
  powers <- lapply(powers, as.integer)
  if(!isCoeffs(coeffs)) {
    stop("Invalid `coeffs` vector.")
  }
  if(length(powers) != length(coeffs)) {
    stop("`powers` and `coeffs` must have the same length.")
  }
  qspray_from_list(qspray_maker(powers, as.character(coeffs)))
}

qspray_from_list <- function(qspray_as_list) {
  powers <- qspray_as_list[["powers"]]
  if(is.null(powers)) {
    new(
      "qspray", 
      powers = list(), coeffs = character(0L)
    )
  } else {
    new(
      "qspray", 
      powers = powers, coeffs = qspray_as_list[["coeffs"]]
    )
  }
}

stringToQspray <- function(p){
  stopifnot(isString(p))
  p <- gsub("\\)\\s*-\\s*(\\d*/*\\d*)\\s*", ")+-\\1", p)
  p <- gsub("^-\\s*x", "-1x", trimws(p, "left"))
  terms <- strsplit(p, "+", fixed = TRUE)[[1L]]
  csts <- !grepl("x", terms)
  terms[csts] <- paste0(terms[csts], "x^(0")
  ss <- transpose(strsplit(terms, "x^(", fixed = TRUE))
  coeffs <- trimws(unlist(ss[[1L]], recursive = FALSE))
  coeffs[coeffs == ""] <- "1"
  powers <- sub(")", "", unlist(ss[[2L]], recursive = FALSE), fixed = TRUE)
  powers <- lapply(strsplit(powers, ","), as.integer)
  list(
    "powers" = powers, "coeffs" = coeffs
  )
}
