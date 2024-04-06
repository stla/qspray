#' @title Print a 'qspray' object
#' @description Print a \code{qspray} object given a function which prints 
#'   the monomials.
#'
#' @param showMonomial a function which takes as argument a sequence of 
#'   exponents and which returns a string representing the corresponding 
#'   monomial
#' @param compact a Boolean value
#'
#' @return A function which prints a \code{qspray} object.
#' @export
#' 
#' @seealso \code{\link{showQsprayCanonical}}, \code{\link{showMonomial}} 
showQspray <- function(showMonomial, compact = FALSE) {
  function(qspray) {
    if(length(qspray@coeffs) == 0L) {
      return("0")
    }
    qspray <- orderedQspray(qspray)
    nterms <- length(qspray@coeffs)
    constantTerm <- getConstantTerm(qspray)
    monomials <- vapply(qspray@powers, showMonomial, FUN.VALUE = character(1L))
    coeffs <- gmp::as.bigq(qspray@coeffs)
    plus <- vapply(coeffs, function(x) x >= 0L, FUN.VALUE = logical(1L))
    plusSign <- ifelse(compact, "+", " + ")
    minusSign <- ifelse(compact, "-", " - ")
    signs <- c(ifelse(plus[-1L], plusSign, minusSign), "")
    abscoeffs <- as.character(abs(coeffs))
    terms <- paste0(
      ifelse(abscoeffs == "1", "", paste0(abscoeffs, "*")), monomials
    )
    if(constantTerm != 0L) {
      terms[nterms] <- as.character(abs(constantTerm))
    }
    leader <- if(plus[1L]) "" else "-"
    paste0(c(leader, c(rbind(terms, signs))), collapse = "")
  }
}

#' @title Print a monomial
#' @description Print a monomial like \code{"x1.x3^2"}.
#'
#' @param var a string, usually a letter such as \code{"x"} or \code{"X"}, to 
#'   denote the unindexed variables
#'
#' @return A function which takes as argument a sequence of exponents and 
#'   which prints the corresponding monomial.
#' @export
#'
#' @examples
#' showMonomialCanonical("X")(c(1, 0, 2))
showMonomialCanonical <- function(var) {
  function(exponents) {
    exponents <- exponents[exponents != 0L]
    paste0(vapply(seq_along(exponents), function(i) {
      e <- exponents[i]
      if(e == 1L) {
        sprintf("%s%d", var, i)
      } else {
        sprintf("%s%d^%d", var, i, e)
      }
    }, character(1L)), collapse = ".")
  }
}

#' @title Print a monomial
#' @description Print a monomial like \code{"x^(1, 0, 2)"}.
#'
#' @param var a string, usually a letter such as \code{"x"} or \code{"X"}, to 
#'   denote the variable
#'
#' @return A function which takes as argument a sequence of exponents and 
#'   which prints the corresponding monomial.
#' @export
#'
#' @examples
#' showMonomial("X")(c(1, 0, 2))
showMonomial <- function(var = "x") {
  function(exponents) {
    paste0(sprintf("%s^(", var), exponents, ")") 
  }
}

#' @title Print a 'qspray' object
#' @description Print a \code{qspray} object given a string for the variable.
#'
#' @param var a string, usually a letter such as \code{"x"} or \code{"X"}, 
#'   which denotes the unindexed variables
#' @param ... arguments passed to \code{\link{showQspray}}, such as 
#'   \code{compact=TRUE}
#'
#' @return A function which prints a \code{qspray} object.
#' @export
#'
#' @examples
#' qspray <- rQspray()
#' showQsprayCanonical("X")(qspray)
showQsprayCanonical <- function(var, ...) {
  showQspray(showMonomialCanonical(var), ...)
}
