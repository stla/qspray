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
#' @importFrom gmp as.bigq
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
    monomials <- vapply(
      qspray@powers, showMonomial, FUN.VALUE = character(1L)
    )
    coeffs <- as.bigq(qspray@coeffs)
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
#' @description Prints a monomial like \code{"x1.x3^2"}.
#'
#' @param x a string, usually a letter such as \code{"x"} or \code{"X"}, to 
#'   denote the unindexed variables
#'
#' @return A function which takes as argument a sequence of exponents and 
#'   which prints the corresponding monomial.
#' @export
#'
#' @examples
#' showMonomialCanonical("X")(c(1, 0, 2))
showMonomialCanonical <- function(x) {
  function(exponents) {
    paste0(vapply(which(exponents != 0L), function(i) {
      e <- exponents[i]
      if(e == 1L) {
        sprintf("%s%s", x, i)
      } else {
        sprintf("%s%s^%d", x, i, e)
      }
    }, character(1L)), collapse = ".")
  }
}

showMonomialUnivariate <- function(x) {
  function(e) {
    if(length(e) == 0L) {
      ""
    } else if(e == 1L) {
      x
    } else {
      sprintf("%s^%d", x, e)
    }
  }
}

#' @title Print a univariate polynomial
#' @description Prints a polynomial by printing monomials like \code{"x^5"}.
#'
#' @param x a string, usually a letter such as \code{"x"} or \code{"X"}, to 
#'   denote the variable
#'
#' @return A function which prints a univariate \code{qspray} object.
#' @export
showQsprayUnivariate <- function(x) {
  showQspray(showMonomialUnivariate(x = x))
}


#' @title Print a monomial
#' @description Print a monomial like \code{"x^(1, 0, 2)"}.
#'
#' @param x a string, usually a letter such as \code{"x"} or \code{"X"}, to 
#'   denote the variable
#'
#' @return A function which takes as argument a sequence of exponents and 
#'   which prints the corresponding monomial.
#' @export
#'
#' @examples
#' showMonomial("X")(c(1, 0, 2))
showMonomial <- function(x = "x") {
  function(exponents) {
    paste0(sprintf("%s^(", x), exponents, ")") 
  }
}

#' @title Print a 'qspray' object
#' @description Print a \code{qspray} object given a string for the variable.
#'
#' @param x a string, usually a letter such as \code{"x"} or \code{"X"}, 
#'   which denotes the unindexed variables
#' @param ... arguments passed to \code{\link{showQspray}}, such as 
#'   \code{compact=TRUE}
#'
#' @return A function which prints a \code{qspray} object.
#' @export
#' 
#' @note The \code{show} method for \code{qspray} objects uses
#'   \code{showQsprayCanonical("x")} by default.
#'   But this can be controlled as follows. If a
#'   \code{qspray} object has an attribute \code{"x"}, then the value
#'   of this attribute will replace \code{"x"} in the \code{show} output.
#'
#' @examples
#' qspray <- rQspray()
#' showQsprayCanonical("X")(qspray)
showQsprayCanonical <- function(x, ...) {
  showQspray(showMonomialCanonical(x), ...)
}

#' @title Set show option to a 'qspray' object
#' @description Set show option to a \code{qspray} object
#' 
#' @param x a \code{qspray} object
#' @param which which option to set; this can be \code{"x"}, 
#'   \code{"showMonomial"}, or \code{"showQspray"}
#' @param value the value of the option
#'
#' @return This returns the updated \code{qspray}.
#' @export
#'
#' @examples
#' qspray <- rQspray()
#' showQsprayOption(qspray, "x") <- "a"
#' qspray
`showQsprayOption<-` <- function(x, which, value) {
  which <- match.arg(which, c("x", "showMonomial", "showQspray"))
  showOpts <- attr(x, "showOpts") %||% TRUE
  attr(showOpts, which) <- value
  if(which == "x") {
    univariate <- numberOfVariables(x) == 1L
    attr(showOpts, "showQspray") <- if(univariate) {
      showQsprayUnivariate(x = value)
    } else {
      showQsprayCanonical(x = value)
    }
  } else if(which == "showMonomial") {
    attr(showOpts, "showQspray") <- showQspray(showMonomial = value)
  } else {
    attr(showOpts, "showQspray") <- value
  }
  attr(x, "showOpts") <- showOpts
  x
}

getShowQspray <- function(qspray) {
  showOpts <- attr(qspray, "showOpts")
  attr(showOpts, "showQspray") %||%
    attr(attr(showOpts, "showSymbolicQspray"), "showQspray") %||%
    showQsprayCanonical(
      x = attr(showOpts, "x") %||% attr(showOpts, "x") %||% "x",
      qspray = qspray
    )
}
