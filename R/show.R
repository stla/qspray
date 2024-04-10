#' @title Print a 'qspray' object
#' @description Prints a \code{qspray} object given a function which prints 
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
#' @seealso \code{\link{showQsprayX1X2X3}}, \code{\link{showMonomial}}, 
#'   \code{\link{showQsprayOption<-}}.
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
    plusSign  <- ifelse(compact, "+", " + ")
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
#' @seealso \code{\link{showQsprayX1X2X3}}, 
#'   \code{\link{showMonomialUnivariate}}, \code{\link{showQsprayOption<-}}. 
#' 
#' @examples
#' showMonomialX1X2X3("X")(c(1, 0, 2))
showMonomialX1X2X3 <- function(x) {
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

#' @title Print a monomial
#' @description Prints a monomial like \code{"xz^2"}.
#'
#' @param letters a vector of strings, usually some letters such as \code{"x"} 
#'   and \code{"y"}, to denote the variables
#'
#' @return A function which takes as argument a sequence of exponents and 
#'   which prints the corresponding monomial.
#' @export
#'
#' @seealso \code{\link{showQsprayXYZ}}, 
#'   \code{\link{showMonomialX1X2X3}}, \code{\link{showQsprayOption<-}}. 
#' 
#' @examples
#' showMonomialXYZ()(c(1, 0, 2))
showMonomialXYZ <- function(letters = c("x", "y", "z")) {
  function(exponents) {
    paste0(vapply(which(exponents != 0L), function(i) {
      e <- exponents[i]
      if(e == 1L) {
        letters[i]
      } else {
        sprintf("%s^%d", letters[i], e)
      }
    }, character(1L)), collapse = "")
  }
}

#' @title Print a polynomial
#' @description Prints a polynomial by printing monomials like \code{"x^2yz"}.
#'
#' @param letters a vector of strings, usually some letters such as \code{"x"}
#'   and \code{"y"}, to denote the variables
#' @param ... arguments passed to \code{\link{showQspray}}, such as 
#'   \code{compact=TRUE}
#'
#' @return A function which prints a \code{qspray} object.
#' @export
#' 
#' @seealso \code{\link{showMonomialXYZ}}, \code{\link{showQspray}}, 
#'   \code{\link{showQsprayOption<-}}.
showQsprayXYZ <- function(letters = c("x", "y", "z"), ...) {
  showQspray(showMonomialXYZ(letters), ...)
}

#' @title Print a monomial
#' @description Prints a monomial like \code{"x^(1, 0, 2)"}. This way of 
#'   showing a monomial was used by default in previous versions of this 
#'   package.
#'
#' @param x a string, usually a letter such as \code{"x"} or \code{"X"}, to 
#'   denote the variable
#'
#' @return A function which takes as argument a sequence of exponents and 
#'   which prints the corresponding monomial.
#' @export
#'
#' @seealso \code{\link{showMonomialX1X2X3}}, \code{\link{showQspray}}, 
#'   \code{\link{showQsprayOption<-}}. 
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
#'   which denotes the non-indexed variables
#' @param ... arguments passed to \code{\link{showQspray}}, such as 
#'   \code{compact=TRUE}
#'
#' @return A function which prints a \code{qspray} object.
#' @export
#' 
#' @note The way to print \code{qspray} objects can be controlled with the 
#'  help of the function \code{\link{showQsprayOption<-}}.
#'
#' @examples
#' qspray <- rQspray()
#' showQsprayX1X2X3("X")(qspray)
showQsprayX1X2X3 <- function(x, ...) {
  showQspray(showMonomialX1X2X3(x), ...)
}

#' @title Set a show option to a 'qspray' object
#' @description Set a show option to a \code{qspray} object
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
#' ( qspray <- rQspray() )
#' showQsprayOption(qspray, "x") <- "a"
#' qspray
#' # this is identical to:
#' showQsprayOption(qspray, "showMonomial") <- showMonomialX1X2X3("a")
#' # and also identical to:
#' showQsprayOption(qspray, "showQspray") <- showQsprayX1X2X3("a")
`showQsprayOption<-` <- function(x, which, value) {
  which <- match.arg(which, c("x", "showMonomial", "showQspray"))
  showOpts <- attr(x, "showOpts") %||% TRUE
  attr(showOpts, which) <- value
  if(which == "x") {
    univariate <- isUnivariate(x)
    attr(showOpts, "showQspray") <- if(univariate) {
      showQsprayXYZ(letters = value)
    } else {
      showQsprayX1X2X3(x = value)
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
    showQsprayX1X2X3(
      x = attr(showOpts, "x") %||% "x"
    )
}
