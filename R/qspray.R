#' @useDynLib qspray, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods setMethod setClass new canCoerce as
#' @importFrom gmp as.bigq factorialZ asNumeric
#' @importFrom purrr transpose
#' @include qspray.R
NULL

# if(getRversion() >= "2.15.1") {
#   globalVariables("as.ratioOfQsprays")
# }

setClass(
  "qspray",
  slots = c(powers = "list", coeffs = "character")
)

setMethod(
  "show", "qspray", 
  function(object) {
    if(is.null(attr(object, "showOpts"))) {
      trivariate <- numberOfVariables(object) <= 3
      showQsprayOption(object, "showQspray") <- 
        if(trivariate) showQsprayXYZ() else showQsprayX1X2X3("x")
    }
    f <- getShowQspray(object)#attr(attr(object, "showOpts"), "showQspray")
    cat(f(object), "\n")
  }
)

setGeneric(
  "showCoefficient", function(x) {
    NULL
  }
)

#' @name showCoefficient
#' @aliases showCoefficient,qspray-method 
#' @docType methods
#' @title Function which prints the coefficients
#' @description This method is for internal usage. For any \code{qspray} 
#'   object, it returns a function which takes as argument a coefficient 
#'   and which returns the string formed by this coefficient surrounded by
#'   parentheses.
#'
#' @param x a \code{qspray} object 
#'
#' @return A function which associates a character string to the coefficients.
#' @export
setMethod(
  "showCoefficient", "qspray", 
  function(x) {
    function(coeff) {
      sprintf("(%s)", as.character(coeff))
    }
  }
)


as.qspray.character <- function(x) {
  stopifnot(isFraction(x))
  if(as.bigq(x) == 0L) {
    new("qspray", powers = list(), coeffs = character(0L))
  } else {
    new("qspray", powers = list(integer(0L)), coeffs = x)
  }
}

as_qspray_gmp <- function(x) {
  if(x == 0L) {
    new("qspray", powers = list(), coeffs = character(0L))
  } else {
    new("qspray", powers = list(integer(0L)), coeffs = as.character(x))
  }
}

as.qspray.numeric <- function(x) {
  stopifnot(isInteger(x))
  if(x == 0L) {
    new("qspray", powers = list(), coeffs = character(0L))
  } else {
    new("qspray", powers = list(integer(0L)), coeffs = as.character(x))
  }
}

setGeneric(
  "as.qspray", function(x) {
    NULL
  }
)

#' @name as.qspray
#' @aliases as.qspray,character-method as.qspray,qspray-method as.qspray,numeric-method as.qspray,bigz-method as.qspray,bigq-method
#' @exportMethod as.qspray
#' @docType methods
#' @title Coercion to a 'qspray' object
#'
#' @param x a \code{qspray} object or an object yielding a quoted integer or a 
#'   quoted fraction after an application of \code{as.character}
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' as.qspray(2)
#' as.qspray("1/3")
setMethod(
  "as.qspray", "character",
  function(x) {
    as.qspray.character(x)
  }
)

#' @rdname as.qspray
setMethod(
  "as.qspray", "qspray",
  function(x) {
    x
  }
)

#' @rdname as.qspray
setMethod(
  "as.qspray", "numeric",
  function(x) {
    as.qspray.numeric(x)
  }
)

#' @rdname as.qspray
setMethod(
  "as.qspray", "bigz",
  function(x) {
    as_qspray_gmp(x)
  }
)

#' @rdname as.qspray
setMethod(
  "as.qspray", "bigq",
  function(x) {
    as_qspray_gmp(x)
  }
)

#' @name qspray-unary
#' @title Unary operators for qspray objects
#' @description Unary operators for qspray objects.
#' @aliases +,qspray,missing-method -,qspray,missing-method
#' @param e1 object of class \code{qspray}
#' @param e2 nothing
#' @return A \code{qspray} object.
setMethod(
  "+", 
  signature(e1 = "qspray", e2 = "missing"), 
  function(e1, e2) e1
)
#' @rdname qspray-unary
setMethod(
  "-", 
  signature(e1 = "qspray", e2 = "missing"), 
  function(e1, e2) {
    qspray <- new(
      "qspray", 
      powers = e1@powers, coeffs = as.character(-as.bigq(e1@coeffs))
    )
    passShowAttributes(e1, qspray)
  }
)


qspray_arith_character <- function(e1, e2) {
  qspray <- switch(
    .Generic,
    "+" = e1 + as.qspray.character(e2),
    "-" = e1 - as.qspray.character(e2),
    "*" = e1 * as.qspray.character(e2),
    "/" = e1 * as.qspray.character(paste0("1/", e2)),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
  passShowAttributes(e1, qspray)
}

#' @importFrom utils installed.packages
qspray_arith_qspray <- function(e1, e2) {
  qspray <- switch(
    .Generic,
    "+" = qspray_from_list(
      qspray_add(e1@powers, e1@coeffs, e2@powers, e2@coeffs)
    ),
    "-" = qspray_from_list(
      qspray_subtract(e1@powers, e1@coeffs, e2@powers, e2@coeffs)
    ),
    "*" = qspray_from_list(
      qspray_mult(e1@powers, e1@coeffs, e2@powers, e2@coeffs)
    ),
    "/" = {
      if(canCoerce(e1, "ratioOfQsprays")) {
        as(e1, "ratioOfQsprays") / as(e2, "ratioOfQsprays") 
      } else {
        x <- "loaded"
        stop(
          "Division of 'qspray' objects is possible only with the ",
          "'ratioOfQsprays' package, and this package is not ", x, "."
        )
      }
    },
    # "/" = tryCatch({
    #     as.ratioOfQsprays(e1) / as.ratioOfQsprays(e2)
    # }, error = function(e) {
    #   x <- ifelse(
    #     "ratioOfQsprays" %in% rownames(installed.packages()),
    #     "loaded",
    #     "installed"
    #   )
    #   stop(
    #     "Division of 'qspray' objects is possible only with the ",
    #     "'ratioOfQsprays' package, and this package is not ", x, "."
    #   )
    # }),
    stop(gettextf(
      "Binary operator %s not defined for qspray objects.", dQuote(.Generic)
    ))
  )
  passShowAttributes(e2, qspray)
}

qsprayPower <- function(e1, n) {
  stopifnot(isPositiveInteger(n))
  qspray <- qspray_power(e1@powers, e1@coeffs, n)
  passShowAttributes(e1, qspray)
}

qspray_arith_gmp <- function(e1, e2) {
  qspray <- switch(
    .Generic,
    "+" = e1 + as_qspray_gmp(e2),
    "-" = e1 - as_qspray_gmp(e2),
    "*" = e1 * as_qspray_gmp(e2),
    "/" = e1 * as_qspray_gmp(1L/e2),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
  passShowAttributes(e1, qspray)
  qspray
}

qspray_arith_numeric <- function(e1, e2) {
  qspray <- switch(
    .Generic,
    "+" = e1 + as.qspray.numeric(e2),
    "-" = e1 - as.qspray.numeric(e2),
    "*" = e1 * as.qspray.numeric(e2),
    "/" = e1 / as.character(e2),
    "^" = qspray_from_list(qsprayPower(e1, e2)),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
  passShowAttributes(e1, qspray)
}

character_arith_qspray <- function(e1, e2) {
  qspray <- switch(
    .Generic,
    "+" = as.qspray.character(e1) + e2,
    "-" = as.qspray.character(e1) - e2,
    "*" = as.qspray.character(e1) * e2,
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
  passShowAttributes(e2, qspray)
}

gmp_arith_qspray <- function(e1, e2) {
  qspray <- switch(
    .Generic,
    "+" = as_qspray_gmp(e1) + e2,
    "-" = as_qspray_gmp(e1) - e2,
    "*" = as_qspray_gmp(e1) * e2,
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
  passShowAttributes(e2, qspray)
}

numeric_arith_qspray <- function(e1, e2) {
  qspray <- switch(
    .Generic,
    "+" = as.qspray.numeric(e1) + e2,
    "-" = as.qspray.numeric(e1) - e2,
    "*" = as.qspray.numeric(e1) * e2,
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
  passShowAttributes(e2, qspray)
}

setMethod(
  "Arith", 
  signature(e1 = "qspray", e2 = "qspray"), 
  qspray_arith_qspray
)

setMethod(
  "Arith", 
  signature(e1 = "qspray", e2 = "character"), 
  qspray_arith_character
)

setMethod(
  "Arith", 
  signature(e1 = "qspray", e2 = "bigq"), 
  qspray_arith_gmp
)

setMethod(
  "Arith", 
  signature(e1 = "qspray", e2 = "bigz"), 
  qspray_arith_gmp
)

setMethod(
  "Arith", 
  signature(e1 = "character", e2 = "qspray"), 
  character_arith_qspray
)

setMethod(
  "Arith", 
  signature(e1 = "bigq", e2 = "qspray"), 
  gmp_arith_qspray
)

setMethod(
  "Arith", 
  signature(e1 = "bigz", e2 = "qspray"), 
  gmp_arith_qspray
)

setMethod(
  "Arith", 
  signature(e1 = "qspray", e2 = "numeric"), 
  qspray_arith_numeric
)

setMethod(
  "Arith", 
  signature(e1 = "numeric", e2 = "qspray"), 
  numeric_arith_qspray
)


setMethod(
  "Compare",
  signature(e1 = "qspray", e2 = "qspray"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = qspray_equality(e1@powers, e1@coeffs, e2@powers, e2@coeffs),
      "!=" = !qspray_equality(e1@powers, e1@coeffs, e2@powers, e2@coeffs),
      stop(gettextf(
        "Comparison operator %s not defined for qspray objects.", dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "qspray", e2 = "character"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.qspray(e2),
      "!=" = e1 != as.qspray(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.", 
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "qspray", e2 = "numeric"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.qspray(e2),
      "!=" = e1 != as.qspray(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.", 
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "qspray", e2 = "bigz"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.qspray(e2),
      "!=" = e1 != as.qspray(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.", 
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "qspray", e2 = "bigq"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = e1 == as.qspray(e2),
      "!=" = e1 != as.qspray(e2),
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.", 
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "character", e2 = "qspray"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.qspray(e1) == e2,
      "!=" = as.qspray(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.", 
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "numeric", e2 = "qspray"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.qspray(e1) == e2,
      "!=" = as.qspray(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.", 
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "bigz", e2 = "qspray"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.qspray(e1) == e2,
      "!=" = as.qspray(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.", 
        dQuote(.Generic)
      ))
    )
  }
)
setMethod(
  "Compare",
  signature(e1 = "bigq", e2 = "qspray"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = as.qspray(e1) == e2,
      "!=" = as.qspray(e1) != e2,
      stop(gettextf(
        "Comparison operator %s not defined for these two objects.", 
        dQuote(.Generic)
      ))
    )
  }
)
