#' @importFrom qspray qsprayDivision as.qspray qone qzero
#' @importFrom resultant gcd
#' @importFrom methods setMethod setClass new show
#' @importFrom gmp as.bigq 
#' @include ratioOfQsprays.R
NULL

setClass(
  "ratioOfQsprays",
  slots = c(numerator = "qspray", denomiator = "qspray")
)

showRatioOfQsprays <- function(roq) {
  if(roq@numerator == 0L) {
    return("0")
  }
  sprintf(
    "[%s] / [%s]", 
    trimws(capture.output(show(roq@numerator)),   which = "right"), 
    trimws(capture.output(show(roq@denominator)), which = "right")
  )
}

setMethod(
  "show", "ratioOfQsprays", 
  function(object) {
    cat(showRatioOfQsprays(object), "\n")
  }
)

identifyQspray <- function(qspray) {
  new("ratioOfQsprays", numerator = qspray, denominator = qone())
}

setGeneric(
  "as.ratioOfQsprays", function(x) {
    NULL
  }
)

#' @name as.ratioOfQsprays
#' @aliases as.ratioOfQsprays,character-method as.ratioOfQsprays,ratioOfQsprays-method as.ratioOfQsprays,qspray-method as.ratioOfQsprays,numeric-method as.ratioOfQsprays,bigz-method as.ratioOfQsprays,bigq-method
#' @exportMethod as.ratioOfQsprays
#' @docType methods
#' @title Coercion to a 'ratioOfQsprays' object
#'
#' @param x a \code{ratioOfQsprays} object, a \code{qspray} object, or an 
#'   object yielding a quoted integer or a quoted fraction after an application 
#'   of \code{as.character}
#'
#' @return A \code{ratioOfQsprays} object.
#' @export
#'
#' @examples
#' library(qspray)
#' as.ratioOfQsprays(2)
#' as.ratioOfQsprays("1/3")
#' as.ratioOfQsprays(5*qlone(1) + qlone(2)^2)
setMethod(
  "as.ratioOfQsprays", "character",
  function(x) {
    identifyQspray(as.qspray(x))
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "ratioOfQsprays",
  function(x) {
    x
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "qspray",
  function(x) {
    identifyQspray(x)
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "numeric",
  function(x) {
    identifyQspray(as.qspray(x))
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "bigz",
  function(x) {
    identifyQspray(as.qspray(x))
  }
)

#' @rdname as.ratioOfQsprays
setMethod(
  "as.ratioOfQsprays", "bigq",
  function(x) {
    identifyQspray(as.qspray(x))
  }
)

#' @name ratioOfQsprays-unary
#' @title Unary operators for ratioOfQsprays objects
#' @description Unary operators for ratioOfQsprays objects.
#' @aliases +,ratioOfQsprays,missing-method -,ratioOfQsprays,missing-method
#' @param e1 object of class \code{ratioOfQsprays}
#' @param e2 nothing
#' @return A \code{ratioOfQsprays} object.
setMethod(
  "+", 
  signature(e1 = "ratioOfQsprays", e2 = "missing"), 
  function(e1, e2) e1
)
#' @rdname ratioOfQsprays-unary
setMethod(
  "-", 
  signature(e1 = "ratioOfQsprays", e2 = "missing"), 
  function(e1, e2) {
    new(
      "ratioOfQsprays", 
      powers = e1@powers, coeffs = as.character(-as.bigq(e1@coeffs))
    )
  }
)

simplifyRatioOfQsprays <- function(roq) {
  num <- roq@numerator
  den <- roq@denominator
  g <- gcd(num, den)
  new(
    "ratioOfQsprays",
    numerator   = qsprayDivision(num, g),
    denominator = qsprayDivision(den, g)
  )
}

ratioOfQsprays_arith_ratioOfQsprays <- function(e1, e2) {
  num1 <- e1@numerator
  den1 <- e1@denominator
  num2 <- e2@numerator
  den2 <- e2@denominator
  switch(
    .Generic,
    "+" = simplifyRatioOfQsprays(
      new(
        "ratioOfQsprays",
        numerator   = num1 * den2 + num2 * den1,
        denominator = den1 * den2
      )
    ),
    "-" = simplifyRatioOfQsprays(
      new(
        "ratioOfQsprays",
        numerator   = num1 * den2 - num2 * den1,
        denominator = den1 * den2
      )
    ),
    "*" = simplifyRatioOfQsprays(
      new(
        "ratioOfQsprays",
        numerator   = num1 * num2,
        denominator = den1 * den2
      )
    ),
    "/" = simplifyRatioOfQsprays(
      new(
        "ratioOfQsprays",
        numerator   = den1 * den2,
        denominator = num1 * num2
      )
    ),
    stop(gettextf(
      "Binary operator %s not defined for ratioOfQsprays objects.", 
      dQuote(.Generic)
    ))
  )
}

ratioOfQsprays_arith_qspray <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.ratioOfQsprays(e2),
    "-" = e1 - as.ratioOfQsprays(e2),
    "*" = e1 * as.ratioOfQsprays(e2),
    "/" = e1 / as.ratioOfQsprays(e2),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

qspray_arith_ratioOfQsprays <- function(e1, e2) {
  ratioOfQsprays_arith_qspray(e2, e1)
}

ratioOfQsprays_arith_character <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.ratioOfQsprays(e2),
    "-" = e1 - as.ratioOfQsprays(e2),
    "*" = e1 * as.ratioOfQsprays(e2),
    "/" = e1 * as.ratioOfQsprays(paste0("1/", e2)),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

ratioOfQspraysPower <- function(ratioOfQsprays, n) {
  if(n == 0L) {
    as.ratioOfQsprays(1L)
  } else if(n > 0L) {
    simplifyRatioOfQsprays(
      new(
        "ratioOfQsprays",
        numerator   = ratioOfQsprays@numerator^n,
        denominator = ratioOfQsprays@denominator
      )
    )
  } else {
    simplifyRatioOfQsprays(
      new(
        "ratioOfQsprays",
        numerator   = ratioOfQsprays@denominator^(-n),
        denominator = ratioOfQsprays@numerator
      )
    )
  }
}

ratioOfQsprays_arith_gmp <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as_ratioOfQsprays(e2),
    "-" = e1 - as_ratioOfQsprays(e2),
    "*" = e1 * as_ratioOfQsprays(e2),
    "/" = e1 * as_ratioOfQsprays(1L/e2),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

ratioOfQsprays_arith_numeric <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.ratioOfQsprays(e2),
    "-" = e1 - as.ratioOfQsprays(e2),
    "*" = e1 * as.ratioOfQsprays(e2),
    "/" = e1 / as.character(e2),
    "^" = ratioOfQspraysPower(e1, e2),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

character_arith_ratioOfQsprays <- function(e1, e2) {
  ratioOfQsprays_arith_character(e2, e1)
}

gmp_arith_ratioOfQsprays <- function(e1, e2) {
  ratioOfQsprays_arith_gmp(e2, e1)
}

numeric_arith_ratioOfQsprays <- function(e1, e2) {
  ratioOfQsprays_arith_numeric(e2, e1)
}

setMethod(
  "Arith", 
  signature(e1 = "ratioOfQsprays", e2 = "ratioOfQsprays"), 
  ratioOfQsprays_arith_ratioOfQsprays
)

setMethod(
  "Arith", 
  signature(e1 = "ratioOfQsprays", e2 = "qspray"), 
  ratioOfQsprays_arith_qspray
)

setMethod(
  "Arith", 
  signature(e1 = "ratioOfQsprays", e2 = "character"), 
  ratioOfQsprays_arith_character
)

setMethod(
  "Arith", 
  signature(e1 = "ratioOfQsprays", e2 = "bigq"), 
  ratioOfQsprays_arith_gmp
)

setMethod(
  "Arith", 
  signature(e1 = "ratioOfQsprays", e2 = "bigz"), 
  ratioOfQsprays_arith_gmp
)

setMethod(
  "Arith", 
  signature(e1 = "qspray", e2 = "ratioOfQsprays"), 
  qspray_arith_ratioOfQsprays
)

setMethod(
  "Arith", 
  signature(e1 = "character", e2 = "ratioOfQsprays"), 
  character_arith_ratioOfQsprays
)

setMethod(
  "Arith", 
  signature(e1 = "bigq", e2 = "ratioOfQsprays"), 
  gmp_arith_ratioOfQsprays
)

setMethod(
  "Arith", 
  signature(e1 = "bigz", e2 = "ratioOfQsprays"), 
  gmp_arith_ratioOfQsprays
)

setMethod(
  "Arith", 
  signature(e1 = "ratioOfQsprays", e2 = "numeric"), 
  ratioOfQsprays_arith_numeric
)

setMethod(
  "Arith", 
  signature(e1 = "numeric", e2 = "ratioOfQsprays"), 
  numeric_arith_ratioOfQsprays
)

setMethod(
  "Compare",
  signature(e1 = "ratioOfQsprays", e2 = "ratioOfQsprays"),
  function(e1, e2) {
    num <- (e1 - e2)@numerator
    switch(
      .Generic,
      "==" = num == qzero(),
      "!=" = num != qzero(),
      stop(gettextf(
        "Comparison operator %s not defined for qspray objects.", dQuote(.Generic)
      ))
    )
  }
)
