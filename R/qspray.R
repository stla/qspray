#' @useDynLib qspray, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods setMethod setClass new
#' @importFrom gmp as.bigq factorialZ
#' @include qspray.R
NULL

#' Title
#'
#' @param M xx
#'
#' @return xx
#' @export
detQ <- function(M) {
  detQ_rcpp(M)
}

setClass(
  "qspray",
  slots = c(powers = "list", coeffs = "character")
)

showQspray <- function(qspray) {
  powers <- vapply(qspray@powers, toString, FUN.VALUE = character(1L))
  coeffs <- as.bigq(qspray@coeffs)
  plus <- vapply(coeffs, function(x) x >= 0, FUN.VALUE = logical(1L))
  signs <- c(ifelse(plus[-1L], " + ", " - "), "")
  abscoeffs <- as.character(abs(coeffs))
  terms <- paste0(
    ifelse(abscoeffs == "1", "", paste0(abscoeffs, "*")), 
    "x^(", powers, ")"
  )
  leader <- if(plus[1L]) "" else "-"
  paste0(c(leader, c(rbind(terms, signs))), collapse = "")
}

setMethod(
  "show", "qspray", 
  function(object) {
    cat(showQspray(object), "\n")
  }
)

qspray_from_list <- function(qspray_as_list) {
  powers <- qspray_as_list[["powers"]]
  if(is.null(powers)) {
    new(
      "qspray", 
      powers = list(integer(0L)), coeffs = "0"
    )
  } else {
    new(
      "qspray", 
      powers = powers, coeffs = qspray_as_list[["coeffs"]]
    )
  }
}

#' Title
#'
#' @param powers ss
#' @param coeffs ss
#'
#'
#' @return xx
#' @export
qsprayMaker <- function(powers, coeffs) {
  stopifnot(is.list(powers))
  check_powers <- all(vapply(powers, isExponents, FUN.VALUE = logical(1L)))
  if(!check_powers) {
    stop("Invalid `powers` list.")
  }
  powers <- lapply(powers, as.integer)
  if(isCoeffs(coeffs)) {
    stop("Invalid `coeffs` vector.")
  }
  if(length(powers) != length(coeffs)) {
    stop("`powers` and `coeffs` must have the same length.")
  }
  qspray_from_list(qspray_maker(powers, as.character(coeffs)))
}

#' Title
#'
#' @param n xx
#'
#' @return xx
#' @export
lone <- function(n) {
  stopifnot(isNonnegativeInteger(n))
  powers <- integer(n)
  powers[n] <- 1L
  new("qspray", powers = list(powers), coeffs = "1")
}


#' Title
#'
#' @param qspray xx
#' @param values xx
#'
#' @return xx
#' @export
evalQspray <- function(qspray, values) {
  powers <- qspray@powers
  coeffs <- as.bigq(qspray@coeffs)
  values <- as.bigq(values)
  out <- as.bigq(0L)
  for(i in seq_along(powers)) {
    exponents <- powers[[i]]
    term <- as.bigq(1L)
    for(j in seq_along(exponents)) {
      term <- term * values[j]^exponents[j]
    }
    out <- out + coeffs[i] * term
  }
  out
}


as.qspray.character <- function(x) {
  stopifnot(isFraction(x))
  new("qspray", powers = list(integer(0L)), coeffs = x)
}

as_qspray_gmp <- function(x) {
  new("qspray", powers = list(integer(0L)), coeffs = as.character(x))
}

as.qspray.numeric <- function(x) {
  stopifnot(isInteger(x))
  new("qspray", powers = list(integer(0L)), coeffs = as.character(x))
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
#' @title Coercion to a qspray
#'
#' @param x a \code{qspray} object or an object yielding an integer or a 
#'   fraction after an application of \code{as.character}
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' x
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
    new("qspray", powers = e1@powers, coeffs = as.character(as.bigq(e1@coeffs)))
  }
)

qspray_arith_qspray <- function(e1, e2) {
  switch(
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
    stop(gettextf(
      "Binary operator %s not defined for lazy vectors.", dQuote(.Generic)
    ))
  )
}

qsprayPower <- function(qspray, n) {
  stopifnot(isPositiveInteger(n))
  qspray_power(qspray@powers, qspray@coeffs, n)
}

qspray_arith_character <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.qspray.character(e2),
    "-" = e1 - as.qspray.character(e2),
    "*" = e1 * as.qspray.character(e2),
    "/" = e1 * as.qspray.character(paste0("1/", e2)),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

qspray_arith_gmp <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as_qspray_gmp(e2),
    "-" = e1 - as_qspray_gmp(e2),
    "*" = e1 * as_qspray_gmp(e2),
    "/" = e1 / as.character(e2),
    stop(gettextf(
      "Binary operator %s not defined for qspray objects.", dQuote(.Generic)
    ))
  )
}

qspray_arith_numeric <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.qspray.numeric(e2),
    "-" = e1 - as.qspray.numeric(e2),
    "*" = e1 * as.qspray.numeric(e2),
    "/" = e1 / as.character(e2),
    "^" = qspray_from_list(qsprayPower(e1, e2)),
    stop(gettextf(
      "Binary operator %s not defined for qspray objects.", dQuote(.Generic)
    ))
  )
}

character_arith_qspray <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as.qspray.character(e1) + e2,
    "-" = as.qspray.character(e1) - e2,
    "*" = as.qspray.character(e1) * e2,
    stop(gettextf(
      "Binary operator %s not defined for qspray objects.", dQuote(.Generic)
    ))
  )
}

gmp_arith_qspray <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as_qspray_gmp(e1) + e2,
    "-" = as_qspray_gmp(e1) - e2,
    "*" = as_qspray_gmp(e1) * e2,
    stop(gettextf(
      "Binary operator %s not defined for qspray objects.", dQuote(.Generic)
    ))
  )
}

numeric_arith_qspray <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as.qspray.numeric(e1) + e2,
    "-" = as.qspray.numeric(e1) - e2,
    "*" = as.qspray.numeric(e1) * e2,
    stop(gettextf(
      "Binary operator %s not defined for qspray objects.", dQuote(.Generic)
    ))
  )
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

#' Title
#'
#' @param P xx
#' @param S xx
#'
#' @return xx
#' @export
integratePolynomialOnSimplex <- function(P, S) {
  S <- as.bigq(S)
  n <- ncol(S)
  v <- t(S[n+1L, ])
  B <- t(S[1L:n, ]) - do.call(function(...) cbind(...), replicate(n, v))
  gens <- lapply(1L:n, function(i) lone(i))
  newvars <- vector("list", n)
  for(i in 1L:n) {
    newvar <- v[i]
    Bi <- B[i, ]
    for(j in 1L:n) {
      newvar <- newvar + Bi[j] * gens[[j]]
    }
    newvars[[i]] <- newvar
  }
  Q <- as.bigq(0L)
  exponents <- P@powers
  coeffs    <- P@coeffs 
  for(i in 1L:length(exponents)) {
    powers <- exponents[[i]]
    term <- as.bigq(1L)
    for(j in 1L:length(powers)) {
      term <- term * newvars[[j]]^powers[j] 
    }
    Q <- Q + coeffs[i] * term
  }
  s <- as.bigq(0L)
  exponents <- Q@powers
  coeffs    <- Q@coeffs 
  for(i in 1L:length(exponents)) {
    coef <- as.bigq(coeffs[i])
    powers <- exponents[[i]]
    d <- sum(powers)
    if(d == 0L) {
      s <- s + coef
      next
    }
    coef <- coef * prod(factorialZ(powers))
    s <- s + coef / prod((n+1L):(n+d))
  }
  abs(as.bigq(detQ(as.character(B)))) *  s / factorialZ(n)
}
