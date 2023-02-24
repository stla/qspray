#' @useDynLib qspray, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods setMethod setClass new
#' @importFrom gmp as.bigq factorialZ asNumeric
#' @importFrom purrr transpose
#' @include qspray.R
NULL

setClass(
  "qspray",
  slots = c(powers = "list", coeffs = "character")
)

showQspray <- function(qspray) {
  if(length(qspray@coeffs) == 0L) {
    return("0")
  }
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

#' @title Make a 'qspray' object
#' @description Make a \code{qspray} object from a list of exponents and a 
#'   vector of coefficients.
#'
#' @param powers list of positive integer vectors
#' @param coeffs a vector such that each element of \code{as.character(coeffs)} 
#'   is a quoted integer or a quoted fraction; it must have the same length 
#'   as the \code{powers} list
#' @param string if not \code{NULL}, this argument takes precedence over 
#'   \code{powers} and \code{vertices}; it must be a string representing a 
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

#' @title Evaluate a 'qspray' object
#' @description Evaluation of the multivariate polynomial represented by a 
#'   \code{qspray} object.
#'
#' @param qspray a \code{qspray} object
#' @param values_re vector of the real parts of the values; each element of 
#'   \code{as.character(values_re)} must be quoted integer or a quoted fraction
#' @param values_im vector of the imaginary parts of the values; each element of 
#'   \code{as.character(values_im)} must be quoted integer or a quoted fraction
#'
#' @return A \code{bigq} number if \code{values_im=NULL}, a pair of \code{bigq} 
#'   numbers otherwise: the real part and the imaginary part of the result.
#' @export
#' @examples 
#' x <- qlone(1); y <- qlone(2)
#' P <- 2*x + "1/2"*y
#' evalQspray(P, c("2", "5/2", "99999")) # "99999" will be ignored
evalQspray <- function(qspray, values_re, values_im = NULL) {
  powers <- qspray@powers
  if(length(powers) == 0L) {
    return(as.bigq(0L))
  }
  n <- arity(qspray)
  if(length(values_re) < n) {
    stop("Insufficient number of values.")
  }
  values_re <- as.character(values_re)
  check <- all(vapply(values_re, isFraction, logical(1L)))
  if(!check) {
    stop("Invalid vector `values_re`.")
  }
  if(!is.null(values_im)) {
    stopifnot(length(values_re) == length(values_im))
    values_im <- as.character(values_im)
    check <- all(vapply(values_im, isFraction, logical(1L)))
    if(!check) {
      stop("Invalid vector `values_im`.")
    }
    result <- evalQxspray(powers, qspray@coeffs, values_re, values_im)
    return(as.bigq(result))
  }
  coeffs <- as.bigq(qspray@coeffs)
  values <- as.bigq(values_re)
  out <- 0
  for(i in seq_along(powers)) {
    exponents <- powers[[i]]
    term <- 1
    for(j in seq_along(exponents)) {
      term <- term * values[j]^exponents[j]
    }
    out <- out + coeffs[i] * term
  }
  out
}

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
      "Binary operator %s not defined for qspray objects.", dQuote(.Generic)
    ))
  )
}

qsprayPower <- function(qspray, n) {
  stopifnot(isPositiveInteger(n))
  qspray_power(qspray@powers, qspray@coeffs, n)
}

#' @title Partial derivative
#' @description Partial derivative of a qspray polynomial.
#'
#' @param qspray object of class \code{qspray}
#' @param i integer, the dimension to differentiate with respect to
#' @param derivative integer, how many times to differentiate
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' library(qspray)
#' x <- qlone(1)
#' y <- qlone(2)
#' qspray <- 2*x  + 3*x*y
#' derivQspray(qspray, 1)
derivQspray <- function(qspray, i, derivative = 1) {
  stopifnot(inherits(qspray, "qspray"))
  stopifnot(isNonnegativeInteger(i))
  stopifnot(isPositiveInteger(derivative))
  if(i > arity(qspray)) {
    return(as.qspray(0))
  }
  n    <- integer(length = i)
  n[i] <- as.integer(derivative)
  drv  <- qspray_deriv(qspray@powers, qspray@coeffs, n)
  qspray_from_list(drv)
}

#' @title Partial differentiation
#' @description Partial differentiation of a qspray polynomial.
#'
#' @param qspray object of class \code{qspray}
#' @param orders integer vector, the orders of the differentiation
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' library(qspray)
#' x <- qlone(1)
#' y <- qlone(2)
#' qspray <- x + 2*y  + 3*x*y
#' dQspray(qspray, c(1, 1))
#' derivQspray(derivQspray(qspray, 1), 2)
dQspray <- function(qspray, orders) {
  stopifnot(inherits(qspray, "qspray"))
  for(i in seq_along(orders)) {
    stopifnot(isPositiveInteger(orders[i]))
  }
  orders <- removeTrailingZeros(orders)
  if(length(orders) > arity(qspray)) {
    return(as.qspray(0))
  }
  n    <- as.integer(orders)
  drv  <- qspray_deriv(qspray@powers, qspray@coeffs, n)
  qspray_from_list(drv)
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
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
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
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
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
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
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
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
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
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
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

#' @title Integral of a multivariate polynomial over a simplex
#' @description Returns the exact value of the integral of a multivariate 
#'   polynomial with rational coefficients over a simplex whose vertices have 
#'   rational coordinates.
#'
#' @param P a \code{qspray} object
#' @param S the simplex, a \code{(n+1)xn} matrix such that each entry of the  
#'   matrix \code{as.character(S)} is a quoted integer or a quoted fraction
#'
#' @return A \code{bigq} number, the exact value of the integral.
#' @export
#' @importFrom RationalMatrix Qdet
#' @examples 
#' library(qspray)
#' x <- qlone(1); y <- qlone(2)
#' P <- x/2 + x*y
#' S <- rbind(c("0", "0"), c("1", "0"), c("1", "1")) # a triangle
#' integratePolynomialOnSimplex(P, S)
integratePolynomialOnSimplex <- function(P, S) {
  storage.mode(S) <- "character"
  check <- all(vapply(S, isFraction, logical(1L)))
  if(!check) {
    stop("Invalid entries in the matrix `S`.")
  }
  n <- ncol(S)
  if(nrow(S) != n+1L) {
    stop("The matrix `S` does not represent a simplex.")
  }
  S <- as.bigq(S)
  v <- t(S[n+1L, ])
  B <- t(S[1L:n, ]) - do.call(function(...) cbind(...), replicate(n, v))
  gens <- lapply(1L:n, function(i) qlone(i))
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
    for(j in seq_along(powers)) {
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
  abs(as.bigq(Qdet(as.character(B)))) *  s / factorialZ(n)
}

#' @title Monomial symmetric function
#' @description Returns a monomial symmetric function as a polynomial.
#'
#' @param m integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#'   positive integers
#'
#' @return A \code{qspray} object.
#' @importFrom DescTools Permn
#' @export
#'
#' @examples
#' library(qspray)
#' MSFpoly(3, c(3, 1))
MSFpoly <- function(m, lambda) {
  stopifnot(isNonnegativeInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0]
  if(length(lambda) > m) return(as.qspray(0))
  kappa <- numeric(m)
  kappa[seq_along(lambda)] <- lambda
  perms <- Permn(kappa)
  n <- nrow(perms)
  powers <- lapply(1L:n, function(i) perms[i, ])
  coeffs <- rep("1", n)
  qsprayMaker(powers, coeffs)
}

#' @title Elementary symmetric function
#' @description Returns an elementary symmetric function as a polynomial.
#'
#' @param m integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#'   positive integers
#'
#' @return A \code{qspray} object.
#' @importFrom DescTools Permn
#' @export
#'
#' @examples
#' library(qspray)
#' ESFpoly(3, c(3, 1))
ESFpoly <- function(m, lambda) {
  stopifnot(isNonnegativeInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0]
  if(any(lambda > m)) return(as.qspray(0))
  out <- 1
  for(k in seq_along(lambda)) {
    kappa <- integer(m)
    kappa[seq_len(lambda[k])] <- rep(1L, lambda[k])
    perms <- Permn(kappa)
    n <- nrow(perms)
    powers <- lapply(1L:n, function(i) perms[i, ])
    ek <- qsprayMaker(powers = powers, coeffs = rep("1", n))
    out <- out * ek
  }
  out
}
