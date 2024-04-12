# orderedQspray <- function(qspray, d) {
#   powers <- qspray@powers
#   Mpowers <- do.call(rbind, lapply(powers, grow, n = d))
#   ordr <- lexorder(Mpowers)
#   list(
#     "powers" = Mpowers[ordr, , drop = FALSE], 
#     "coeffs" = qspray@coeffs[ordr]
#   )
# }

lexLeading <- function(M, i = 1L, b = seq_len(nrow(M))) {
  if(nrow(M) == 1L || i > ncol(M)) {
    b[1L]
  } else {
    col_i <- M[, i]
    mx <- max(col_i) == col_i
    lexLeading(M[mx, , drop = FALSE], i + 1L, b[mx])
  }
}

leading <- function(qspray, d) {
  powers <- qspray@powers
  Mpowers <- do.call(rbind, lapply(powers, grow, n = d))
  i <- lexLeadingArma(Mpowers)
  list("powers" = Mpowers[i, ], "coeff" = qspray@coeffs[i])
}

leadingTerm <- function(qspray, d) {
  l <- leading(qspray, d)
  list("powers" = l[["powers"]], "coeff" = as.bigq(l[["coeff"]]))
}

coeffInverse <- function(coeff) {
  if(grepl("^-", coeff)) {
    coeff <- sub("^-", "", coeff)
    paste0("-", paste0(strsplit(coeff, "/")[[1L]][c(2L, 1L)], collapse = "/"))
  } else {
    paste0(strsplit(coeff, "/")[[1L]][c(2L, 1L)], collapse = "/")
  }
}

# S polynomial ####
S <- function(f, g) {
  d <- max(arity(f), arity(g))
  leading_f <- leading(f, d)
  leading_g <- leading(g, d)
  lpows_f <- leading_f[["powers"]]
  lpows_g <- leading_g[["powers"]]
  lcoef_f <- leading_f[["coeff"]]
  lcoef_g <- leading_g[["coeff"]]
  gamma <- pmax(lpows_f, lpows_g)
  beta_f <- gamma - lpows_f
  beta_g <- gamma - lpows_g
  w_f <- qsprayMaker(list(beta_f), coeffInverse(lcoef_f))
  w_g <- qsprayMaker(list(beta_g), coeffInverse(lcoef_g))
  w_f * f - w_g * g
}


# division ####

divides <- function(g, f) { # whether term g divides term f
  all(g[["powers"]] <= f[["powers"]])
}

quotient <- function(f, g) { # quotient of term f divided by term g, assuming 
  powers <- f[["powers"]] - g[["powers"]]                       # g divides f
  coeff  <- f[["coeff"]] / g[["coeff"]]
  qsprayMaker(powers = list(powers), coeffs = coeff)
}

termAsQspray <- function(term) {
  qsprayMaker(powers = list(term[["powers"]]), coeffs = term[["coeff"]])
}

#' @title Division of a qspray polynomial
#' @description Division of a qspray polynomial by a list of qspray 
#'   polynomials. See the reference for the definition.
#' 
#' @param qspray the dividend, a \code{qspray} object 
#' @param divisors the divisors, a list of \code{qspray} objects
#'
#' @return The remainder of the division, a \code{qspray} object. 
#' @export
#'
#' @references 
#' Michael Weiss, 2010. 
#' \href{https://math.nyu.edu/degree/undergrad/ug_research/Weiss_SURE_Paper.pdf}{Computing Gröbner Bases in Python with Buchberger’s Algorithm}.
#' 
#' @examples
#' # a univariate example
#' library(qspray)
#' x <- qlone(1)
#' f <- x^4 - 4*x^3 + 4*x^2 - x # 0 and 1 are trivial roots
#' g <- x * (x - 1)
#' qdivision(f, list(g)) # should be zero
qdivision <- function(qspray, divisors) {
  stopifnot(is.list(divisors))
  if(qspray == qzero()) {
    return(qzero())
  }
  d <- max(vapply(divisors, arity, integer(1L)))
  LTdivisors <- lapply(divisors, leading, d = d)
  BBdivision(qspray, divisors, LTdivisors)
  
  # # we store the successive leading terms in LTs_f
  # d <- max(arity(qspray), max(vapply(divisors, arity, integer(1L))))
  # oqspray <- orderedQspray(qspray, d)
  # opowers <- oqspray[["powers"]]
  # ocoeffs <- as.bigq(oqspray[["coeffs"]])
  # LTs_f <- lapply(seq_along(ocoeffs), function(i) {
  #   list("powers" = opowers[i, ], "coeff" = ocoeffs[i])
  # })
  # 
  # ndivisors <- length(divisors)
  # nterms <- length(qspray@coeffs)
  # 
  # qgs <- list() # to store the products q*g_i, in order to check at the end
  # quotients <- list() # to store the quotients
  # 
  # cur <- qspray
  # for(k in 1L:nterms) {
  #   LT_cur <- LTs_f[[k]]
  #   i <- 1L
  #   while(i <= ndivisors) {
  #     g <- divisors[[i]]
  #     LT_g <- leadingTerm(g, d)
  #     while(divides(LT_g, LT_cur)) {
  #       q <- quotient(LT_cur, LT_g)
  #       quotients <- append(quotients, q)
  #       qgs <- append(qgs, q * g)
  #       cur <- cur - q * g
  #       if(cur == qzero()) {
  #         if(check) {
  #           sum_qgs <- qzero()
  #           for(i in seq_along(qgs)) {
  #             sum_qgs <- sum_qgs + qgs[[i]]
  #           }
  #           stopifnot(sum_qgs == qspray)
  #         }
  #         remainder <- qzero()
  #         if(d == 1L) {
  #           qtnt <- qzero()
  #           for(i in seq_along(quotients)) {
  #             qtnt <- qtnt + quotients[[i]]
  #           }
  #           attr(remainder, "quotient") <- qtnt
  #         }
  #         return(remainder)
  #       }
  #       LT_cur <- leadingTerm(cur, d)
  #     }
  #     i <- i + 1L
  #   }
  # }
  # # check
  # if(check) {
  #   sum_qgs <- qzero()
  #   for(i in seq_along(qgs)) {
  #     sum_qgs <- sum_qgs + qgs[[i]]
  #   }
  #   stopifnot(qspray == sum_qgs + cur)
  # }
  # # return remainder
  # remainder <- cur
  # if(d == 1L) {
  #   qtnt <- qzero()
  #   for(i in seq_along(quotients)) {
  #     qtnt <- qtnt + quotients[[i]]
  #   }
  #   attr(remainder, "quotient") <- qtnt
  # }
  # remainder
}

# internal division for Buchberger algorithm
BBdivision <- function(qspray, divisors, LTdivisors) {
  if(qspray == qzero()) {
    return(qzero())
  }
  # we store the successive leading terms in LTs_f
  d <- max(arity(qspray), max(vapply(divisors, arity, integer(1L))))
  # oqspray <- orderedQspray(qspray, d)
  # opowers <- oqspray[["powers"]]
  # ocoeffs <- oqspray[["coeffs"]]
  # LTs_f <- lapply(seq_along(ocoeffs), function(i) {
  #   list("powers" = opowers[i, ], "coeff" = ocoeffs[i])
  # })
  gs <- lapply(divisors, function(qspr) {
    list("powers" = qspr@powers, "coeffs" = qspr@coeffs)
  })
  Powers <- qspray@powers
  coeffs <- qspray@coeffs
  outList <- BBdivisionRcpp(
    Powers, coeffs, gs, LTdivisors, d
  )
  qspray_from_list(outList)  
}

combn2 <- function(j, s) {
  allCombs <- rbind(
    do.call(c, lapply(1L:(j-1L), function(i) 1L:i)),
    rep(2L:j, times = 1L:(j-1L))  
  )
  allCombs[, (s+1L):ncol(allCombs), drop = FALSE]
}

#' @title Gröbner basis
#' @description Returns a Gröbner basis following Buchberger's algorithm 
#'   using the lexicographical order.
#' @param G a list of qspray polynomials, the generators of the ideal
#' @param minimal Boolean, whether to return a minimal basis
#' @param reduced Boolean, whether to return the reduced basis
#' @return A Gröbner basis of the ideal generated by \code{G}, given as a list 
#'   of qspray polynomials.
#' @export
#' @references 
#' Cox, Little & O'Shea. 
#' \emph{Ideals, Varieties, and Algorithms. 
#' An Introduction to Computational Algebraic Geometry and Commutative Algebra.}
#' Fourth edition, Springer 2015.
#' @examples
#' library(qspray)
#' f <- qsprayMaker(string = "x^(3) - 2 x^(1,1)")
#' g <- qsprayMaker(string = "x^(2,1) - 2 x^(0,2) + x^(1)")
#' groebner(list(f, g), FALSE, FALSE)
#' # other example
#' \donttest{x <- qlone(1); y <- qlone(2); z <- qlone(3)
#' f1 <- x^2 + y + z^2 - 1
#' f2 <- x^2 + y + z - 1
#' f3 <- x + y^2 + z - 1
#' gb <- groebner(list(f1, f2, f3))
#' lapply(gb, prettyQspray, vars = c("x", "y", "z"))}
groebner <- function(G, minimal = TRUE, reduced = TRUE) {
  # d <- max(vapply(G, arity, integer(1L)))
  d <- max(vapply(G, numberOfVariables, integer(1L)))
  LT_G <- lapply(G, leading, d = d)
  Ss <- list()
  j <- length(G)
  combins <- combn2(j, 0L)
  i <- 1L
  l <- ncol(combins)
  while(i <= l) {
    combin <- combins[, i]
    Sfg <- S(G[[combin[1L]]], G[[combin[2L]]])
    d <- max(d, arity(Sfg))
    Sbar_fg <- BBdivision(Sfg, G, LT_G)
    if(Sbar_fg != qzero()) {
      G <- append(G, Sbar_fg)
	    d <- max(d, arity(Sbar_fg))
      LT_G <- append(LT_G, list(leading(Sbar_fg, d)))
      j <- j + 1L
      combins <- combn2(j, i)
      l <- ncol(combins)
      i <- 1L
    } else {
      i <- i + 1L
    }
  }
  #
  if(minimal || reduced) {
    d <- max(vapply(G, arity, integer(1L)))
    indices <- seq_along(G)
    toRemove <- drop <- integer(0L)
    for(i in indices) {
      LT_f <- leadingTerm(G[[i]], d)
      drop <- c(toRemove, i)
      for(j in setdiff(indices, drop)) {
        if(divides(leadingTerm(G[[j]], d), LT_f)) {
          toRemove <- c(toRemove, i)
          break
        }
      }
    }
    # normalization
    if(length(toRemove) > 0L) {
      G <- G[-toRemove]
      for(i in seq_along(G)) {
        G[[i]] <- G[[i]] / leadingTerm(G[[i]], d)[["coeff"]]
      }
    }
    # reduction
    if(reduced && length(G) > 1L) {
      indices <- seq_along(G)
      for(i in indices) {
        keep <- setdiff(indices, i)
        G[[i]] <- qdivision(G[[i]], G[keep])
      }
    }
  }
  G
}

#' @title Implicitization with Gröbner bases
#' @description Implicitization of a system of parametric equations 
#'   (see examples).
#'
#' @param nvariables number of variables
#' @param parameters character vector of the names of the parameters, or 
#'   \code{NULL} if there's no parameter
#' @param equations list of qspray polynomials representing the parametric 
#'   equations
#' @param relations list of qspray polynomials representing the relations 
#'   between the variables and the parameters, or \code{NULL} if there is none
#'
#' @return A list of qspray polynomials.
#' @export
#' @importFrom utils tail
#'
#' @examples
#' library(qspray)
#' # ellipse example ####
#' # variables 
#' cost <- qlone(1)
#' sint <- qlone(2)
#' # parameters
#' a <- qlone(3)
#' b <- qlone(4)
#' #
#' nvariables <- 2
#' parameters <- c("a", "b")
#' equations <- list(
#'   "x" = a * cost,
#'   "y" = b * sint
#' )
#' relations <- list(
#'   cost^2 + sint^2 - 1
#' )
#' # 
#' eqs <- implicitization(nvariables, parameters, equations, relations)
implicitization <- function(nvariables, parameters, equations, relations) {
  stopifnot(isPositiveInteger(nvariables))
  stopifnot(is.null(parameters) || isStringVector(parameters))
  stopifnot(is.list(equations), length(equations) > 1L)
  stopifnot(is.null(relations) || is.list(relations))
  #
  nequations <- length(equations)
  nrelations <- length(relations)
  nqlone <- max(vapply(equations, arity, integer(1L)))
  coordinates <- lapply((nqlone+1L):(nqlone+nequations), qlone)
  generators <- relations
  for(i in seq_along(equations)) {
    generators <- append(generators, coordinates[[i]] - equations[[i]])
  }
  #
  gb <- groebner(generators)
  isfree <- function(i) {
    all(vapply(gb[[i]]@powers, function(pows) {
      length(pows) == 0L || 
        (length(pows > nvariables) && all(pows[1L:nvariables] == 0L))
    }, logical(1L)))
  }
  free <- c(FALSE, vapply(2L:length(gb), isfree, logical(1L)))
  #
  results <- gb[free]
  for(i in seq_along(results)) {
    el <- results[[i]]
    coeffs <- el@coeffs
    powers <- el@powers
    for(j in seq_along(powers)) {
      powers[[j]] <- tail(powers[[j]], -nvariables)
    }
    results[[i]] <- qsprayMaker(powers, coeffs)
  }
  vars <- c(parameters, names(equations))
  messages <- lapply(results, prettyQspray, vars = vars)
  for(msg in messages) {
    message(msg)
  }
  invisible(results)
}

#' @title Whether a 'qspray' is a polynomial of some given 'qsprays'
#' @description Checks whether a \code{qspray} polynomial can be written as 
#'   a polynomial of some given \code{qspray} polynomials. If \code{TRUE}, 
#'   this polynomial is returned.
#' 
#' @param qspray a \code{qspray} object
#' @param qsprays a list of \code{qspray} objects
#'
#' @return A Boolean value indicating whether the polynomial defined by 
#'   \code{qspray} can be written as a polynomial of the polynomials defined 
#'   by the \code{qspray} objects given in the \code{qsprays} list. If this is 
#'   \code{TRUE}, this polynomial is returned as an attribute named 
#'   \code{"polynomial"}.
#' @export
#'
#' @examples
#' library(qspray)
#' P <- function(X, Y) X^2*Y + 2*X + 3
#' x <- qlone(1); y <- qlone(2); z <- qlone(3)
#' q1 <- x + y
#' q2 <- x*z^2 + 4
#' qspray <- P(q1, q2)
#' ( check <- isPolynomialOf(qspray, list(q1, q2)) )
#' POLYNOMIAL <- attr(check, "polynomial")
#' composeQspray(POLYNOMIAL, list(q1, q2)) == qspray # should be TRUE
isPolynomialOf <- function(qspray, qsprays) {
  n <- max(vapply(qsprays, numberOfVariables, integer(1L)))
  if(numberOfVariables(qspray) > n) {
    return(FALSE)
  }
  i_ <- seq_len(length(qsprays))
  G <- lapply(i_, function(i) qsprays[[i]] - qlone(n + i))
  B <- groebner(G, TRUE, FALSE)
  constantTerm <- getCoefficient(qspray, integer(0L))
  g <- qdivision(qspray - constantTerm, B)
  check <- all(vapply(g@powers, function(pwr) {
    length(pwr) > n && all(pwr[1L:n] == 0L)
  }, logical(1L)))
  if(!check) {
    return(FALSE)
  }
  powers <- lapply(g@powers, function(pwr) {
    pwr[-(1L:n)]
  })
  P <- qsprayMaker(powers, g@coeffs) + constantTerm
  out <- TRUE
  attr(out, "polynomial") <- P
  out
}
