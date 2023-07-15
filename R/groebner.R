grow <- function(powers, n) {
  c(powers, integer(n - length(powers)))
}

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
  i <- lexLeading(Mpowers)
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
#' @param check Boolean, whether to check the division; this argument will be 
#'   removed in a future version
#'
#' @return The remainder of the division, a \code{qspray} object. Moreover, 
#'   if \code{qspray} is univariate, the quotient is attached to the remainder 
#'   as an attribute.
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
#' ( r <- qdivision(f, list(g)) ) # should be zero
#' attr(r, "quotient")
qdivision <- function(qspray, divisors, check = FALSE) {
  
  stopifnot(is.list(divisors))
  
  if(qspray == qzero()) {
    return(qzero())
  }
  
  d <- max(arity(qspray), max(vapply(divisors, arity, integer(1L))))

  ndivisors <- length(divisors)
  nterms <- length(qspray@coeffs)
  LTs_f <- vector("list", nterms) # to store the successive leading terms
  qgs <- list() # to store the products q*g_i, in order to check at the end
  quotients <- list()
  
  cur <- qspray
  for(k in 1L:nterms) {
    # take the next leading term of f
    tmp <- qspray
    for(j in seq_len(k-1L)) {
      tmp <- tmp - termAsQspray(LTs_f[[j]])
    }
    
    # cat("LTs_f[[k]] - k=", k, "\n")
    
    LTs_f[[k]] <- leadingTerm(tmp, d) -> LT_cur
    i <- 1L
    while(i <= ndivisors) {
      g <- divisors[[i]]
      
      # cat("LT_g - i=", i, "\n")
      
      LT_g <- leadingTerm(g, d)
      while(divides(LT_g, LT_cur)) {
        
        # print("quotient(LT_cur, LT_g)")
        
        q <- quotient(LT_cur, LT_g)
        # quotients <- append(quotients, q)
        # qgs <- append(qgs, q * g)
        cur <- cur - q * g
        if(cur == qzero()) {
          if(check) {
            sum_qgs <- qzero()
            for(i in seq_along(qgs)) {
              sum_qgs <- sum_qgs + qgs[[i]]
            }
            stopifnot(sum_qgs == qspray)
          }
          remainder <- qzero()
          # if(d == 1L) {
          #   qtnt <- qzero()
          #   for(i in seq_along(quotients)) {
          #     qtnt <- qtnt + quotients[[i]]
          #   }
          #   attr(remainder, "quotient") <- qtnt
          # }
          return(remainder)
        }
        # print("LT_cur")
        LT_cur <- leadingTerm(cur, d)
      }
      i <- i + 1L
    }
  }
  # check
  if(check) {
    sum_qgs <- qzero()
    for(i in seq_along(qgs)) {
      sum_qgs <- sum_qgs + qgs[[i]]
    }
    stopifnot(qspray == sum_qgs + cur)
  }
  # return remainder
  remainder <- cur
  # if(d == 1L) {
  #   qtnt <- qzero()
  #   for(i in seq_along(quotients)) {
  #     qtnt <- qtnt + quotients[[i]]
  #   }
  #   attr(remainder, "quotient") <- qtnt
  # }
  remainder
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
#' @importFrom utils combn
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
  Ss <- list()
  j <- length(G)
  combins <- combn(j, 2L)
  i <- 1L
  indices <- 1L:ncol(combins)
  while(i <= length(indices)) {
    combin <- combins[, indices[i]]
    id <- paste0(combin[1L], "-", combin[2L])
    if(id %in% names(Ss)) {
      print("that should not happen")
      #Sfg <- Ss[[id]]
      Sbar_fg <- qzero()
    } else {
	#print("calc Sfg")
      Sfg <- S(G[[combin[1L]]], G[[combin[2L]]])
      Ss_new <- list(Sfg)
      names(Ss_new) <- id
      Ss <- c(Ss, Ss_new)
	    # print("calc division")
      Sbar_fg <- qdivision(Sfg, G)
    }
    i <- i + 1L
    if(Sbar_fg != qzero()) {
      i <- 1L
      G <- append(G, Sbar_fg)
      j <- j + 1L
      #print(j)
      combins <- combn(j, 2L)
      allids <- paste0(combins[1L, ], "-", combins[2L, ])
      indices <- which(!is.element(allids, names(Ss)))
	#cat("nindices: ", length(indices), "\n") 
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
    #
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
