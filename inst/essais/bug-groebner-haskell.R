library(qspray)

xorderedQspray <- function(qspray, d) {
  powers <- qspray@powers
  Mpowers <- do.call(rbind, lapply(powers, qspray:::grow, n = d))
  ordr <- qspray:::lexorder(Mpowers)
  list(
    "powers" = Mpowers[ordr, , drop = FALSE],
    "coeffs" = qspray@coeffs[ordr]
  )
}

DIVISION <- function(qspray, divisors, check = FALSE) {
  # we store the successive leading terms in LTs_f
  d <- max(numberOfVariables(qspray), max(vapply(divisors, numberOfVariables, integer(1L))))
  oqspray <- xorderedQspray(qspray, d)
  opowers <- oqspray[["powers"]]
  ocoeffs <- gmp::as.bigq(oqspray[["coeffs"]])
  LTs_f <- lapply(seq_along(ocoeffs), function(i) {
    list("powers" = opowers[i, ], "coeff" = ocoeffs[i])
  })

  ndivisors <- length(divisors)
  nterms <- length(qspray@coeffs)

  qgs <- list() # to store the products q*g_i, in order to check at the end
  quotients <- list() # to store the quotients

  cur <- qspray
  for(k in 1L:nterms) {
    LT_cur <- LTs_f[[k]]
    i <- 1L
    while(i <= ndivisors) {
      g <- divisors[[i]]
      LT_g <- qspray:::leadingTerm(g, d)
      while(qspray:::divides(LT_g, LT_cur)) {
        q <- qspray:::quotient(LT_cur, LT_g)
        quotients <- append(quotients, q)
        qgs <- append(qgs, q * g)
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
          if(d == 1L) {
            qtnt <- qzero()
            for(i in seq_along(quotients)) {
              qtnt <- qtnt + quotients[[i]]
            }
            attr(remainder, "quotient") <- qtnt
          }
          return(remainder)
        }
        LT_cur <- qspray:::leadingTerm(cur, d)
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
  if(d == 1L) {
    qtnt <- qzero()
    for(i in seq_along(quotients)) {
      qtnt <- qtnt + quotients[[i]]
    }
    attr(remainder, "quotient") <- qtnt
  }
  remainder
}


groebnerTest <- function(G, minimal = TRUE, reduced = TRUE) {
  d <- max(vapply(G, numberOfVariables, integer(1L)))
  LT_G <- lapply(G, qspray:::leading, d = d)
  Ss <- list()
  j <- length(G)
  combins <- qspray:::combn2(j, 0L)
  i <- 1L
  l <- ncol(combins)
  while(j <= 100 && i <= l) {
    combin <- combins[, i]
    Sfg <- qspray:::S(G[[combin[1L]]], G[[combin[2L]]])
    d <- max(d, numberOfVariables(Sfg))
    #Sbar_fg <- DIVISION(Sfg, G)
    Sbar_fg <- qspray:::BBdivision(Sfg, G, LT_G)
    if(!isQzero(Sbar_fg)) {
      G <- append(G, Sbar_fg)
      # if(length(G) == 8L) {
      #   return(list("G" = G, "Sfg" = Sfg))
      # }
      d <- max(d, numberOfVariables(Sbar_fg))
      LT_G <- append(LT_G, list(qspray:::leading(Sbar_fg, d)))
      j <- j + 1L
      combins <- qspray:::combn2(j, i)
      l <- ncol(combins)
      i <- 1L
    } else {
      i <- i + 1L
    }
  }
  G
}

x1 <- qlone(1); x2 <- qlone(2); x3 <- qlone(3); x4 <- qlone(4)
x5 <- qlone(5); x6 <- qlone(6); x7 <- qlone(7); x8 <- qlone(8)

G <- list(
  ESFpoly(4, 1) - x5,
  x1*x2 + x1*x3 + x1*x4 + x2*x3 - x6,
  ESFpoly(4, 3) - x7,
  ESFpoly(4, 4) - x8
)

gbasis <- groebnerTest(G, F, F)
length(gbasis)

reduceGB <- function(G) {
  d <- max(vapply(G, numberOfVariables, integer(1L)))
  indices <- seq_along(G)
  toRemove <- drop <- integer(0L)
  for(i in indices) {
    LT_f <- qspray:::leadingTerm(G[[i]], d)
    drop <- c(toRemove, i)
    for(j in setdiff(indices, drop)) {
      if(qspray:::divides(qspray:::leadingTerm(G[[j]], d), LT_f)) {
        toRemove <- c(toRemove, i)
        break
      }
    }
  }
  # normalization
  if(length(toRemove) > 0L) {
    G <- G[-toRemove]
    for(i in seq_along(G)) {
      G[[i]] <- G[[i]] / qspray:::leadingTerm(G[[i]], d)[["coeff"]]
    }
  }
  G
}

rgbasis <- reduceGB(gbasis)
length(rgbasis)


## 
library(giacR)
giac <- Giac$new()
command <- 
  "gbasis([x1 + x2 + x3 + x4 - x5, x1*x2 + x1*x3 + x1*x4 + x2*x3 - x6, x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - x7, x1*x2*x3*x4 - x8],[x1,x2,x3,x4,x5,x6,x7,x8])"
giac$execute(command)

