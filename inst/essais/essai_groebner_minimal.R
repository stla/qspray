library(qspray)

f <- qsprayMaker(string = "x^(3) - 2 x^(1,1)")
g <- qsprayMaker(string = "x^(2,1) - 2 x^(0,2) + x^(1)")
f <- qsprayMaker(string = "-2 x^(1,1) + x^(1)")
g <- qsprayMaker(string = "x^(3,1) - 2 x^(2) + x^(0,1)")
G <- groebner(list(f, g))

f1 <- qsprayMaker(string = "x^(2) + x^(0,1) + x^(0,0,2) - 1")
f2 <- qsprayMaker(string = "x^(2) + x^(0,1) + x^(0,0,1) - 1")
f3 <- qsprayMaker(string = "x^(1) + x^(0,2) + x^(0,0,1) - 1")
G <- groebner(list(f1, f2, f3))


# R.<x ,y, z> = PolynomialRing(QQ, 3, order = "lex")
# I = R.ideal([x^2 + y + z^2 - 1, x^2 + y + z - 1, x + y^2 + z - 1])
# I.groebner_basis()
# # minimale
# P_1 = z + y^2 + x - 1
# P_2 = z^2 - z
# P_3 = -z + 2y^2z + y^4 + y + z^2 - 2y^2
# 
# # réduite:
# P1, P2, -2y^2 + y + y^4 + 2*y^2z

d <- 3L
l <- length(G)
toRemove <- integer(0L) -> drop
for(i in seq_len(l)) {
  LT_f <- qspray:::leadingTerm(G[[i]], d)
  drop <- c(toRemove, i)
  print(drop)
  removal <- FALSE
  for(j in setdiff(seq_len(l), drop)) {
    if(qspray:::divides(qspray:::leadingTerm(G[[j]], d), LT_f)) {
      toRemove <- c(toRemove, i)
      removal <- TRUE
      break
    }
  }
}
if(length(toRemove) > 0L) {
  G <- G[-toRemove]
  Gnorm <- vector("list", length(G))
  for(i in seq_along(Gnorm)) {
    Gnorm[[i]] <- G[[i]] / qspray:::leadingTerm(G[[i]], d)$coeff
  }
}
Gnorm

# reduced base - useless if #Gnorm=1
G <- Gnorm
indices <- seq_along(G)
bar <- function(f, i) {
  keep <- setdiff(indices, i)
  qspray:::divisionRemainder(f, G[keep])  
}

for(i in indices) {
  G[[i]] <- bar(G[[i]], i)
}
G