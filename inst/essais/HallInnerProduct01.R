library(qspray)
library(gmp)

spray1 <- jack::SchurPol(3, c(3, 1))
spray2 <- jack::SchurPol(4, c(3, 1))
HallInnerProduct(spray1, spray2)

alpha <- as.bigq(2L)
t <- as.integer(alpha)
spray1 <- jack::JackPol(3, c(2, 1), alpha)
spray2 <- jack::JackPol(3, c(2, 1), alpha)
HallInnerProduct(spray1, spray2, alpha = t)
spray1 <- jack::ZonalPol(3, c(2,1))
spray2 <- jack::ZonalPol(3, c(2,1))
HallInnerProduct(spray1, spray2, alpha = t)
t <- alpha
3*t^3/(t^2 + 3/2*t + 1/2)
(2*t^3 + t^2)/(t + 2) 
t^3/6 + t^2/2 + t/3

poly <- PSFpoly(4, c(4))
HallInnerProduct(poly, poly, alpha = t)
4*t
poly <- PSFpoly(4, c(3,1))
HallInnerProduct(poly, poly, alpha = t)
3*t^2
poly <- PSFpoly(4, c(2,2))
HallInnerProduct(poly, poly, alpha = t)
8*t^2
poly <- PSFpoly(4, c(2,1,1))
HallInnerProduct(poly, poly, alpha = t)
4*t^3
poly <- PSFpoly(4, c(1,1,1,1))
HallInnerProduct(poly, poly, alpha = t)
24 * t^4

HallInnerProduct(
  jack::JackPol(3, c(2,1), alpha), PSFpoly(3, c(2,1)), alpha = t
)
t + 2 # no

# 
library(jack)

skewJackPol <- function(n, lambda, mu, alpha) {
  Jlambda <- JackPolCPP(n, lambda, alpha)
  Jmu     <- JackPolCPP(n, mu, alpha)
  nus <- partitions::parts(sum(lambda) - sum(mu))
  terms <- apply(nus, 2L, function(nu) {
    Jnu <- JackPolCPP(n, nu, alpha)
    coeff <- HallInnerProduct(Jlambda, Jmu * Jnu, alpha) / # rq: Jlambda décomposé plusieurs fois !
      HallInnerProduct(Jnu, Jnu, alpha) 
    coeff * Jnu
  }, simplify = FALSE)
  Reduce(`+`, terms)
}

skewJackPol(3, c(2,1), c(1), 2L)


