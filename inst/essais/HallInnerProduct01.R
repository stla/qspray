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

HallInnerProduct(PSFpoly(4, c(1,1,1,1)), PSFpoly(4, c(1,1,1,1)), alpha = t)
24 * t^4

HallInnerProduct(
  jack::JackPol(3, c(2,1), alpha), PSFpoly(3, c(2,1)), alpha = t
)
t + 2 # no
