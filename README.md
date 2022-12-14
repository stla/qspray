The ‘qspray’ package
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/qspray/workflows/R-CMD-check/badge.svg)](https://github.com/stla/qspray/actions)
[![R-CMD-check-valgrind](https://github.com/stla/qspray/actions/workflows/R-CMD-check-valgrind.yaml/badge.svg)](https://github.com/stla/qspray/actions/workflows/R-CMD-check-valgrind.yaml)
<!-- badges: end -->

*R package for multivariate polynomials with rational coefficients.*

This package is strongly inspired by Robin Hankin’s **spray** package.
The C++ implementations are very similar.

``` r
library(qspray)
```

The easiest way to define a multivariate polynomial with **qspray** is
to start by introducing the generating variables with the help of the
`qlone` function and then to combine them with arithmetic operations:

``` r
x <- qlone(1); y <- qlone(2); z <- qlone(3)
pol <- 4*x^2 + "1/2"*y - 5*x*y*z
pol
## -5*x^(1, 1, 1) + 4*x^(2) + 1/2*x^(0, 1)
```

Or maybe you prefer to define the polynomial by giving it as a string:

``` r
qsprayMaker(string = "4 x^(2) + 1/2 x^(0, 1) - 5 x^(1, 1, 1)")
## -5*x^(1, 1, 1) + 1/2*x^(0, 1) + 4*x^(2)
```

As you want, but this method is not highly robust.

Some arithmetic on this polynomial:

``` r
-pol
## -5*x^(1, 1, 1) + 4*x^(2) + 1/2*x^(0, 1)
2 * pol
## -10*x^(1, 1, 1) + 8*x^(2) + x^(0, 1)
pol / 2
## -5/2*x^(1, 1, 1) + 2*x^(2) + 1/4*x^(0, 1)
"5/3" * pol
## -25/3*x^(1, 1, 1) + 20/3*x^(2) + 5/6*x^(0, 1)
pol + 5
## 1/2*x^(0, 1) + 5*x^() + 4*x^(2) - 5*x^(1, 1, 1)
pol - "2/5"
## 1/2*x^(0, 1) - 2/5*x^() + 4*x^(2) - 5*x^(1, 1, 1)
pol^2
## 25*x^(2, 2, 2) + 16*x^(4) - 40*x^(3, 1, 1) + 1/4*x^(0, 2) + 4*x^(2, 1) - 5*x^(1, 2, 1)
```

Two polynomials can be added and multiplied:

``` r
pol1 <- pol
pol2 <- pol
pol1 + pol2
## x^(0, 1) + 8*x^(2) - 10*x^(1, 1, 1)
pol1 - pol2
## 0
pol1 * pol2
## 25*x^(2, 2, 2) - 40*x^(3, 1, 1) + 16*x^(4) - 5*x^(1, 2, 1) + 4*x^(2, 1) + 1/4*x^(0, 2)
```

Use `evalQspray` to evaluate a polynomial for some values of the
variables:

``` r
evalQspray(pol, c("1", "2", "3/2"))
## Big Rational ('bigq') :
## [1] -10
```

Alternatively, you can convert the polynomial to a function:

``` r
f <- as.function(pol)
f("1", "2", "3/2")
## [1] "-10"
```

You can pass the strings you want as the arguments of this functions:

``` r
f("x", "y", "z")
## [1] "(8*x^2-10*x*y*z+y)/2"
```

The package also provides a function which returns the exact value of
the integral of a polynomial with rational coefficients over a simplex
with rational vertices:

``` r
# variables
x <- qlone(1); y <- qlone(2); z <- qlone(3)
# polynomial
P <- x^4 + y + 2*x*y^2 - 3*z
# simplex (tetrahedron) vertices
v1 <- c(1, 1, 1)
v2 <- c(2, 2, 3)
v3 <- c(3, 4, 5)
v4 <- c(3, 2, 1)
# simplex
S <- rbind(v1, v2, v3, v4)
# integral
integratePolynomialOnSimplex(P, S)
## Big Rational ('bigq') :
## [1] 1387/42
```
