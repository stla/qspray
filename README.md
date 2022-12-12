The ‘qspray’ package
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/qspray/workflows/R-CMD-check/badge.svg)](https://github.com/stla/qspray/actions)
<!-- badges: end -->

*R package for multivariate polynomials with rational coefficients.*

``` r
library(qspray)
```

Define a polynomial:

``` r
x <- lone(1); y <- lone(2); z <- lone(3)
pol <- 4*x^2 + "1/2"*y - 5*x*y*z
pol
## -5x^(1, 1, 1) + 4x^(2) + 1/2x^(0, 1)
```

Some arithmetic on this polynomial:

``` r
-pol
## -5x^(1, 1, 1) + 4x^(2) + 1/2x^(0, 1)
2 * pol
## -10x^(1, 1, 1) + 8x^(2, 0) + x^(0, 1)
# TODO: pol / 2
"5/3" * pol
## -25/3x^(1, 1, 1) + 20/3x^(2, 0) + 5/6x^(0, 1)
pol + 5
## 1/2x^(0, 1) + 5x^() + 4x^(2) - 5x^(1, 1, 1)
pol - "2/5"
## 1/2x^(0, 1) - 2/5x^() + 4x^(2) - 5x^(1, 1, 1)
pol^2
## 25x^(2, 2, 2) + 16x^(4, 0, 0) - 40x^(3, 1, 1) + 4x^(2, 1, 0) + 1/4x^(0, 2, 0) - 5x^(1, 2, 1)
```

Two polynomials can be added and multiplied:

``` r
pol1 <- pol
pol2 <- pol
pol1 + pol2
## x^(0, 1) + 8x^(2) - 10x^(1, 1, 1)
pol1 - pol2
## 0x^()
pol1 * pol2
## 25x^(2, 2, 2) - 40x^(3, 1, 1) + 16x^(4, 0) - 5x^(1, 2, 1) + 4x^(2, 1) + 1/4x^(0, 2)
```

Use `evalQspray` to evaluate a polynomial for some values of the
variables:

``` r
evalQspray(pol, c("1", "2", "3/2"))
## Big Rational ('bigq') :
## [1] -10
```
