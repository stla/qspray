The ‘qspray’ package
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/qspray/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/qspray/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check-valgrind](https://github.com/stla/qspray/actions/workflows/R-CMD-check-valgrind.yaml/badge.svg)](https://github.com/stla/qspray/actions/workflows/R-CMD-check-valgrind.yaml)
<!-- badges: end -->

*R package for multivariate polynomials with rational coefficients.*

------------------------------------------------------------------------

This package is strongly inspired by Robin Hankin’s **spray** package.
The C++ implementations are very similar.

``` r
library(qspray)
```

The **qspray** package provides the `qspray` objects, which represent
multivariate polynomials whose coefficients are rational numbers.

## Creating a `qspray` and arithmetic

The easiest way to build a multivariate polynomial with **qspray** is to
start by introducing the generating variables with the help of the
`qlone` function and then to combine them with arithmetic operations:

``` r
x <- qlone(1); y <- qlone(2); z <- qlone(3)
( pol <- 4*x^2 + "1/2"*y - 5*x*y*z/3 )
## 4*x^2 - 5/3*x.y.z + 1/2*y
```

I often like to use a function like this:

``` r
f <- function(x, y, z) {
  4*x^2 + y/2 - 5*x*y*z/3
}
f(x, y, z)
## 4*x^2 - 5/3*x.y.z + 1/2*y
```

Or maybe you prefer to define the polynomial by giving it as a string:

``` r
qsprayMaker(string = "4 x^(2) + 1/2 x^(0, 1) - 5/3 x^(1, 1, 1)")
## 4*x^2 - 5/3*x.y.z + 1/2*y
```

As you want, but this method is not highly robust.

Some arithmetic on this polynomial:

``` r
-pol
## -4*x^2 + 5/3*x.y.z - 1/2*y
2 * pol
## 8*x^2 - 10/3*x.y.z + y
pol / 2
## 2*x^2 - 5/6*x.y.z + 1/4*y
"5/3" * pol
## 20/3*x^2 - 25/9*x.y.z + 5/6*y
pol + 5
## 4*x^2 - 5/3*x.y.z + 1/2*y + 5
pol - gmp::as.bigq("2/5")
## 4*x^2 - 5/3*x.y.z + 1/2*y - 2/5
pol^2
## 16*x^4 - 40/3*x^3.y.z + 25/9*x^2.y^2.z^2 + 4*x^2.y - 5/3*x.y^2.z + 1/4*y^2
```

Two polynomials can be added and multiplied:

``` r
pol1 <- pol
pol2 <- pol
pol1 + pol2
## 8*x^2 - 10/3*x.y.z + y
pol1 - pol2
## 0
pol1 * pol2
## 16*x^4 - 40/3*x^3.y.z + 25/9*x^2.y^2.z^2 + 4*x^2.y - 5/3*x.y^2.z + 1/4*y^2
```

## Evaluating a `qspray`

Use `evalQspray` to evaluate a polynomial for some values of the
variables:

``` r
evalQspray(pol, c("1", "2", "3/2"))
## Big Rational ('bigq') :
## [1] 0
```

Alternatively, you can convert the polynomial to a function:

``` r
g <- as.function(pol)
g("1", "2", "3/2")
## [1] "0"
```

You can pass the strings you want as the arguments of this function:

``` r
g("x", "y", "z")
## [1] "(24*x^2-10*x*z*y+3*y)/6"
```

If you want numerical approximations in the results, use the option
`N=TRUE`:

``` r
h <- as.function(pol, N = TRUE)
h("1", "2", "3/2")
## [1] 2e-29
h("x", "y", "z")
## expression(4 * x^2 - 1.6666666666 * x * y * z + 0.5 * y)
```

You can also perform “partial evaluation” of a `qspray`, that is to say
replacing only certain variables. This is done by using the function
`substituteQspray` and indicating the variables to be kept with `NA`:

``` r
substituteQspray(pol, c("1", NA, "3/2"))
## -2*y + 4
f(gmp::as.bigq(1), y, gmp::as.bigq("3/2"))
## -2*y + 4
g("1", "y", "3/2")
## [1] "2*(2-y)"
h("1", "y", "3/2")
## expression(-2 * y + 4)
```

## Showing a `qspray`

You can control the way of printing a `qspray` with the help of the
function `showQsprayOption<-`. By default, the monomials of a `qspray`
are printed in the style of `x^2.y.z^3` if there are at most three
variables, otherwise they are printed like `x1^2.x2.x3^3`:

``` r
set.seed(3141)
( qspray <- rQspray() ) # a random qspray
## -2*x^4.y^3.z^4 - 4*y^2.z^2
qspray + qlone(4)^99
## -2*x1^4.x2^3.x3^4 - 4*x2^2.x3^2 + x4^99
```

If you want to always use the second way, you can do:

``` r
showQsprayOption(qspray, "x") <- "x"
qspray
## -2*x1^4.x2^3.x3^4 - 4*x2^2.x3^2
```

If you want to restore the way `qspray` objects were printed in previous
versions, you can do

``` r
showQsprayOption(qspray, "showMonomial") <- showMonomialOld()
qspray
## -2*x^(4, 3, 4) - 4*x^(0, 2, 2)
```

There are three possible show options that can be passed to
`showQsprayOption`:

- The most general show option is `"showQspray"`. A `showQspray`
  function, that is to say a function appropriate for the `"showQspray"`
  option, must be a function which transforms a `qspray` to a string.
  The package provides some helper functions to built such functions,
  like `showQsprayXYZ` and `showQsprayX1X2X3`. With `showQsprayXYZ`, you
  can choose the letters you want to denote the variables:

``` r
f <- showQsprayXYZ(c("A", "B", "C"))
f(qspray)
## [1] "-2*A^4.B^3.C^4 - 4*B^2.C^2"
```

With `showQsprayX1X2X3`, you choose only one letter for the variables
which will be indexed:

``` r
f <- showQsprayX1X2X3("X")
f(qspray)
## [1] "-2*X1^4.X2^3.X3^4 - 4*X2^2.X3^2"
```

Once you have constructed such a function, you pass it as a show option
by doing `showQsprayOption(qspray, "showQspray") <- f`.

- The second possible show option is `"showMonomial"`, to control the
  way the monomials are printed. Actually in the two above examples of
  `showQsprayXYZ` and `showQsprayX1X2X3` we only changed the way the
  monomials are printed. Indeed, these two commands are equivalent:

``` r
showQsprayOption(qspray, "showQspray") <- showQsprayXYZ(c("A", "B", "C"))
showQsprayOption(qspray, "showMonomial") <- showMonomialXYZ(c("A", "B", "C"))
```

and these two commands are equivalent as well:

``` r
showQsprayOption(qspray, "showQspray") <- showQsprayX1X2X3("X")
showQsprayOption(qspray, "showMonomial") <- showMonomialX1X2X3("X")
```

But the `showQspray` functions allow finer control, e.g. they allow to
control the multiplication symbol which separates a coefficient and a
monomial within a term.

- Finally there is the show option `"x"`. Setting this option to a
  letter `x`:

``` r
showQsprayOption(qspray, "x") <- x
```

is equivalent to:

``` r
showQsprayOption(qspray, "showMonomial") <- showMonomialX1X2X3(x)
```

But `showMonomialX1X2X3` also allows to control the way the individual
powers are collapsed, e.g. `"x^2.y.z^3"` (the default) or `"x^2*y*z^3"`,
or `"x^2yz^3"`. If the dot is nice for you, use the `"x"` option, that’s
less code to type.

By the way, a `qspray` object is an S4 object with two slots: `powers`
and `coeffs`. The `powers` slot is a list of vector of exponents and the
`coeffs` slot is a character vector, whose each element is coercable to
a `bigq` number by an application of the function `gmp::as.bigq`. The
`showMonomial` functions act only on the `powers` slot.

When an arithmetic operation is performed between two `qspray` objects,
the show options of the first one are passed to the result, *if
possible*:

``` r
qspray + qlone(4)^99
## -2*X1^4.X2^3.X3^4 - 4*X2^2.X3^2 + X4^99
```

For example, this is not possible if you specify only three letters for
the variables and you perform an operation with a `qspray` involving the
fourth variable:

``` r
showQsprayOption(qspray, "showMonomial") <- showMonomialXYZ(c("a", "b", "c"))
qspray
## -2*a^4.b^3.c^4 - 4*b^2.c^2
qspray + qlone(4)^99
## -2*a1^4.a2^3.a3^4 - 4*a2^2.a3^2 + a4^99
```

## Exact integration over a simplex

The package provides a function which returns the exact value of the
integral of a polynomial with rational coefficients over a simplex whose
vertices have rational Cartesian coordinates:

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

## Transforming a `qspray`

You can get a derivative of a `qspray`:

``` r
derivQspray(P, i = 2) # derivative w.r.t y
## 4*x.y + 1
```

You can permute the variables of a `qspray`:

``` r
swapVariables(P, 1, 3)
## -3*x + 2*y^2.z + y + z^4
```

## Gröbner bases

Finally, let us mention the `groebner` function, which computes a
Gröbner basis of the ideal generated by a list of polynomials:

``` r
f <- qsprayMaker(string = "x^(3) - 2 x^(1,1)")
g <- qsprayMaker(string = "x^(2,1) - 2 x^(0,2) + x^(1)")
groebner(list(f, g))
## [[1]]
## x - 2*y^2 
## 
## [[2]]
## y^3
```

## Packages using ‘qspray’

There are packages depending on the **qspray** package (some of them are
not on CRAN yet):

- [jack](https://github.com/stla/jackR) (Jack polynomials)
- [resultant](https://github.com/stla/resultant) (resultant,
  subresultants, and greatest common divisor of two `qspray` objects)
- [ratioOfQsprays](https://github.com/stla/ratioOfQsprays) (fractions of
  `qspray` objects)
- [symbolicQspray](https://github.com/stla/symbolicQspray) (multivariate
  polynomials whose coefficients are `ratioOfQsprays` objects)
