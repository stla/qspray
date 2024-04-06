# qspray 2.2.0

- The `show` method has been changed. The monomial previously printed as `x^(2,0,3)` is not printed as `x1^2.x3^3`.

- New helper functions to construct a custom `show` method.

- New function `substituteQspray`, to substitute some variables in a polynomial.

- New function `composeQspray`, to get the polynomial obtained by substituting the variables of a polynomial with polynomials.
 
- New function `isSymmetricPolynomial`, to check whether a `qspray` polynomial is symmetric. It is inefficient, but useful for testing purpose. Use `MSPcombination` for efficiency.

- There was an error in `qdivision` (without severe consequences).

- The function `as.function.qspray` now has a Boolean argument `N`; if set to `TRUE`, the function returns a numerical approximation of the result.

- New function `PSFpoly`, which computes the power sum polynomials.

- New function `MSPcombination`, to get a symmetric polynomial as a linear combination of the monomial symmetric polynomials. It is also useful to check whether a polynomial is symmetric.

- New function `prettySymmetricQspray` which prints a symmetric polynomial as a linear combination of the monomial symmetric polynomials.

- New function `sprayDivision`, returning the quotient and the remainder of the division of a polynomial by a polynomial.

- New function `HallInnerProduct`, which computes the Hall inner product between symmetric polynomials. Very inefficient (will be hopefully improved in the future).


# qspray 2.1.1

- Unit tests.

- The polynomials are now printed with the lexicographic order of the monomials.


# qspray 2.1.0

- New function `implicitization`, to transform a system of parametric equations to an implicit equation. This is based on Gröbner bases. This doesn't always work (use the 'giacR' package if needed).

- The division step in the Buchberger algorithm is now performed in C++.

- The `groebner` function has been considerably improved, but it still can be very slow.


# qspray 2.0.0

- The implementation of the opposite polynomial was wrong (it did nothing).

- The pretty form of a polynomial returned by `prettyQspray` has been improved: now there are spaces around '+' and '-'.

- The division of a polynomial by a `bigq` number did not always work.

- New function `groebner`, to compute a Gröbner basis.


# qspray 1.1.1

Fixed a small error in the C++ code.


# qspray 1.1.0

Differentiation of `qspray` polynomials (functions `derivQspray` and `dQspray`).


# qspray 1.0.0

The package does not depend on 'RcppCGAL' neither 'RcppEigen' anymore. 
It now imports 'RationalMatrix' (to compute the determinant of a rational 
matrix).


# qspray 0.1.1

Removed a useless function in the C++ code.


# qspray 0.1.0

First release.
