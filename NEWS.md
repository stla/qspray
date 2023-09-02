# qspray 2.1.0.9000

- Unit tests.

- The polynomials are now printed with the lexicographic order of the monomials.


# qspray 2.1.0

- New function `implicitization`, to transform a system of parametric equations to an implicit equation. This is based on Gröbner bases. This doesn't always 
work (use the 'giacR' package if needed).

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
