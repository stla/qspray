test_that("Symmetric polynomial", {
  e1 <- ESFpoly(3, 1)
  e2 <- ESFpoly(3, 2)
  e3 <- ESFpoly(3, 3)
  f <- e1 + 2*e2 + 3*e3
  expect_true(isSymmetricPolynomial(f))
})