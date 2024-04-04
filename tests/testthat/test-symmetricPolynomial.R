test_that("Symmetric polynomial", {
  e1 <- ESFpoly(4, 1)
  e2 <- ESFpoly(4, 2)
  e3 <- ESFpoly(4, 3)
  f <- e1 + 2*e2 + 3*e3
  expect_true(isSymmetricPolynomial(f))
})

test_that("MSPcombination", {
  qspray <- PSFpoly(4, c(3, 1)) + ESFpoly(4, c(2, 2)) + 4L
  expect_no_error(MSPcombination(qspray, check = TRUE))
})