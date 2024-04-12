test_that("Groebner basis", {
  x <- qlone(1); y <- qlone(2); z <- qlone(3)
  p1 <- x^2 + y + z - 1
  p2 <- x + y^2 + z - 1
  p3 <- x + y + z^2 - 1
  #
  gb <- groebner(list(p1, p2, p3))
  #
  expect_length(gb, 4L)
  #
  f1 <- as.function(gb[[1L]], N = TRUE)
  f2 <- as.function(gb[[2L]], N = TRUE)
  f3 <- as.function(gb[[3L]], N = TRUE)
  f4 <- as.function(gb[[4L]], N = TRUE)
  #
  x1 <- y1 <- z1 <- sqrt(2) - 1
  expect_equal(f1(x1, y1, z1), 0)
  expect_equal(f2(x1, y1, z1), 0)
  expect_equal(f3(x1, y1, z1), 0)
  expect_equal(f4(x1, y1, z1), 0)
  #
  x2 <- y2 <- z2 <- -sqrt(2) - 1
  expect_equal(f1(x2, y2, z2), 0)
  expect_equal(f2(x2, y2, z2), 0)
  expect_equal(f3(x2, y2, z2), 0)
  expect_equal(f4(x2, y2, z2), 0)
})

test_that("isPolynomialOf", {
  P <- function(X, Y) X^2*Y + 2*X + 3
  x <- qlone(1); y <- qlone(2); z <- qlone(3)
  q1 <- x + y + z
  q2 <- x^2*z^2 + 4*y + 1
  qspray <- P(q1, q2)
  check <- isPolynomialOf(qspray, list(q1, q2))
  expect_true(check)
  POLYNOMIAL <- attr(check, "polynomial")
  expect_true(POLYNOMIAL == P(qlone(1), qlone(2)))
  expect_true(composeQspray(POLYNOMIAL, list(q1, q2)) == qspray)
})

