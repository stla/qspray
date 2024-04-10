test_that("show", {
  set.seed(3141L)
  q1 <- rQspray()
  expect_identical(Print(q1), "-2*x^4y^3z^4 - 4*y^2z^2 ")
  showQsprayOption(q1, "x") <- "A"
  expect_identical(Print(q1), "-2*A1^4.A2^3.A3^4 - 4*A2^2.A3^2 ")
  q2 <- rQspray()
  expect_identical(
    Print(q1+q2), 
    "-2*A1^4.A2^3.A3^4 + 5*A1^4.A2 - A2^2.A3^3 - 4*A2^2.A3^2 - 3 "
  )
  
  q3 <- 1 + qlone(1) 
  expect_identical(Print(q3), "x + 1 ")
  showQsprayOption(q3, "x") <- "A"
  expect_identical(Print(q3), "A + 1 ")
  expect_identical(Print(q3 + qlone(1)), "2*A + 1 ")
  expect_identical(Print(q3 + qlone(2)), "x + y + 1 ")
})