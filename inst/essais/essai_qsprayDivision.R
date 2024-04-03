library(qspray)
x <- qlone(1)
y <- qlone(2)
z <- qlone(3)

B <- x*y^2 + z*x^2 + 1
Q <- x^2*y^2*z^2 - 3
R <- x*y
A <- B*Q + R

divis <- qspray:::qsprayDivisionRcpp(
  A@powers, A@coeffs, B@powers, B@coeffs, 3L
) 

qspray:::qspray_from_list(divis[["Q"]])
qspray:::qspray_from_list(divis[["R"]])