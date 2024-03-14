library(qspray)
x <- qlone(1)
y <- qlone(2)
z <- qlone(3)
p <- x^2 + y^2 + x*y*z - 1
qspray:::substituteQspray(p, c(NA, NA, "3/2"))