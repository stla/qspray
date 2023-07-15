library(qspray)
x <- qlone(1)
y <- qlone(2)
cost <- qlone(3)
sint <- qlone(4)
a <- qlone(5)

p1 <- x - a*cost
p2 <- y - 3*sint
p3 <- cost^2 + sint^2 - 1

groebner(list(p1, p2, p3))
