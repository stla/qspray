# # Enneper 

library(qspray)
u <- qlone(4); v <- qlone(5); x <- qlone(1); y <- qlone(2); z <- qlone(3)
f1 <- x - (3*u + 3*u*v^2 - u^3)
f2 <- y - (3*v + 3*u^2*v - v^3)
f3 <- z - (3*u^2 - 3*v^2)
gb <- groebner(list(f1, f2, f3))
gb
