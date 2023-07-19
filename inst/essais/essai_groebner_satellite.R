library(qspray)

# satellite curve
cost  <- qlone(1)
sint  <- qlone(2)
cos2t <- qlone(3)
sin2t <- qlone(4)
A     <- qlone(5)
B     <- qlone(6)
nvariables <- 4
parameters <- c("A", "B")
equations <- list(
  "x" = A*cost*cos2t - sint*sin2t,
  "y" = A*sint*cos2t + cost*sin2t,
  "z" = B*cos2t
)
relations <- list(
  cost^2 + sint^2 - 1,
  cos2t - cost^2+sint^2, 
  sin2t - 2*sint*cost,
  A^2 + B^2 - 1
)


implicitization(nvariables, parameters, equations, relations)
