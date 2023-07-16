# Enneper 
library(qspray)
u <- qlone(1)
v <- qlone(2)
nvariables <- 2
parameters <- NULL
equations <- list(
  "x" = 3*u + 3*u*v^2 - u^3,
  "y" = 3*v + 3*u^2*v - v^3,
  "z" = 3*u^2 - 3*v^2
)
relations <- NULL

implicitization(nvariables, parameters, equations, relations)
