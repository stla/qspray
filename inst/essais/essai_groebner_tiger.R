library(qspray)

cosu <- qlone(1)
sinu <- qlone(2)
cosv <- qlone(3)
sinv <- qlone(4)
cosw <- qlone(5)
sinw <- qlone(6)

nvariables <- 6
parameters <- NULL
equations <- list(
  "x1" = cosu * (2 + cosw),
  "x2" = sinu * (2 + cosw),
  "x3" = cosv * (2 + sinw),
  "x4" = sinv * (2 + sinw)
)
relations <- list(
  cosu^2 + sinu^2 - 1,
  cosv^2 + sinv^2 - 1,
  cosw^2 + sinw^2 - 1
)

implicitization(nvariables, parameters, equations, relations)
