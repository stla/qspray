library(qspray)
x <- qlone(1)
y <- qlone(2)
qspray <- 2*x  + 3*x*y
derivQspray(qspray, 1)

derivQspray(qspray, 1, 2)
derivQspray(qspray, 2, 2)
 