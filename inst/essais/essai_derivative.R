library(qspray)
x <- qlone(1)
y <- qlone(2)
qspray <- 2*x  + 3*x*y
derivQspray(qspray, 3)

derivQspray(qspray, 1, 2)
derivQspray(qspray, 2, 2)
 
qspray <- x + 2*y  + 3*x*y
dQspray(qspray, c(1, 1))
derivQspray(derivQspray(qspray, 1), 2)


