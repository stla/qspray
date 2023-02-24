library(qspray)
x <- qlone(1)
y <- qlone(2)
qspray <- 2*x  + 3*x*y
derivQspray(qspray, 3) # 0

derivQspray(qspray, 2, 1) # 3*x^(1)
dQspray(qspray, c(0, 1)) # 3*x^(1)

derivQspray(qspray, 1, 2) # 0
derivQspray(qspray, 2, 2) # 0
 
qspray <- x + 2*y  + 3*x*y
dQspray(qspray, c(1, 1)) # 3*x^()
derivQspray(derivQspray(qspray, 1), 2) # 3*x^()


library(spray)
x <- lone(1, 2)
y <- lone(2, 2)
S <- 2*x + 3*x*y
deriv(S, 2, 2)
