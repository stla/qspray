library(ratioOfQsprays)

set.seed(3141)
( q1 <- rQspray() )
( q2 <- rQspray() )
q1 + q2
q1 / q2

showQsprayOption(q1, "x") <- "A"
q1 + q2
q1 / q2

( q1 <- qlone(1) )
showQsprayOption(q1, "x") <- "A"
q1 / (1 + q1)
q1 + q2
q1 / q2

( q1 <- qlone(1) + qlone(2) )
showQsprayOption(q1, "x") <- "A"
q1
q1 + q2
q1 / q2
