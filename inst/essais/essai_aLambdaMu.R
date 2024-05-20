library(qspray)
library(syt)
library(partitions)

aLambdaMu <- function(lambda, mu) {
  partitions <- parts(sum(lambda))
  KN1 <- apply(partitions, 2L, function(kappa) {
    KostkaNumber(kappa, lambda)
  })
  KN2 <- apply(partitions, 2L, function(kappa) {
    KostkaNumber(kappa, mu)
  })
  sum(KN1 * KN2)
}

bLambdaMu <- function(lambda, mu) {
  partitions <- parts(sum(lambda))
  KN1 <- apply(partitions, 2L, function(kappa) {
    KostkaNumber(kappa, lambda)
  })
  KN2 <- apply(partitions, 2L, function(kappa) {
    KostkaNumber(conjugate(kappa), mu)
  })
  sum(KN1 * KN2)
}

aLambdaMu(c(4, 1, 1), c(3, 3))
HallInnerProduct(CSHFpoly(6, c(4, 1, 1)), CSHFpoly(6, c(3, 3)))

bLambdaMu(c(1, 1, 1, 1, 1, 1), c(3, 3))
HallInnerProduct(CSHFpoly(6, c(1, 1, 1, 1, 1, 1)), ESFpoly(6, c(3,3)))


#          for j=2:m
#              for i=1:j-1
#                  Ci(i,j) = -Ci(i,i:j-1)*Ci(i:j-1,j);
#              end
#          end
km <- jack::KostkaNumbers(4, "2")
library(gmp)
Ci <- as.bigq(km)
for(j in 2:5) {
  for(i in 1:(j-1)) {
    Ci[i, j] <- -sum(c(Ci[i, i:(j-1)]) * c(Ci[i:(j-1), j]))
  }
}
Ci
RationalMatrix::Qinverse(km)

invTriMatrix <- function(A) {
  d <- nrow(A)
  if(d == 1L) {
    return(as.matrix(1/A[1,1]))
  } else {
    B <- invTriMatrix(A[1L:(d-1L), 1L:(d-1L)])
    newColumn <- as.bigq(integer(d-1L))
    for(i in 1L:(d-1L)) {
      newColumn[i] <- -sum(c(B[i, i:(d-1L)]) * c(A[i:(d-1L), d])) / A[d,d]
    }
    newRow <- as.bigq(integer(d))
    newRow[d] <- 1/A[d,d]
    B <- rbind(cbind(B, newColumn), newRow)
    return(B)
  }
}
km <- jack::KostkaNumbers(4, "2")
invTriMatrix(as.bigq(km)*2)
RationalMatrix::Qinverse(as.bigq(km)*2)
