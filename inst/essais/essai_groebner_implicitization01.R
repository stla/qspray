library(qspray)

cost <- qlone(1)
sint <- qlone(2)
a    <- qlone(3)
b    <- qlone(4)

variables <- c("cost", "sint")
parameters <- c("a", "b")
equations <- list(
  "x" = a * cost,
  "y" = b * sint
)
relations <- list(
  cost^2 + sint^2 - 1
)

nequations <- length(equations)
nrelations <- length(relations)
n1 <- max(vapply(equations, qspray:::arity, integer(1L))) # 4
coordinates <- lapply((n1+1L):(n1+nequations), qlone)

generators <- relations
for(i in seq_along(equations)) {
  generators <- append(generators, coordinates[[i]] - equations[[i]])
}

gb <- groebner(generators)

nvariables <- length(variables)
isfree <- function(i) {
  all(vapply(gb[[i]]@powers, function(pows) {
    length(pows) == 0L || 
      (length(pows > nvariables) && all(pows[1L:nvariables] == 0L))
  }, logical(1L)))
}

free <- c(FALSE, vapply(2L:length(gb), isfree, logical(1L)))

results <- gb[free]
vars <- c(parameters, names(equations))
for(i in seq_along(results)) {
  el <- results[[i]]
  coeffs <- el@coeffs
  powers <- el@powers
  for(j in seq_along(powers)) {
    powers[[j]] <- tail(powers[[j]], -nvariables)
  }
  results[[i]] <- qsprayMaker(powers, coeffs)
}

lapply(results, prettyQspray, vars = vars)

