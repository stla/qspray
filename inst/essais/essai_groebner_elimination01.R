library(qspray)
library(Ryacas)

x <- qlone(1); y <- qlone(2); z <- qlone(3)
f1 <- x^2 + y + z - 1
f2 <- x + y^2 + z - 1
f3 <- x + y + z^2 - 1

gb <- groebner(list(f1, f2, f3))
( gbstrings <- lapply(gb, prettyQspray, vars = c("x", "y", "z")) )

yac_str("Factor(z^6 - z^2 + 4*z^3 - 4*z^4)")
yac_str("factor1 := 2*z + z^2 - 1")
yac_str("roots := PSolve(2*z + z^2 - 1, z)")
yac_str("root1 := Nth(roots, 1)") # Sqrt(2) - 1
yac_str("root2 := Nth(roots, 2)") # -Sqrt(2) - 1

yac_str(sprintf("p3 := %s", gbstrings[[3L]]))
yac_str("PSolve(Eliminate(z, Sqrt(2)-1, p3), y)")
yac_str("N(PSolve(Eliminate(z, Sqrt(2)-1, p3), y))")

yac_str(sprintf("p2 := %s", gbstrings[[2L]]))
yac_str("PSolve(Eliminate(z, Sqrt(2)-1, p2), y)")
yac_str("N(PSolve(Eliminate(z, Sqrt(2)-1, p2), y))") # 2nd root is not root of p3

yac_str(sprintf("p1 := %s", gbstrings[[1L]]))
yac_str(
  "PSolve(Eliminate(y, Sqrt(2)-1, Eliminate(z, Sqrt(2)-1, p1)), x)"
)
yac_str(
  "N(PSolve(Eliminate(y, Sqrt(2)-1, Eliminate(z, Sqrt(2)-1, p1)), x))"
)
