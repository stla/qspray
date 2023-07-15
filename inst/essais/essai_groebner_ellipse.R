library(qspray)
cost <- qlone(1)
sint <- qlone(2)
x <- qlone(3)
y <- qlone(4)
p1 <- cost^2 + sint^2 - 1
p2 <- x - 2*cost
p3 <- y - 3*sint
( gb <- groebner(list(p1, p2, p3)) )

# 
cost <- qlone(1)
sint <- qlone(2)
a <- qlone(3)
x <- qlone(4)
y <- qlone(5)

p1 <- cost^2 + sint^2 - 1
p2 <- x - a*cost
p3 <- y - 3*sint

( gb <- groebner(list(p1, p2, p3)) )

# input: 2L = nb de variables t_i
#  1L = nb de constantes (ici a)
# ... non
# input: equations = (a*cost, 3*sint), relations ('constraints') p1 (au moins une)
#  et le user doit mettre les t_i (ici cost et sint) en premiers qlone
#  problème si relation entre les constantes ?.. je ne crois pas
prettyQspray(gb[[5]], vars = c("degage1", "degage2", "a", "x", "y"))

vapply(gb[[5]]@powers, function(pows) {
  length(pows) == 0L || (length(pows > 2L) && all(pows[1L:2L] == 0L))
}, logical(1L))


# ça me semble plus simple de mettre les t_i à la fin :
a <- qlone(1)
x <- qlone(2)
y <- qlone(3)
cost <- qlone(4)
sint <- qlone(5)

p1 <- cost^2 + sint^2 - 1
p2 <- x - a*cost
p3 <- y - 3*sint

( gb <- groebner(list(p1, p2, p3)) )

# => ça ne marche pas !!