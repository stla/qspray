library(qspray)

#' @title Implicitization
#' @description Implicitization of a system of parametric equations.
#'
#' @param nvariables number of variables
#' @param parameters character vector of the names of the parameters, or 
#'   \code{NULL} if there's no parameter
#' @param equations list of qspray polynomials representing the parametric 
#'   equations
#' @param relations list of qspray polynomials representing the relations 
#'   between the variables and the parameters, or \code{NULL} if there is none
#'
#' @return A list of qspray polynomials.
#' @export
#'
#' @examples
#' library(qspray)
#' # ellipse example ####
#' # variables 
#' cost <- qlone(1)
#' sint <- qlone(2)
#' # parameters
#' a <- qlone(3)
#' b <- qlone(4)
#' #
#' nvariables <- 2
#' parameters <- c("a", "b")
#' equations <- list(
#'   "x" = a * cost,
#'   "y" = b * sint
#' )
#' relations <- list(
#'   cost^2 + sint^2 - 1
#' )
#' # 
#' eqs <- implicitization(nvariables, parameters, equations, relations)
implicitization <- function(nvariables, parameters, equations, relations) {
  stopifnot(isPositiveInteger(nvariables))
  stopifnot(is.null(parameters) || isStringVector(parameters))
  stopifnot(is.list(equations), length(equations) > 1L)
  stopifnot(is.null(relations) || is.list(relations))
  #
  nequations <- length(equations)
  nrelations <- length(relations)
  nqlone <- max(vapply(equations, arity, integer(1L)))
  coordinates <- lapply((nqlone+1L):(nqlone+nequations), qlone)
  generators <- relations
  for(i in seq_along(equations)) {
    generators <- append(generators, coordinates[[i]] - equations[[i]])
  }
  #
  gb <- groebner(generators)
  isfree <- function(i) {
    all(vapply(gb[[i]]@powers, function(pows) {
      length(pows) == 0L || 
        (length(pows > nvariables) && all(pows[1L:nvariables] == 0L))
    }, logical(1L)))
  }
  free <- c(FALSE, vapply(2L:length(gb), isfree, logical(1L)))
  #
  results <- gb[free]
  for(i in seq_along(results)) {
    el <- results[[i]]
    coeffs <- el@coeffs
    powers <- el@powers
    for(j in seq_along(powers)) {
      powers[[j]] <- tail(powers[[j]], -nvariables)
    }
    results[[i]] <- qsprayMaker(powers, coeffs)
  }
  vars <- c(parameters, names(equations))
  messages <- lapply(results, prettyQspray, vars = vars)
  for(msg in messages) {
    message(msg)
  }
  invisible(results)
}

#
# variables
cost <- qlone(1)
sint <- qlone(2)
# parameters
a    <- qlone(3)
b    <- qlone(4)
#
nvariables <- 2
parameters <- c("a", "b")
equations <- list(
  "x" = a * cost,
  "y" = b * sint
)
relations <- list(
  cost^2 + sint^2 - 1
)

# Enneper 
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

# satellite curve
cost  <- qlone(1)
sint  <- qlone(2)
cos2t <- qlone(3)
sin2t <- qlone(4)
A     <- qlone(5)
B     <- qlone(6)
nvariables <- 4
parameters <- c("A", "B")
equations <- list(
  "x" = A*cost*cos2t - sint*sin2t,
  "y" = A*sint*cos2t + cost*sin2t,
  "z" = B*cos2t
)
relations <- list(
  cost^2 + sint^2 - 1,
  cos2t - cost^2+sint^2, 
  sin2t - 2*sint*cost,
  A^2 + B^2 - 1
)

# implicitization(variables, parameters, equations, relations)


