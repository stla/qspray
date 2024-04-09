rox <- function() {
  roxygen2::roxygenise()
}

myinstall <- function() {
  try(pkgload::unload("qspray"))
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "devtools::install(quick = TRUE, keep_source = TRUE)"
    )
  } else {
    devtools::install(quick = TRUE, keep_source = TRUE)
  }
}

mydocument <- function() {
  try(pkgload::unload("qspray"))
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "roxygen2::roxygenise(load_code = roxygen2::load_installed)" 
    )
  } else {
    roxygen2::roxygenise(load_code = roxygen2::load_installed)
  }
}
