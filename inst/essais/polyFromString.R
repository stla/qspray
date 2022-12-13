stringToPol <- function(p){
  p <- gsub("\\)\\s*-\\s*(\\d*/*\\d*)\\s*", ")+-\\1", p)
  p <- gsub("^-\\s*x", "-1x", trimws(p, "left"))
  terms <- strsplit(p, "+", fixed = TRUE)[[1L]]
  csts <- !grepl("x", terms)
  terms[csts] <- paste0(terms[csts], "x^(0")
  ss <- purrr::transpose(strsplit(terms, "x^(", fixed = TRUE))
  coeffs <- trimws(unlist(ss[[1L]], recursive = FALSE))
  coeffs[coeffs == ""] <- "1"
  powers <- sub(")", "", unlist(ss[[2L]], recursive = FALSE), fixed = TRUE)
  powers <- lapply(strsplit(powers, ","), as.integer)
  list(
    "coeffs" = coeffs, 
    "powers" = powers
  )
}

stringToPol("-x^(0, 2, 0) + 3/2 x^(1, 0) - 4")