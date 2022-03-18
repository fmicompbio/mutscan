## usethis namespace: start
#' @useDynLib mutscan, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

# global variables definition for R CMD check
# (avoid NOTE about "no visible binding for global variable")
#' @importFrom utils globalVariables
utils::globalVariables(c("."))
