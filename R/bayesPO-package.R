## usethis namespace: start
#' @useDynLib bayesPO, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppProgress
## usethis namespace: end
#' @importFrom methods show
#' @importFrom methods initialize
#' @importFrom methods is

.onUnload <- function (libpath) {
  library.dynam.unload("bayesPO", libpath)
}
