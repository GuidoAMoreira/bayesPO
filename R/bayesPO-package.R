## usethis namespace: start
#' @useDynLib bayesPO, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppProgress
## usethis namespace: end


.onUnload <- function (libpath) {
  library.dynam.unload("bayesPO", libpath)
}
