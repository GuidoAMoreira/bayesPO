#' Class for covariates importance matrices
#'
#' Objects of this class is the output of the "covariates_importance" object
#' from the \code{\link{bayesPO_fit-class}}. It can be plotted which uses
#' the \code{\link[graphics]{graphics}} package.
#'
#' @details
#' Objects of this class have two matrices. One for the intensity importance
#' and one for the observability importance.
#' @name covariates_importance-class
NULL

#' @name covariates_importance-class
#' @param x The \code{covariates_importance} object.
#' @param component Either \code{"intensity"}, \code{"observability"} or
#' \code{"both"}.
#' @param ... Ignored.
#' @method print covariates_importance
#' @export
print.covariates_importance <- function(x, component = "intensity", ...){
  component <- tolower(component)
  stopifnot(component %in% c("intensity", "observability", "both"))
  if (component != "both") print(colMeans(x[[component]])) else {
    cat("Intensity covariates importance:\n\n")
    print.covariates_importance(x, "intensity")
    cat("\n------------------------\n\nObservability covariates importance:\n\n")
    print.covariates_importance(x, "observability")
  }

  invisible(x)
}

#' @name covariates_importance-class
#' @param x The \code{covariates_importance} object.
#' @param component Either \code{"intensity"}, \code{"observability"} or
#' \code{"both"}.
#' @param y Ignored.
#' @param ... Other parameters passed to \code{\link[graphics]{plot}}.
#' @method plot covariates_importance
#' @importFrom stats density
#' @importFrom graphics plot
#' @export
plot.covariates_importance <- function(x, component = "intensity", y, ...) {
  component <- tolower(component)
  stopifnot(component %in% c("intensity", "observability", "both"))
  if (component != "both") {
    plots <- ncol(x[[component]])
    large_row_nmbr <- floor(sqrt(plots))
    col_nmbrs <- large_row_nmbr:(large_row_nmbr + 2)
    pp <- par(mfrow = c(large_row_nmbr,
                        col_nmbrs[min(which(col_nmbrs * large_row_nmbr >= plots))]))
    for (p in 1:plots){
      graphics::plot(stats::density(x[[component]][, p]),
                     xlab = paste(colnames(x[[component]])[p], "importance"),
                     main = colnames(x[[component]])[p],...)
    }
    par(pp)
  } else {
    plot.covariates_importance(x, "intensity")
    cat("\nHit <Return> to see next plot: ")
    line <- readline()
    plot.covariates_importance(x, "observability")
  }
}

#' @name covariates_importance-class
#' @param height The \code{covariates_importance} object.
#' @param component Either \code{"intensity"}, \code{"observability"} or
#' \code{"both"}.
#' @param ... Other parameters passed to \code{\link[graphics]{barplot}}.
#' @method barplot covariates_importance
#' @importFrom graphics barplot
#' @export
barplot.covariates_importance <- function(x, component = "intensity", y, ...) {
  component <- tolower(component)
  stopifnot(component %in% c("intensity", "observability", "both"))
  if (component != "both") {
    graphics::barplot(colMeans(x[[component]]), ...)
  } else {
    plot.covariates_importance(x, "intensity")
    cat("\nHit <Return> to see next plot: ")
    line <- readline()
    plot.covariates_importance(x, "observability")
  }
}
