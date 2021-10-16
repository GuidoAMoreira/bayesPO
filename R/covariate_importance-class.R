#' Class for covariates importance matrices
#'
#' Objects of this class is the output of the "covariates_importance" object
#' from the \code{\link{bayesPO_fit-class}}. It can be plotted which uses
#' the \code{\link[graphics]{graphics}} package. The \code{print} method
#' gives a point-wise estimation, the same seen in the \code{bacplot} method.
#' Both \code{plot} and \code{boxplot} methods use the posterior distribution
#' of the importance.
#'
#' @details
#' Objects of this class have two matrices where the Monte Carlo samples on the
#' rows and parameters on the columns. One matrix is for the intensity
#' importance and the other for the observability importance.
#' @name covariates_importance-class
NULL

#' @rdname covariates_importance-class
#' @param x The \code{covariates_importance} object.
#' @param component Either \code{"intensity"}, \code{"observability"} or
#' \code{"both"}.
#' @param ... Ignored.
#' @return The invisible object.
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

#' @rdname covariates_importance-class
#' @param x The \code{covariates_importance} object.
#' @param component Either \code{"intensity"}, \code{"observability"} or
#' \code{"both"}.
#' @param y Either \code{"interval"} or \code{"density"}. The formal gives
#' vertical credible intervals, and the latter gives separate density plots
#' with the specified quantiles as vertical lines.
#' @param quantiles A 2- or 3-simensional vector with the desired quantiles
#' specified. If 3-dimensiona, the middle point is drawn as a dot when the
#' \code{y} parameter is set as \code{"interval"}.
#' @param ... Other parameters passed to \code{\link[graphics]{plot}}.
#' @return Nothing is returned. Plot is called and drawn on the configured
#' device.
#' @method plot covariates_importance
#' @importFrom stats density
#' @importFrom graphics plot axis segments points par
#' @importFrom methods hasArg
#' @importFrom tools toTitleCase
#' @importFrom stats quantile
#' @export
plot.covariates_importance <- function(x, component = "intensity",
                                       y = "importance",
                                       quantiles = c(0.025, 0.5, 0.975), ...) {
  component <- tolower(component)
  y <- tolower(y)
  stopifnot(component %in% c("intensity", "observability", "both"),
            y %in% c("importance", "density"),
            length(quantiles) %in% c(2, 3), all(quantiles > 0 & quantiles < 1))

  quantiles <- sort(quantiles)
  if (component != "both") {
    intervals <- apply(x[[component]], 2, stats::quantile,
                       quantiles[c(1, length(quantiles))])
    plots <- ncol(x[[component]])
    if (!methods::hasArg("ylim")) yl <- c(min(intervals) * 0.9,
                                        max(intervals) * 1.1) else yl <- get("ylim")
    if (y == "importance") {
      graphics::plot(NA, xlim = c(1, plots), xaxt = "n",
                     main = paste(tools::toTitleCase(component), "variables"),
           ylim = yl,
           xlab = ifelse(!methods::hasArg("xlab"), "Covariates", get("xlab")),
           ylab = ifelse(!methods::hasArg("ylab"), "Importance", get("ylab")), ...)
      graphics::axis(1, at = 1:plots, colnames(x[[component]]))
      graphics::segments(1:plots, intervals[1, ], 1:plots, intervals[2, ])
      if (length(quantiles) == 3)
        graphics::points(1:plots, apply(x[[component]], 2, stats::quantile,
                                        quantiles[2]), pch = 19)
    } else {
      large_row_nmbr <- floor(sqrt(plots))
      col_nmbrs <- large_row_nmbr:(large_row_nmbr + 2)
      pp <- graphics::par(mfrow = c(large_row_nmbr,
                          col_nmbrs[min(which(col_nmbrs * large_row_nmbr >= plots))]))
      on.exit(graphics::par(pp))
      for (p in 1:plots){
        graphics::plot(stats::density(x[[component]][, p]),
                       xlab = paste(colnames(x[[component]])[p], "importance"),
                       main = colnames(x[[component]])[p],...)
      }
    }
  } else {
    plot.covariates_importance(x, "intensity", y)
    cat("\nHit <Return> to see next plot: ")
    line <- readline()
    plot.covariates_importance(x, "observability", y)
  }
}

#' @rdname covariates_importance-class
#' @param height The \code{covariates_importance} object.
#' @param component Either \code{"intensity"}, \code{"observability"} or
#' \code{"both"}.
#' @param ... Other parameters passed to \code{\link[graphics]{barplot}}.
#' @method barplot covariates_importance
#' @return A barplot. See \code{barplot} for details. If component is selected
#' as \code{"both"}, only the second barplot is returned.
#' @seealso \code{\link[graphics]{barplot}}.
#' @importFrom graphics barplot
#' @export
barplot.covariates_importance <- function(height, component = "intensity", y, ...) {
  component <- tolower(component)
  stopifnot(component %in% c("intensity", "observability", "both"))
  if (component != "both") {
    graphics::barplot(colMeans(height[[component]]), ...)
  } else {
    barplot.covariates_importance(height, "intensity")
    cat("\nHit <Return> to see next plot: ")
    line <- readline()
    barplot.covariates_importance(height, "observability")
  }
}

#' @rdname covariates_importance-class
#' @param x The \code{covariates_importance} object.
#' @param component Either \code{"intensity"}, \code{"observability"} or
#' \code{"both"}.
#' @param ... Other parameters passed to \code{\link[graphics]{boxplot}}.
#' @method boxplot covariates_importance
#' @return A boxplot. See \code{boxplot} for details. If component is selected
#' as \code{"both"}, only the second boxplot is returned.
#' @seealso \code{\link[graphics]{boxplot}}.
#' @importFrom graphics boxplot
#' @export
boxplot.covariates_importance <- function(x, component = "intensity", ...) {
  component <- tolower(component)
  stopifnot(component %in% c("intensity", "observability", "both"))
  if (component != "both") {
    graphics::boxplot(x[[component]], ...)
  } else {
    boxplot.covariates_importance(x, "intensity")
    cat("\nHit <Return> to see next plot: ")
    line <- readline()
    boxplot.covariates_importance(x, "observability")
  }
}
