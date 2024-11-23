#' Simulated Presence-Only data
#'
#' This small dataset was simulated using the package's model, with the probit
#' link in both intensity and observability. It is used to showcase the model
#' in the vignette.
#'
#' @format ## `bayesPO_sim`
#' This object contains all data necessary to fit a model and also to plot the
#' data. It is a list with components:
#' \describe{
#' \item{po}{A matrix with 420 rows and 2 columns. The first column is the
#' intensity covariates for the sightings, the second is the observability one.}
#' \item{bkg}{A matrix with 2500 rows and 2 columns. The first column is the
#' intensity covariates for the background, the second is the observability one.}
#' \item{grid}{The grid coordinates for the background covariates.}
#' \item{po_sightings}{A 605 long logical vector that identified the
#' occurrences coordinates about whether they are observed or not.}
#' \item{occurrences_points}{A matrix with the coordinates where all
#' occurrences happened (whether they were observed or not.)}
#' }
#' @source The data was simulated.
"bayesPO_sim"
