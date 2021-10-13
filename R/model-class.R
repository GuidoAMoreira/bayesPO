#' @include initial-class.R prior-class.R
NULL

#' Class that defines a model for the bayesPO package.
#'
#' The model includes the presence-only data, all selected variables, the link
#' functions for \eqn{q} and \eqn{p}, the initial values and the prior
#' distribution.
#' @field po The matrix containing the covariates values for the data.
#' @field intensityLink A string informing about the chosen link for the
#' intensity covariates. Current acceptable choice is only \code{"logit"}.
#' @field intensitySelection A vector containing the indexes of the selected
#' intensity columns in the \code{po} matrix.
#' @field observabilityLink A string informing about the chosen link for the
#' observability covariates. Current acceptable choice is only \code{"logit"}.
#' @field observabilitySelection A vector containing the indexes of the selected
#' observability columns in the \code{po} matrix.
#' @field init A list with objects of class \code{bayesPO_initial} indicating
#' the initial values for each chain. The length of this list tells the program
#' how many chains are requested to be run.
#' @field prior An object of class \code{bayesPO_prior} which indicates the
#' joint prior distribution for the model parameters.
#' @field iSelectedColumns If the intensity covariates selection was made with
#' the name of the columns, they are stored in this slot.
#' @field iSelectedColumns If the observability covariates selection was made
#' with the name of the columns, they are stored in this slot.
#' @seealso \code{\link{bayesPO_initial-class}} and
#' \code{\link{bayesPO_prior-class}} and \code{\link{bayesPO_model}}
#' @export
#' @exportClass bayesPO_model
methods::setClass("bayesPO_model",
                  methods::representation(po = "matrix",
                        intensityLink = "character",
                        intensitySelection = "numeric",
                        observabilityLink = "character",
                        observabilitySelection = "numeric",
                        init = "list",
                        prior = "bayesPO_prior",
                        iSelectedColumns = "character",
                        oSelectedColumns = "character"
                        ),
         validity = function(object){
           s <- function(n) methods::slot(object, n)
           if (length(s("intensityLink"))>1)
             stop("Argument intensityLink must have length 1.")
           if (length(s("observabilityLink"))>1)
             stop("Argument observabilityLink must have length 1.")
           validLinks = c("logit")
           if (!(s("intensityLink") %in% validLinks))
             stop(paste0(s("intensityLink")," is not a valid link. Accepted options for current version are\n",validLinks))
           if (!(s("observabilityLink") %in% validLinks)) stop(paste0(s("observabilityLink")," is not a valid link. Accepted options for current version are\n", validLinks))
           if (any(s("intensitySelection") != as.integer(s("intensitySelection")))) stop("intensitySelection must be integers.")
           if (any(s("observabilitySelection") != as.integer(s("observabilitySelection")))) stop("observabilitySelection must be integers.")
           nbSel = length(s("intensitySelection")) + 1; ndSel = length(s("observabilitySelection")) + 1
           for (i in s("init")){
             if (!is(i,"bayesPO_initial")) stop("Initial values must be constructed with the initial function.")
             if (length(methods::slot(i,"beta")) != nbSel) stop(paste("\nInitial values for beta has the wrong size. Expected size:", nbSel))
             if (length(methods::slot(i,"delta")) != ndSel) stop(paste("\nInitial values for delta has the wrong size. Expected size:", ndSel))
           }
           if (length(methods::slot(methods::slot(s("prior"),"beta"),"mu")) != nbSel) stop(paste("Prior for beta has wrong sized parameters. Expected size:", nbSel))
           if (length(methods::slot(methods::slot(s("prior"),"delta"),"mu")) != ndSel) stop(paste("Prior for delta has wrong sized parameters. Expected size:", ndSel))
           TRUE
         })

#' Initialize method for bayesPO_model objects
#'
#' Fills the object with the necessary slots. Provides some text informing
#' the loaded data.
#' @param .Object An empty bayesPO_model object.
#' @param po The matrix containing the covariates in the observed locations.
#' @param intensityLink A string containing the intensity Link function.
#' @param intensitySelection The positions of selected columns for the intensity
#' covariates. Can be empty if the selection was from a character vector.
#' @param observabilityLink A string containing the observability Link function.
#' @param observabilitySelection The positions of selected columns for the
#' observability covariates. Can be empty if the selection was from a character
#' vector.
#' @param init A bayesPO_initial object.
#' @param prior A bayesPO_prior object.
#' @param iSelectedColumns A character vector containing the selected
#' intensity covariates. Can be empty if the selection was from the columns
#' positions.
#' @param oSelectedColumns A character vector containing the selected
#' observability covariates. Can be empty if the selection was from the columns
#' positions.
#' @export
#' @exportMethod initialize
methods::setMethod("initialize", "bayesPO_model", function(.Object, po, intensityLink, intensitySelection,
                                                observabilityLink, observabilitySelection,
                                                init, prior, iSelectedColumns, oSelectedColumns){
  cat("Loading data with", nrow(po), "observed points.\n")
  methods::slot(.Object, "po") <- po
  methods::slot(.Object, "intensityLink") <- intensityLink; methods::slot(.Object, "observabilityLink") = observabilityLink
  methods::slot(.Object, "intensitySelection") <- intensitySelection; methods::slot(.Object, "observabilitySelection") = observabilitySelection
  methods::slot(.Object, "init") <- init
  methods::slot(.Object, "prior") <- prior
  methods::slot(.Object, "iSelectedColumns") <- iSelectedColumns
  methods::slot(.Object, "oSelectedColumns") <- oSelectedColumns
  methods::validObject(.Object)
  cat("Data loaded successfully with ", length(intensitySelection), " intensity variables and ",
      length(observabilitySelection), " observability variables selected.\n", length(init), " chains ",
      ifelse(length(init) > 1, "were","was"), " initialized.\n", sep = "")
  if (!length(iSelectedColumns)) cat("Intensity covariates selected with column indexes. Make sure the background covariates are in the same position.\n")
  if (!length(oSelectedColumns)) cat("Observability covariates selected with column indexes. Make sure the background covariates are in the same position.\n")
  .Object
})

#' @rdname bayesPO_model-class
#'
#' @param x The bayesPO_model object.
#' @export
#' @exportMethod names
methods::setMethod("names", "bayesPO_model",
          function(x) c("po", "intensityLink", "intensitySelection",
                        "observabilityLink", "observabilitySelection",
                        "initial values", "prior"))

#' @rdname bayesPO_model-class
#'
#' @param x The bayesPO_model object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$","bayesPO_model",function(x,name){
  if (name %in% c("initial values", "initial")) methods::slot(x, "init") else methods::slot(x, name)
})

#' @rdname bayesPO_model-class
#'
#' @param x The bayesPO_model object.
#' @param name The requested slot.
#' @param value New value.
#' @export
#' @exportMethod $<-
methods::setMethod("$<-","bayesPO_model",function(x, name, value){
  if (name %in% c("initial values", "initial")) methods::slot(x, "init") = value else methods::slot(x, name) = value
  methods::validObject(x)
  x
})


