#' @include initial-class.R prior-class.R model-class.R
NULL

#' Fit Presence-Only data using a Bayesian Poisson Process model
#'
#' The model uses a data augmentation scheme to avoid performing approximations
#' on the likelihood function.
#' @param data A list containing a matrix with the covariates for the observed
#' points, a matrix with the covariates for the background, the region's area
#' and the column choices for the intensity and observability sets. See details.
#' @param intensityLink A string with the chosen link function for the
#' intensity process. Currently only \code{"logit"} is accepted.
#' @param intensityPrior A string with the chosen prior distribution for the
#' intensity effects \code{beta}. Currently only \code{"normal"} is accepted.
#' @param observabilityLink A string with the chosen link function for the
#' observability process. Currently only \code{"logit"} is accepted.
#' @param observabilityPrior A string with the chosen prior distribution for the
#' observability effects \code{delta}. Currently only \code{"normal"} is accepted.
#' @param init The initial values for the chains. The string \code{"random"} can
#' be used to use randomly chosen values. See details for instructions on how
#' to provide specific values.
#' @param prior Parameters for the prior. It must be a \code{list} containing
#' a list for each of the parameters \code{beta}, \code{delta} and
#' \code{lambdaStar}. See details for more information.
#' @param mcmc_setup MCMC configuration in the form of a list containing
#' elements burnin, thin and n_iter. See details for more information.
#' @param chains Number of independent chains to run. If more than 1 are
#' requested and argument \code{init} is not \code{"random"}, initial values
#' must be provided for all chains. See details for more instructions about
#' this.
#' @return An object of class \code{"bayesPO_fit"}. It has its own print and
#' summary methods. It has methods for as.array, as.matrix and as.data.frame,
#' thus can be used with the bayesplot package (as.array) or ggplot2
#' (as.data.frame).
#' @details
#' ## Input instructions
#' A few notes about properly feeding information to the model are provided
#' here.
#' ### data
#' The data input must be a list with the following elements:
#' \itemize{
#' \item{observed}{A matrix containing the covariate values for the
#' presence-only observed points.}
#' \item{background}{A matrix containing the covariate values for the
#' background. The values must span over the entire studied region.}
#' \item{area}{A single positive value informing the studied region's area.
#' Changing this value will not (greatly) affect the estimation of \code{beta}
#' and \code{delta}, but it will affect \code{lambdaStar} and predictions using
#' the model. This means that a correct value must be provided for adequate
#' predictions.}
#' \item{intensity}{A vector indicating which columns explain the species'
#' intensity, that is, its occurrence. It can be either a numeric vector or
#' a character vector. If numeric, it indicates which column positions are
#' to be included as intensity covariates. Arguments \code{observed} and
#' \code{background} must have identical columns positions for the covariates.
#' If character, it indicates the names of the columns to be included as
#' intensity covariates. Both \code{observed} and \code{background} must have
#' the indicated columns, but they do not need to be in the same position.
#' Parameters \code{beta} measure their effects in the model.}
#' \item{observability}{A vector indicating which columns explain the species'
#' observability, that is, the probability that its occurrences are observed.
#' It can be either a numeric vector or a character vector. If numeric, it
#' indicates which column positions are to be included as observability
#' covariates. Arguments \code{observed} and \code{background} must have
#' identical columns positions for the covariates. If character, it indicates
#' the names of the columns to be included as observability covariates.
#' Both \code{observed} and \code{background} must have the indicated columns,
#' but they do not need to be in the same position. Parameters \code{delta}
#' measure their effects in the model.}
#' }
#'
#' ### init
#' The initial values for the MCMC chains can be chosen to be random with the
#' \code{"random"} string. However, if the user wishes to input specific values,
#' they must be included in a list inside a list. The outer list is to group
#' up the initial values for the different chains.
#'
#' The inner ones inform the model the initial values for each chain, so the
#' first list includes the values for the first chain, and so forth. If any of
#' the inner lists does not contain one of them, a random value is included.
#'
#' For \code{beta} and \code{delta}, a single value can be given, at which case
#' the same value is used for the entire vector. Otherwise, the size of
#' \code{beta} and \code{delta} must match the selections of one plus intensity
#' and observability covariates, respectively. The extra value is meant for the
#' intercepts, which go first in the order, meaning that the first initial
#' value of each vector is associated with the intercepts.
#'
#' For \code{lambdaStar} only a single value per chain must be provided.
#'
#' ### prior
#' Different prior choices imply different parameter choices.
#'
#' For the \code{beta} and \code{delta} "normal" choice, a mean vector and
#' covariance matrix must be provided. The \code{beta} mean vector must be one
#' plus the choice size of intensity covariates and the \code{delta} mean
#' vector must be one plus the choice size of observability covariates. The
#' extra one is for the intercept, which comes first, meaning that the first
#' value of the vector is the prior mean of the intercept. The covariance
#' matrices must be symmetric and positive definite.
#'
#' Alternatively, a single value is provided. For the mean, this means that
#' the prior mean is equal to the chosen value for all \code{beta} parameters.
#' The logic is analogous for \code{delta.} If the covariance is a single
#' positive number, then the respective matrix is an identity multiplied by
#' this number.
#'
#' For \code{lambdaStar} the prior is a gamma distribution with parameters a and b
#' for the shape and rate parameters, respectively.
#'
#' ### mcmc_setup
#'
#'
#'
#' @importFrom coda mcmc mcmc.list
#' @export
methods::setGeneric("fit_bayesPO",function(object,background,area = 1,mcmc_setup = list(n_iter = 5000)){standardGeneric("fit_bayesPO")})

#' The fit_bayesPO method for the bayesPO_model class.
#'
#' @param object A bayesPO_model object.
#' @param background A matrix with the covariates values in the region background.
#' @param area A positive number with the region's area
#' @param mcmc_setup A list with the components \code{burnin}, \code{thin} and
#' \code{n_iter} where only the latter is mandatory.
methods::setMethod("fit_bayesPO",signature(object="bayesPO_model",background="matrix"),
          function(object,background,area = 1,mcmc_setup = list(n_iter = 5000)){
  ## Verifying background names if columns are selected by column name. Crewating background selection variables
  backConfig <- checkFormatBackground(object,background)

  ## Verifying area
  if (length(area) > 1) stop("Argument area must have length 1.")
  if (area <= 0) stop("Argument area must be positive.")
  ## Verifying mcmc_setup
  if (is.null(mcmc_setup$burnin)) mcmc_setup$burnin <- 0
  if (is.null(mcmc_setup$thin)) mcmc_setup$thin <- 1
  mcmc_setup <- checkMCMCSetup(mcmc_setup)

  # Helper parameters
  nb <- length(methods::slot(object,"intensitySelection"))+1
  nd <- length(methods::slot(object,"observabilitySelection"))+1
  npar <- nb+nd+4
  chains <- length(methods::slot(object,"init"))

  # Set up parameter names to construct the coda::mcmc objects
  if (length(methods::slot(object,"iSelectedColumns")>0)) betanames <- c("(Intensity intercept)",methods::slot(object,"iSelectedColumns")) else betanames = paste0("beta_",1:nb-1)
  if (length(methods::slot(object,"oSelectedColumns")>0)) deltanames <- c("(Observability intercept)",methods::slot(object,"oSelectedColumns")) else deltanames = paste0("delta_",1:nd-1)
  parnames <- c(betanames,deltanames,
               "lambdaStar","n_U","n_Xprime","log_Posterior"
  )

  time <- Sys.time()
  mcmcRun <- list()
  for (c in 1:chains){
    if (chains > 1) cat("Starting chain ",c,".\n",sep="")
    mcmcRun[[c]] <- do.call(cbind,
                           runBayesPO(methods::slot(methods::slot(object,"init")[[c]],"beta"),
                                      methods::slot(methods::slot(object,"init")[[c]],"delta"),
                                      methods::slot(methods::slot(object,"init")[[c]],"lambdaStar"),
                                      paste0(methods::slot(object,"intensityLink"),"_",
                                             methods::slot(methods::slot(methods::slot(object,"prior"),"beta"),"family")),
                                      paste0(methods::slot(object,"observabilityLink"),"_",
                                             methods::slot(methods::slot(methods::slot(object,"prior"),"delta"),"family")),
                                      methods::slot(methods::slot(methods::slot(object,"prior"),"lambdaStar"),"family"),
                                      retrievePars(methods::slot(methods::slot(object,"prior"),"beta")),
                                      retrievePars(methods::slot(methods::slot(object,"prior"),"delta")),
                                      retrievePars(methods::slot(methods::slot(object,"prior"),"lambdaStar")),
                                      ifelse(is.integer(background),"int_mat","num_mat"), background, area,
                                      ifelse(is.integer(methods::slot(object,"po")),"int_mat","num_mat"),methods::slot(object,"po"),
                                      backConfig$bIS-1,backConfig$bOS-1,
                                      methods::slot(object,"intensitySelection")-1,
                                      methods::slot(object,"observabilitySelection")-1,
                                      mcmc_setup$burnin, mcmc_setup$thin, mcmc_setup$n_iter)
    )
    colnames(mcmcRun[[c]]) <- parnames
    mcmcRun[[c]] <- coda::mcmc(mcmcRun[[c]],thin = mcmc_setup$thin)
    if (chains > 1) cat("Finished chain ",c,".\n\n",sep="")
  }
  if (chains > 1) cat("Total computation time:",format(unclass(Sys.time()-time), digits = 2),attr(Sys.time()-time,"units"),".\n")

  return(methods::new("bayesPO_fit",
         fit = do.call(coda::mcmc.list,mcmcRun),
         original = object,
         backgroundSummary = summary(background),
         area = area,
         parnames = parnames,
         mcmc_setup = mcmc_setup))
})

checkFormatBackground <- function(object,background){
  # Intensity
  if (length(methods::slot(object,"iSelectedColumns")) > 0){
    backIntensitySelection <- c()
    for (col in methods::slot(object,"iSelectedColumns")){
      if (col %in% colnames(background))
        backIntensitySelection <- c(backIntensitySelection, which(col == colnames(background)))
      else stop(paste0("Column ",col," not found in the background covariates"))
    }
  } else
    if (max(methods::slot(object,"intensitySelection")) > ncol(background))
      stop(paste0("Requested background column ",max(methods::slot(object,"intensitySelection")),
                  " not available in background which only has ",ncol(background)," columns."))
    else backIntensitySelection = methods::slot(object,"intensitySelection")

  # Observability
  if (length(methods::slot(object,"oSelectedColumns")) > 0){
    backObservabilitySelection <- c()
    for (col in methods::slot(object,"oSelectedColumns")){
      if (col %in% colnames(background))
        backObservabilitySelection <- c(backObservabilitySelection, which(col == colnames(background)))
      else stop(paste0("Column ",col," not found in the background covariates"))
    }
  } else
    if (max(methods::slot(object,"observabilitySelection")) > ncol(background))
      stop(paste0("Requested background column ",max(methods::slot(object,"observabilitySelection")),
                  " not available in background which only has ",ncol(background)," columns."))
  else
    backObservabilitySelection <- methods::slot(object,"observabilitySelection")

  return(list(bIS = backIntensitySelection, bOS = backObservabilitySelection))
}

checkMCMCSetup <- function(setup){
  if (!("n_iter" %in% names(setup))) stop("MCMC setup must contain n_iter.")
  setup$burnin <- as.numeric(setup$burnin)
  if (is.na(setup$burnin)) stop("Could not turn MCMC configuration parameter burnin to number.")
  if (length(setup$burnin) > 1) stop("MCMC configuration parameter burnin must have length 1.")
  if (setup$burnin != floor(setup$burnin)) stop("MCMC configuration parameter burnin must be a posivite integer")
  if (setup$burnin < 0) stop("MCMC configuration parameter burnin must be a non-negative integer.")
  setup$thin <- as.numeric(setup$thin)
  if (is.na(setup$thin)) stop("Could not turn MCMC configuration parameter thin to number.")
  if (length(setup$thin)>1) stop("MCMC configuration parameter thin must have length 1.")
  if (setup$thin != floor(setup$thin)) stop("MCMC configuration parameter thin must be a posivite integer")
  if (setup$thin <= 0) stop("MCMC configuration parameter thin must be a posivite integer.")
  setup$n_iter = as.numeric(setup$n_iter)
  if (is.na(setup$n_iter)) stop("Could not turn MCMC configuration parameter n_iter to number.")
  if (length(setup$n_iter)>1) stop("MCMC configuration parameter n_iter must have length 1.")
  if (setup$n_iter != floor(setup$n_iter)) stop("MCMC configuration parameter n_iter must be a posivite integer.")
  if (setup$n_iter <= 0) stop("MCMC configuration parameter n_iter must be a posivite integer.")
  if (setup$n_iter<setup$thin) stop("MCMC configuration parameter thin is too large. It cannot be larger than n_iter. No MCMC performed.")

  setup
}


#' @export
bayesPO_model = function(po,intensitySelection,
                         observabilitySelection,
                         intensityLink = "logit",observabilityLink = "logit",
                         initial_values = initial(length(intensitySelection)+1,
                                                  length(observabilitySelection)+1,
                                                  nrow(po),random=TRUE),
                         joint_prior = prior(
                           beta = NormalPrior(
                             rep(0,length(intensitySelection) + 1),
                             10*diag(length(intensitySelection) + 1)),
                           delta = NormalPrior(
                             rep(0,length(observabilitySelection) + 1),
                             10*diag(length(observabilitySelection) + 1)),
                           lambdaStar = GammaPrior(
                             1e-10, 1e-10
                           ))){
  #if (missing(po)) {return(model_builder())} # Easy and interactive model building

  if (is.numeric(intensitySelection)) icharSel = character() else
    if (is.character(intensitySelection)){ # Create columns positions for intensity
      icharSel = intensitySelection
      intensitySelection = c()
      for (selected in icharSel){
        if (!(selected %in% colnames(po))){
          message(paste0("Column ",selected," not found in the observed points covariates. Ignoring it."))
          next
        }
        intensitySelection = c(intensitySelection,which(selected==colnames(po)))
      }
    }
  if (is.numeric(observabilitySelection)) ocharSel = character() else
    if (is.character(observabilitySelection)){ # Create columns positions for intensity
      ocharSel = observabilitySelection
      observabilitySelection = c()
      for (selected in ocharSel){ # Putting indexes together
        if (!(selected %in% colnames(po))){
          message(paste0("Column ",selected," not found in the observed points covariates. Ignoring it."))
          next
        }
        observabilitySelection = c(observabilitySelection,which(selected==colnames(po)))
      }
    }
  if (is(initial_values,"bayesPO_initial"))
    initial_values = list(initial_values)
  if (is.numeric(initial_values))
    initial_values = initial(length(intensitySelection)+1,
                             length(observabilitySelection)+1,
                             nrow(po),random=TRUE) * initial_values
  return(methods::new("bayesPO_model",po=po,intensityLink=intensityLink,intensitySelection=intensitySelection,
               observabilityLink=observabilityLink,observabilitySelection=observabilitySelection,
               init=initial_values,prior=joint_prior,iSelectedColumns=icharSel,oSelectedColumns=ocharSel))
}
