#' @include initial-class.R prior-class.R model-class.R fit-class.R
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
#' elements burnin, thin and iter. See details for more information.
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
methods::setGeneric("fit_bayesPO", function(object, background,
                                            mcmc_setup = list(iter = 5000), ...)
  standardGeneric("fit_bayesPO"))

#' The fit_bayesPO method for the bayesPO_model class.
#'
#' @param object A bayesPO_model object.
#' @param background A matrix with the covariates values in the region background.
#' @param area A positive number with the region's area
#' @param mcmc_setup A list with the components \code{burnin}, \code{thin} and
#' \code{iter} where only the latter is mandatory.
#' @name fit_bayesPO
#' @importFrom parallel detectCores
methods::setMethod("fit_bayesPO", signature(object="bayesPO_model",
                                            background="matrix"),
                   function(object, background, mcmc_setup, area = 1, cores = 1){
                     ## Verifying background names if columns are selected by column name. Crewating background selection variables
                     backConfig <- checkFormatBackground(object, background)
                     cores <- 1

                     # Helper function
                     s <- function(n) methods::slot(object, n)

                     ## Verifying area
                     stopifnot(length(area) == 1, area > 0)
                     ## Verifying mcmc_setup - begin
                     if (is.null(mcmc_setup$burnin)) mcmc_setup$burnin <- 0
                     if (is.null(mcmc_setup$thin)) mcmc_setup$thin <- 1
                     mcmc_setup$burnin <- as.numeric(mcmc_setup$burnin)
                     mcmc_setup$thin <- as.numeric(mcmc_setup$thin)
                     mcmc_setup$iter <- as.numeric(mcmc_setup$iter)
                     stopifnot("iter" %in% names(mcmc_setup), !is.na(mcmc_setup$burnin),
                               length(mcmc_setup$burnin) == 1, mcmc_setup$burnin == floor(mcmc_setup$burnin),
                               mcmc_setup$burnin >= 0, !is.na(mcmc_setup$thin), length(mcmc_setup$thin) == 1,
                               mcmc_setup$thin == floor(mcmc_setup$thin), mcmc_setup$thin > 0,
                               !is.na(mcmc_setup$iter), length(mcmc_setup$iter) == 1,
                               mcmc_setup$iter == floor(mcmc_setup$iter), mcmc_setup$iter > 0,
                               mcmc_setup$iter >= mcmc_setup$thin,
                               cores > 0, cores == floor(cores), length(cores) == 1,
                               cores < parallel::detectCores())
                     ## Verifying mcmc_setup - end

                     # Helper parameters
                     nb <- length(s("intensitySelection")) + 1
                     nd <- length(s("observabilitySelection")) + 1
                     npar <- nb + nd + 4
                     chains <- length(s("init"))

                     # Set up parameter names to construct the coda::mcmc objects
                     if (length(s("iSelectedColumns") > 0))
                       betanames <- c("(Intensity intercept)",
                                      s("iSelectedColumns"))
                     else betanames = paste0("beta_", 1:nb - 1)
                     if (length(s("oSelectedColumns") > 0))
                       deltanames <- c("(Observability intercept)",
                                       s("oSelectedColumns"))
                     else deltanames = paste0("delta_", 1:nd - 1)
                     parnames <- c(betanames,deltanames,
                                   "lambdaStar", "n_U", "n_Xprime", "log_Posterior"
                     )

                     time <- Sys.time()
                     mcmcRun <- list()
                     for (c in 1:chains){
                       if (chains > 1) cat("Starting chain ",c,".\n",sep="")
                       mcmcRun[[c]] <- do.call(cbind,
                                               runBayesPO(
                                                 methods::slot(s("init")[[c]],"beta"),
                                                 methods::slot(s("init")[[c]],"delta"), # Init delta
                                                 methods::slot(s("init")[[c]],"lambdaStar"), # Init lambdaStar
                                                 paste0(s("intensityLink"), "_", # intenisty Link + prior
                                                        methods::slot(
                                                          methods::slot(s("prior"), "beta"), "family")),
                                                 paste0(s("observabilityLink"),"_", # observability Link + prior
                                                        methods::slot(
                                                          methods::slot(s("prior"), "delta"), "family")),
                                                 methods::slot( # lambdaStar prior
                                                   methods::slot(
                                                     s("prior"),"lambdaStar"), "family"),
                                                 retrievePars(methods::slot( # beta prior parameters
                                                   s("prior"),"beta")),
                                                 retrievePars(methods::slot( # delta prior parameters
                                                   s("prior"),"delta")),
                                                 retrievePars( # lambdaStar prior parameters
                                                   methods::slot(
                                                     s("prior"),"lambdaStar")),
                                                 ifelse(is.integer(background), # background class
                                                        "int_mat", "num_mat"),
                                                 background, # background data
                                                 area, # region area
                                                 ifelse(is.integer(s("po")),
                                                        "int_mat", "num_mat"), # po data class
                                                 s("po"), # po data
                                                 backConfig$bIS - 1, # background intensity columns
                                                 backConfig$bOS - 1, # background observability colmns
                                                 s("intensitySelection") - 1, # po intensity columns
                                                 s("observabilitySelection") - 1, # po obserability columns
                                                 mcmc_setup$burnin, # MCMC burn-in
                                                 mcmc_setup$thin, # MCMC thin
                                                 mcmc_setup$iter, # MCMC iterations
                                                 cores)
                       )
                       colnames(mcmcRun[[c]]) <- parnames
                       mcmcRun[[c]] <- coda::mcmc(mcmcRun[[c]], thin = mcmc_setup$thin)
                       if (chains > 1) cat("Finished chain ",c,".\n\n",sep="")
                     }
                     if (chains > 1) cat("Total computation time: ", format(unclass(Sys.time()-time),
                                                                            digits = 2), " ",
                                         attr(Sys.time() - time, "units"), ".\n", sep="")

                     return(methods::new("bayesPO_fit",
                                         fit = do.call(coda::mcmc.list, mcmcRun),
                                         original = object,
                                         backgroundSummary = summary(background),
                                         area = area,
                                         parnames = parnames,
                                         mcmc_setup = mcmc_setup))
                   })

methods::setMethod("fit_bayesPO", signature(object = "bayesPO_fit",
                                            background = "matrix"),
                   function(object, background,
                            mcmc_setup = list(iter = object$mcmc_setup$iter), cores = 1){
                     # Helper function
                     s <- function(n) methods::slot(object, n)
                     so <- function(n) methods::slot(s("original"), n)
                     backConfig <- checkFormatBackground(s("original"), background)

                     # Check background differences
                     cat("Performing error check...\n")
                     stopifnot(all.equal(s("backgroundSummary"), summary(background)),
                               "iter" %in% names(mcmc_setup), !is.na(mcmc_setup$iter),
                               length(mcmc_setup$iter) == 1,
                               mcmc_setup$iter == floor(mcmc_setup$iter), mcmc_setup$iter > 0,
                               mcmc_setup$iter >= mcmc_setup$thin,
                               cores > 0, cores == floor(cores), length(cores) == 1,
                               cores < parallel::detectCores())
                     if ("thin" %in% names(mcmc_setup))
                       stopifnot(mcmc_setup$thin == s("mcmc_setup")$thin)
                     else
                       mcmc_setup$thin = s("mcmc_setup")$thin
                     if ("burnin" %in% names(mcmc_setup) && mcmc_setup$burnin > 0)
                       warning("Burnin is disabled when continuing MCMC procedure.")

                     # Helper parameters
                     betaPos <- 1:length(methods::slot(
                       methods::slot(so("prior"), "beta"), "mu"))
                     deltaPos <- (max(betaPos) + 1):(max(betaPos) +
                                                       length(methods::slot(
                                                         methods::slot(so("prior"),
                                                                       "delta"), "mu")))
                     lambdaPos <- max(deltaPos) + 1
                     chains <- length(s("fit"))
                     lastPoint <- nrow(s("fit")[[1]])
                     lastPoint <- s("fit")[[1]][lastPoint, ]

                     time <- Sys.time()
                     mcmcRun <- list()
                     for (c in 1:chains){
                       if (chains > 1) cat("Starting chain ",c,".\n",sep="")
                       mcmcRun[[c]] <- do.call(cbind,
                                               runBayesPO(
                                                 lastPoint[betaPos], # Starting at last stored point
                                                 lastPoint[deltaPos], # Starting at last stored point
                                                 lastPoint[lambdaPos], # Starting at last stored point
                                                 paste0(so("intensityLink"), "_", # intenisty Link + prior
                                                        methods::slot(
                                                          methods::slot(so("prior"), "beta"), "family")),
                                                 paste0(so("observabilityLink"),"_", # observability Link + prior
                                                        methods::slot(
                                                          methods::slot(so("prior"), "delta"), "family")),
                                                 methods::slot( # lambdaStar prior
                                                   methods::slot(
                                                     so("prior"),"lambdaStar"), "family"),
                                                 retrievePars(methods::slot( # beta prior parameters
                                                   so("prior"),"beta")),
                                                 retrievePars(methods::slot( # delta prior parameters
                                                   so("prior"),"delta")),
                                                 retrievePars( # lambdaStar prior parameters
                                                   methods::slot(
                                                     so("prior"),"lambdaStar")),
                                                 ifelse(is.integer(background), # background class
                                                        "int_mat", "num_mat"),
                                                 background, # background data
                                                 s("area"), # region area
                                                 ifelse(is.integer(so("po")),
                                                        "int_mat", "num_mat"), # po data class
                                                 so("po"), # po data
                                                 backConfig$bIS - 1, # background intensity columns
                                                 backConfig$bOS - 1, # background observability colmns
                                                 so("intensitySelection") - 1, # po intensity columns
                                                 so("observabilitySelection") - 1, # po obserability columns
                                                 0, # MCMC burn-in
                                                 mcmc_setup$thin, # MCMC thin
                                                 mcmc_setup$iter, # MCMC iterations
                                                 cores)
                       )
                       colnames(mcmcRun[[c]]) <- s("parnames")
                       mcmcRun[[c]] <- coda::mcmc(mcmcRun[[c]], thin = mcmc_setup$thin)
                       if (chains > 1) cat("Finished chain ",c,".\n\n",sep="")
                     }
                     if (chains > 1) cat("Total computation time: ", format(unclass(Sys.time()-time),
                                                                            digits = 2), " ",
                                         attr(Sys.time() - time, "units"), ".\n", sep="")

                     return(methods::new("bayesPO_fit",
                                         fit = coda::mcmc.list(
                                           lapply(1:chains, function(i)
                                             coda::mcmc(rbind(s("fit")[[i]],
                                                              mcmcRun[[i]])))
                                         ),
                                         original = s("original"),
                                         backgroundSummary = s("backgroundSummary"),
                                         area = s("area"),
                                         parnames = s("parnames"),
                                         mcmc_setup = list(
                                           burnin = s("mcmc_setup")$burnin,
                                           thin = mcmc_setup$thin,
                                           iter = s("mcmc_setup")$iter + mcmc_setup$iter
                                         )))
                   } )

# Auxiliary function
checkFormatBackground <- function(object,background){
  # Helper function
  s <- function(n) methods::slot(object, n)

  # Intensity
  if (length(s("iSelectedColumns")) > 0){
    backIntensitySelection <- c()
    for (col in s("iSelectedColumns")){
      if (col %in% colnames(background))
        backIntensitySelection <- c(backIntensitySelection, which(col == colnames(background)))
      else stop(paste0("Column ",col," not found in the background covariates"))
    }
  } else
    if (max(s("intensitySelection")) > ncol(background))
      stop(paste0("Requested background column ", max(s(object,"intensitySelection")),
                  " not available in background which only has ",ncol(background)," columns."))
  else backIntensitySelection = s("intensitySelection")

  # Observability
  if (length(s("oSelectedColumns")) > 0){
    backObservabilitySelection <- c()
    for (col in s("oSelectedColumns")){
      if (col %in% colnames(background))
        backObservabilitySelection <- c(backObservabilitySelection, which(col == colnames(background)))
      else stop(paste0("Column ",col," not found in the background covariates"))
    }
  } else
    if (max(s("observabilitySelection")) > ncol(background))
      stop(paste0("Requested background column ", max(s("observabilitySelection")),
                  " not available in background which only has ",ncol(background)," columns."))
  else
    backObservabilitySelection <- s("observabilitySelection")

  return(list(bIS = backIntensitySelection, bOS = backObservabilitySelection))
}

#' @export
bayesPO_model = function(po, intensitySelection,
                         observabilitySelection,
                         intensityLink = "logit", observabilityLink = "logit",
                         initial_values = initial(length(intensitySelection) + 1,
                                                  length(observabilitySelection) + 1,
                                                  nrow(po), random=TRUE),
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
    initial_values = initial(length(intensitySelection) + 1,
                             length(observabilitySelection) + 1,
                             nrow(po), random=TRUE) * initial_values

  return(methods::new("bayesPO_model", po=po, intensityLink = intensityLink,
                      intensitySelection = intensitySelection,
                      observabilityLink = observabilityLink,
                      observabilitySelection = observabilitySelection,
                      init = initial_values, prior = joint_prior,
                      iSelectedColumns = icharSel, oSelectedColumns = ocharSel))
}
