#' @include initial-class.R prior-class.R model-class.R fit-class.R
NULL

#' Fit presence-only data using a Bayesian Poisson Process model
#'
#' The model uses a data augmentation scheme to avoid performing approximations
#' on the likelihood function.
#' @param object Either a \code{bayesPO_model} or \code{bayesPO_fit} object. If
#' a model, then the model is fit according to specifications. If a fit,
#' then the model used to fit the model is recovered and used to continue
#' the MCMC calculations where the previous one left off.
#' @param background A matrix where the rows are the grid cells for the studied
#' region and the columns are the covariates. \code{NA}s must be removed. If
#' the function is being used on a \code{bayesPO_fit} object, the background
#' must be exactly the same as the one used in the original fit.
#' @param mcmc_setup A list containing \code{iter} to inform the model how
#' many iterations are to be run. The list may optionally contain the objects
#' \code{burnin} and \code{thin} to inform these instructions as well.
#' @return An object of class \code{"bayesPO_fit"}.
#' @details The background is kept outside of the
#' @seealso \code{\link{bayesPO_model}} and \code{\link{bayesPO_fit-class}}.
#' @importFrom coda mcmc mcmc.list
#' @export
methods::setGeneric("fit_bayesPO", function(object, background,
                                            mcmc_setup = list(iter = 5000), ...)
  standardGeneric("fit_bayesPO"))

#' @name fit_bayesPO
#' @param area A positive number with the studied region's area.
#' @param cores Currently unused.
#' @importFrom parallel detectCores
#' @export
#' @exportMethod fit_bayesPO
methods::setMethod("fit_bayesPO",
                   methods::signature(object = "bayesPO_model",
                                      background = "matrix"),
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
     if (chains > 1) cat("Starting chain ", c, ".\n",sep="")
     mcmcRun[[c]] <- do.call(cbind,
                             runBayesPO(
                                methods::slot(s("init")[[c]], "beta"),
                                methods::slot(s("init")[[c]], "delta"), # Init delta
                                methods::slot(s("init")[[c]], "lambdaStar"), # Init lambdaStar
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


#' @name fit_bayesPO
#' @param cores Currently unused.
#' @export
#' @exportMethod fit_bayesPO
methods::setMethod("fit_bayesPO", signature(object = "bayesPO_fit",
                                            background = "matrix"),
  function(object, background, mcmc_setup = list(iter = object$mcmc_setup$iter),
           cores = 1){
   # Helper function
   s <- function(n) methods::slot(object, n)
   so <- function(n) methods::slot(s("original"), n)
   backConfig <- checkFormatBackground(s("original"), background)
   cores <- 1

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

#' Build a model to be used in the \code{bayesPO} fitting function
#'
#' Constructor for \code{bayesPO_model-class} objects, built to facilitate
#' the use of the fitting function. The output of this function has the
#' necessary signature for the fit_bayesPO function to start the model fit.
#' @param po A matrix whose rows represent the presence-only data and the
#' columns the covariates observed at each position.
#' @param intensitySelection Either a numeric or character vector and
#' represents the selection of covariates used for the intensity set. If
#' numeric it is the positions of the columns and if character, the names of
#' the columns.
#' @param observabilitySelection Either a numeric or character vector and
#' represents the selection of covariates used for the observability set. If
#' numeric it is the positions of the columns and if character, the names of
#' the columns.
#' @param intensityLink A string to inform what link function the model has
#' with respect to the intensity covariates. Current version accepts 'logit'.
#' @param observabilityLink A string to inform what link function the model has
#' with respect to the observabilitycovariates. Current version accepts 'logit'.
#' @param initial_values Either a single integer or list containing
#' \code{bayesPO_initial-class} objects. The length of the list will inform the
#' model how many independent chains will be run. If an integer, that many
#' initial values will be randomly generated.
#' @param prior A \code{bayesPO_prior} object.
#' @seealso \code{\link{initial}}, \code{\link{prior}} and
#' \code{\link{fit_bayesPO}}.
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

  stopifnot(is.matrix(po),
            is.numeric(intensitySelection) || is.character(intensitySelection),
            is.character(intensityLink), length(intensityLink) == 1,
            is.character(observabilityLink), length(observabilityLink) == 1,
            (is.list(initial_values) && all(do.call(c, lapply(initial_values, is, "bayesPO_initial"))))
            || (is.numeric(initial_values) && length(initial_values) == 1))

  if (is.numeric(intensitySelection)) icharSel = character() else
    if (is.character(intensitySelection)){ # Create columns positions for intensity
      icharSel = intensitySelection
      intensitySelection = c()
      for (selected in icharSel){
        if (!(selected %in% colnames(po))){
          message(paste0("Column ",selected," not found in the observed points covariates. Ignoring it."))
          next
        }
        intensitySelection = c(intensitySelection,which(selected == colnames(po)))
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

  return(methods::new("bayesPO_model", po = po, intensityLink = intensityLink,
                      intensitySelection = intensitySelection,
                      observabilityLink = observabilityLink,
                      observabilitySelection = observabilitySelection,
                      init = initial_values, prior = joint_prior,
                      iSelectedColumns = icharSel, oSelectedColumns = ocharSel))
}
