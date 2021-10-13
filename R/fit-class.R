methods::setOldClass("mcmc.list")
methods::setOldClass("table")
methods::setOldClass("list")

#' @include model-class.R
NULL

#' Class for the result of the MCMC procedure.
#'
#' Objects of this class are the main objects of this package. They contain
#' much information about the fitted model.
#' @slot fit The actual fit from the model. It is an object of class
#' \code{\link[coda]{mcmc.list}}, as generated from the \code{coda} package.
#' @slot original The model used to generate the chains, an object with class
#' \code{bayesPO_model}.
#' @slot backgroundSummary A small summary of the original background
#' covariates. This is to ensure that continuing the chains will use the
#' identical background matrix. Only the summary is kept for storage efficiency.
#' @slot area A positive number indicating the area measure of the region being
#' studied.
#' @slot parnames The names of the parameters. If the model used selects the
#' covariates with column names, they are replicated here. If they are the
#' column indexes, names are generated for identification.
#' @slot mcmc_setup The original mcmc setup used.
#' @seealso \code{\link{fit_bayesPO}}
#' @export
#' @exportClass bayesPO_fit
methods::setClass("bayesPO_fit",
                  representation(fit = "mcmc.list",
                                 original = "bayesPO_model",
                                 backgroundSummary = "table",
                                 area = "numeric",
                                 parnames = "character",
                                 mcmc_setup = "list"))

#### Basic methods ####
#' @name bayesPO_fit-class
#'
#' @param object A bayesPO_fit object.
#' @export
#' @exportMethod show
methods::setMethod("show","bayesPO_fit",function(object){
  cat("Fit of a bayesPO model.\n")

  ## data
  cat("Presence-only dataset size:",nrow(methods::slot(methods::slot(object,"original"),"po")),"\n\n")
  sc <- methods::slot(methods::slot(object,"original"),"iSelectedColumns")
  if (length(sc))
    cat(length(sc)," intensity covariates selected:\n",sc,"\n",sep = "")
  else{
    sc <- methods::slot(methods::slot(object,"original"),"intensitySelection")
    cat(length(sc)," intensity covariates selected. Columns: ",sc,"\n",sep = "")
  }
  sc <- methods::slot(methods::slot(object,"original"),"oSelectedColumns")
  if (length(sc))
    cat(length(sc)," observability covariates selected:\n",sc,"\n",sep="")
  else{
    sc <- methods::slot(methods::slot(object,"original"),"observabilitySelection")
    cat(length(sc)," observability covariates selected. Columns: ",sc,"\n",sep = "")
  }
  cat("\n")

  ## Link function
  links <- c(methods::slot(methods::slot(object,"original"),"intensityLink"),
             methods::slot(methods::slot(object,"original"),"observabilityLink"))
  names(links) <- c("intensity","observability")
  cat("Link functions chosen:\n")
  print(links)
  cat("\n")

  ## Prior
  cat("Prior selection:\n")
  methods::show(methods::slot(methods::slot(object,"original"),"prior"))
  cat("\n")

  ## MCMC configuration
  chains = length(methods::slot(methods::slot(object,"original"),"init"))
  setup = methods::slot(object,"mcmc_setup")
  cat(chains,ifelse(chains>1," chains"," chain")," of MCMC ",
      ifelse(chains>1,"were","was")," configured with ",
      format(setup$burnin, scientific = FALSE)," warmup ",
      ifelse(setup$burnin>1,"iterations","iteration"),
      " and ", format(setup$iter, scientific = FALSE)," valid ",
      ifelse(setup$iter>1,"iterations","iteration"),", storing one in every ",
      ifelse(setup$thin>1,paste(setup$thin,"steps"),"step"),".\n\n",sep="")

  ## Results
  print(round(summary(object),digits = 3))
  cat("\n")

  ## Comments
  cat("The effective sample size represents the sample size of an independent",
      "sample which would yield equivalent estimates as the autocorrelated",
      "MCMC result.\n")
  if (chains > 1)
    cat("Rhat has been calculated from multiple chains. Lower, closer to 1",
        "values indicate better convergence of the chain. For posterior",
        "estimates to be trusted, the upper CI limit should be below 1.1. If",
        "they are not, run more iterations. See help('fit_bayesPO') to see",
        "how to do utilize the iterations already run.\n\n")
  else
    cat("Rhat cannot be estimated with only 1 chain. Run more chains for this",
        "statistic to be displayed.")

  invisible(object)
})

#' @name bayesPO_fit-class
#'
#' @param x A bayesPO_fit object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print", "bayesPO_fit", function(x, ...) methods::show(x))

#' @method print bayesPO_fit
#' @export
print.bayesPO_fit <- function(x,...) methods::show(x)

#' @name bayesPO_fit-class
#'
#' @param object A bayesPO_fit object.
#' @param ... Ignored.
#' @return \strong{\code{summary}}: A matrix with the summary.
#' @export
#' @exportMethod summary
methods::setMethod("summary", "bayesPO_fit", function(object,...) summary.bayesPO_fit(object, ...))

#' @method summary bayesPO_fit
#' @export
summary.bayesPO_fit <- function(object, ...){
  chains <- length(methods::slot(methods::slot(object, "original"), "init"))
  nb <- length(methods::slot(methods::slot(object, "original"), "intensitySelection")) + 1
  nd <- length(methods::slot(methods::slot(object, "original"), "observabilitySelection")) + 1
  npar <- nb + nd + 1 # +1 = lambdaStar
  result <- matrix(0, npar, 6) # Mean, median, sd, lower CI bound, upper CI bound, effective sample size
  colnames(result) <- c("mean", "median", "sd", "2.5%", "97.5%", "eff. sample size")
  rownames(result) <- methods::slot(object, "parnames")[1:npar]
  fitToMatrix <- as.matrix(methods::slot(object, "fit"))[, 1:npar]
  result[,1] <- colMeans(fitToMatrix)
  result[,2] <- apply(fitToMatrix, 2, stats::median)
  result[,3] <- apply(fitToMatrix, 2, stats::sd)
  result[,4] <- apply(fitToMatrix, 2, stats::quantile, 0.025)
  result[,5] <- apply(fitToMatrix, 2, stats::quantile, 0.975)
  result[,6] <- coda::effectiveSize(methods::slot(object,"fit"))[1:npar]
  if (chains > 1){
    cols <- colnames(result)
    result <- cbind(result, coda::gelman.diag(methods::slot(object, "fit"))$psrf[1:npar, ])
    colnames(result) <- c(cols, "Estimated Rhat","Upper CI Rhat")
  }
  result
}

#' @name bayesPO_fit-class
#'
#' @param x A bayesPO_fit object.
#' @export
#' @exportMethod names
methods::setMethod("names", "bayesPO_fit", function(x){
  nn <- c("parameters", "covariates_importance", "mcmc_chains", "model",
          "log_posterior", "eff_sample_size", "area", "initial_values", "mcmc_setup")
  if (length(methods::slot(methods::slot(x, "original"), "init")) > 1)
    nn <- c(nn, "Rhat", "Rhat_upper_CI")
  nn
})

#' @method names bayesPO_fit
#' @export
names.bayesPO_fit <- function(x) names(x)

#' @name bayesPO_fit-class
#'
#' @param x A bayesPO_fit object.
#' @param i The requested slot.
#' @export
#' @exportMethod [[
methods::setMethod("[[", "bayesPO_fit", function(x, i){
  # Helper function
  s <- function(n) methods::slot(x, n)

  nb <- length(methods::slot(s("original"), "intensitySelection")) + 1
  nd <- length(methods::slot(s("original"), "observabilitySelection")) + 1
  npar <- nb + nd + 1 # +1 from lambdaStar

  if (i == "parameters"){
    summ <- summary(x)
    output <- summ[, 1]
    names(output) <- rownames(summ)
  } else
  if (i == "covariates_importance"){
    data <- as.data.frame(x)

    intensity <- t(apply(
      data[2:(which(names(data) == "Observability_Intercept") - 1)], 1,
      function(chain) {c2 <- chain * chain; c2 / sum(c2)}
    ))
    observability <- t(apply(
      data[(which(names(data) == "Observability_Intercept") + 1):(which(names(data) == "lambdaStar") - 1)], 1,
      function(chain) {c2 <- chain * chain; c2 / sum(c2)}
    ))
    colnames(intensity) <- names(data)[2:(which(names(data) == "Observability_Intercept") - 1)]
    colnames(observability) <- names(data)[(which(names(data) == "Observability_Intercept") + 1):(which(names(data) == "lambdaStar") - 1)]
    output <- list(intensity = intensity, observability = observability)
    class(output) <- "covariates_importance"
  }
  if (i == "eff_sample_size"){
    summ <- summary(x)
    output <- summ[, 6]
    names(output) <- rownames(summ)
  } else
  if (i == "Rhat"){
    summ <- summary(x)
    output <- summ[, 7]
    names(output) <- rownames(summ)
  } else
  if (i == "Rhat_upper_CI"){
    summ <- summary(x)
    output <- summ[, 8]
    names(output) <- rownames(summ)
  } else
  if (i == "mcmc_chains") output <- s("fit") else
  if (i == "model") output <- s("original") else
  if (i == "initial_values") output <- methods::slot(s("original"), "init") else
  if (i == "mcmc_setup") output <- s("mcmc_setup") else
  if (i == "log_posterior") output <- as.data.frame(x)$log_Posterior else
  if (i == "area") output <- s("area")

  output
})

#' @name bayesPO_fit-class
#'
#' @param x A bayesPO_fit object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$", "bayesPO_fit", function(x, name) x[[name]])

#' @name bayesPO_fit-class
#' @param x A bayesPO_fit object.
#' @param ... Ignored.
#' @return \strong{\code{as.array}}: The MCMC chains organized in a way ready for the
#' \code{bayesplot} package.
#' @export
#' @exportMethod as.array
methods::setMethod("as.array", "bayesPO_fit", function(x, ...) as.array.bayesPO_fit(x, ...))

namesAid <- function(string){
  new_string <- string

  intInt <- "(Intensity intercept)"
  obsInt <- "(Observability intercept)"
  obsStart <- max(which(string == obsInt), which(string == "delta_0"))
  obsEnd <- which(string == "lambdaStar") - 1

  # Find same covariates
  for (i in 1:(obsStart - 1))
    searching <- grepl(string[i], string[obsStart:obsEnd])
    if (any(searching)){
      new_string[i] <- paste0(string[i], ".int")
      new_string[obsStart - 1 + which(searching)] <- paste0(string[i], ".obs")
    }
  new_string[1] <- ifelse(string[1] == obsInt, "Intensity_Intercept", "beta_0")
  new_string[obsStart] <- ifelse(string[obsStart] == obsInt, "Observability_Intercept", "delta_0")

  new_string
}

#' @name bayesPO_fit-class
#' @param x A bayesPO_fit object.
#' @param ... Ignored in this version.
#' @return \strong{\code{as.array}}: An \code{array} with dimensions I x C x P, where
#' I stands for number of iterations, C for number of chains and P for total
#' number of parameters. P is actually larger than the number of parameters in
#' the model, as the the generated sizes of the latent processes and the
#' log-posterior are also included.
#' @method as.array bayesPO_fit
#' @export
as.array.bayesPO_fit <- function(x, ...){
  nchains <- length(methods::slot(x, "fit"))
  chains <- do.call(rbind, methods::slot(x, "fit"))
  iterations <- nrow(methods::slot(x, "fit")[[1]])
  npar <- ncol(chains)

  ## Format to be used with bayesplot:: functions
  return(
    array(chains,
          dim = c(iterations, nchains, npar),
          dimnames = list(iterations = NULL,
                          chains = paste0("chain:", 1:nchains),
                          parameters = namesAid(methods::slot(x, "parnames"))))
  )
}

#' @export
#' @exportMethod as.matrix
methods::setMethod("as.matrix", "bayesPO_fit", function(x, ...) as.matrix.bayesPO_fit(x, ...))

#' @name bayesPO_fit-class
#' @param x A bayesPO_fit object.
#' @param ... Ignored in this version.
#' @return \strong{\code{as.matrix}}: The dimension of the output is I * C x (P + 2),
#' where I stands for number of iterations, C for number of chains and P for
#' total number of parameters. P is actually larger than the number of
#' parameters in the model, as the generated sizes of the latent processes and
#' the log-posterior are also included.
#'
#' Two extra columns are included to indicate to which chain and to which
#' iteration that draw belongs.
#' @method as.matrix bayesPO_fit
#' @export
as.matrix.bayesPO_fit <- function(x, ...){
  nchains <- length(methods::slot(x, "fit"))
  chains <- do.call(rbind, methods::slot(x, "fit"))
  iterations <- nrow(methods::slot(x, "fit")[[1]])
  parnames <- namesAid(colnames(chains))
  chains <- cbind(chains, rep(factor(1:nchains), each = iterations))
  chains <- cbind(chains, rep(1:iterations, nchains))
  colnames(chains) <- c(parnames, "chain", "iteration")

  return(chains)
}

#' @export
#' @exportMethod as.data.frame
methods::setMethod("as.data.frame","bayesPO_fit",function(x, row.names = NULL, optional = FALSE, ...) as.data.frame.bayesPO_fit(x, row.names = NULL, optional = FALSE, ...))

#' @name bayesPO_fit-class
#' @param x A bayesPO_fit object.
#' @param row.names NULL or a character vector giving the row names for the
#' data frame. Missing values are not allowed.
#' @param optional logical. If TRUE, setting row names and converting column
#' names to syntactic names is optional. See help('as.data.frame') for more.
#' Leaving as \code{FALSE} is recommended.
#' @param ... Ignored in this version.
#' @return \strong{\code{as.data.frame}}: The dimension of the output is I*C x P+2,
#' where I stands for number of iterations, C for number of chains and P for
#' total number of parameters. P is actually larger than the number of
#' parameters in the model, as the generated sizes of the latent processes and
#' the log-posterior are also included.
#'
#' Two extra columns are included to indicate to which chain and to which
#' iteration that draw belongs. This is to facilitate the use of plotting
#' results via the \code{ggplot2} package if desired.
#'
#' If \code{row.names} is left at \code{NULL} then row names are created as
#' CcIi where c is the chain and i is the iteration of that row.
#' @method as.data.frame bayesPO_fit
#' @export
as.data.frame.bayesPO_fit = function(x, row.names = NULL, optional = FALSE, ...){
  nchains <- length(methods::slot(x, "fit"))
  chains <- do.call(rbind, methods::slot(x, "fit"))
  parnames <- namesAid(colnames(chains))
  iterations <- nrow(methods::slot(x, "fit")[[1]])

  colsList <- list()
  for (pp in 1:length(parnames)) colsList[[parnames[pp]]] <- chains[, pp]

  if (!is.null(NULL))
    row_names <- row.names else
      row_names <- paste0(rep(paste0("C", 1:nchains),
                              each=iterations),
                          rep(paste0("I", 1:iterations), nchains))
  output <- do.call(data.frame, c(colsList, list(check.names = TRUE,
                                                 fix.empty.names = TRUE),
                                  list(row.names = row_names)))
  output$chain <- factor(rep(1:nchains, each=iterations))
  output$iteration <- rep(1:iterations, nchains)

  return(output)
}

#### Interaction methods ####
# Adding chains into a single object
#' @name bayesPO_fit-class
#' @param e1 A bayesPO_fit object.
#' @param e2 A bayesPO_fit object with the same background, model (except for
#' initial values), area, parnames and mcmc_setup as \code{e1}.
#' @return \strong{\code{+}}: A new \code{bayesPO_fit} object where the chains
#' are combined into a new multi-chain object.
#' @importFrom methods new
methods::setMethod("+", methods::signature(e1 = "bayesPO_fit", e2 = "bayesPO_fit"),
                   function(e1, e2){
  s1 <- function(n) methods::slot(e1, n)
  s2 <- function(n) methods::slot(e2, n)
  stopifnot(#all.equal(s1("original"), s2("original")),
            all.equal(s1("backgroundSummary"), s2("backgroundSummary")),
            all.equal(s1("area"), s2("area")),
            all.equal(s1("parnames"), s2("parnames")),
            all.equal(s1("mcmc_setup"), s2("mcmc_setup")))

  fff <- coda::mcmc.list(c(s1("fit"), s2("fit")))

  or <- s1("original")
  methods::slot(or, "init") <- c(methods::slot(s1("original"), "init"),
                                 methods::slot(s2("original"), "init"))

  return(methods::new("bayesPO_fit",
                      fit = fff,
                      original = or,
                      backgroundSummary = s1("backgroundSummary"),
                      area = s1("area"),
                      parnames = s1("parnames"),
                      mcmc_setup = s1("mcmc_setup")))
})

# Combining multiple chains
#' @name bayesPO_fit-class
#' @return \strong{\code{c}}: A new \code{bayesPO_fit} object where the chains
#' are combined into a new multi-chain object.
#' @export
#' @exportMethod c
methods::setMethod("c", "bayesPO_fit", function(x, ...) {
  ll <- list(...)
  res <- x
  for (i in 1:length(ll))
    res <- res + ll[[i]]

  res
})

