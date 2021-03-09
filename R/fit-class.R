methods::setOldClass("mcmc.list")
methods::setOldClass("table")
methods::setOldClass("list")

#' Class for the result of the MCMC procedure.
#'
#' @slot fit The actual fit from the model. It is of class \code{\link[coda]{mcmc.list}}, as generated from the \code{coda} package.
#' @slot original The model used to generate the chains, an object with class \code{bayesPO_model}.
#' @slot backgroundSummary A small summary of the original background covariates. This is to ensure that continuing the chains will use the identical background matrix.
#' @slot area A positive number indicating the area measure of the region being studied.
#' @slot parnames The names of the parameters. If the model used selects the covariates with column names, they are replicated here. If they are the column indexes, names are generated for identification.
#' @slot mcmc_setup The original mcmc setup used.
#' @include model-class.R
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
#' Show method for the bayesPO_fit class.
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
      setup$burnin," warmup ",ifelse(setup$burnin>1,"iterations","iteration"),
      " and ",setup$n_iter," valid ",
      ifelse(setup$n_iter>1,"iterations","iteration"),", storing one in every ",
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

#' Print method
#'
#' @param x A bayesPO_fit object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print","bayesPO_fit",function(x,...) methods::show(x))

#' @method print bayesPO_fit
#' @export
print.bayesPO_fit <- function(x,...) methods::show(x)

#' Summary method for the bayesPO_fit class.
#'
#' @param object A bayesPO_fit object.
#' @param ... Ignored.
#' @return A matrix with the summary.
#' @export
#' @exportMethod summary
methods::setMethod("summary","bayesPO_fit",function(object,...) summary.bayesPO_fit(object,...))

#' @export
summary.bayesPO_fit <- function(object,...){
  chains <- length(methods::slot(methods::slot(object,"original"),"init"))
  nb <- length(methods::slot(methods::slot(object,"original"),"intensitySelection"))+1
  nd <- length(methods::slot(methods::slot(object,"original"),"observabilitySelection"))+1
  npar <- nb+nd+1 # +1 = lambdaStar
  result <- matrix(0,npar,6) # Mean, median, sd, lower CI bound, upper CI bound, effective sample size
  colnames(result) <- c("mean","median","sd","2.5%","97.5%","eff. sample size")
  rownames(result) <- methods::slot(object,"parnames")[1:npar]
  fitToMatrix <- as.matrix(methods::slot(object,"fit"))[,1:npar]
  result[,1] <- colMeans(fitToMatrix)
  result[,2] <- apply(fitToMatrix,2,stats::median)
  result[,3] <- apply(fitToMatrix,2,stats::sd)
  result[,4] <- apply(fitToMatrix,2,stats::quantile,0.025)
  result[,5] <- apply(fitToMatrix,2,stats::quantile,0.975)
  result[,6] <- coda::effectiveSize(methods::slot(object,"fit"))[1:npar]
  if (chains > 1){
    cols <- colnames(result)
    result <- cbind(result,coda::gelman.diag(methods::slot(object,"fit"))$psrf[1:npar,])
    colnames(result) <- c(cols,"Estimated Rhat","Upper CI Rhat")
  }
  result
}

#' Names method for the bayesPO_fit class.
#'
#' @param x A bayesPO_fit object.
#' @export
#' @exportMethod names
methods::setMethod("names","bayesPO_fit",function(x){
  nn <- c("parameters","mcmc_chains","model","log_posterior","eff. sample size","area","initial values","mcmc setup")
  if (length(methods::slot(methods::slot(x,"original"),"init")) > 1)
    nn <- c(nn,"Rhat","Rhat_upper_CI")
  nn
})

#' The '[[' method for the bayesPO_fit class.
#'
#' @param x A bayesPO_fit object.
#' @param i The requested slot.
#' @export
#' @exportMethod [[
methods::setMethod("[[","bayesPO_fit",function(x,i){
  nb <- length(methods::slot(methods::slot(x,"original"),"intensitySelection"))+1
  nd <- length(methods::slot(methods::slot(x,"original"),"observabilitySelection"))+1
  npar <- nb+nd+1 # +1 from lambdaStar
  summ <- summary(x)

  if (i == "parameters"){
    output <- summ[,1]
    names(output) <- rownames(summ)
  } else
  if (i == "eff. sample size"){
    output <- summ[,6]
    names(output) <- rownames(summ)
  } else
  if (i == "Rhat"){
    output <- summ[,7]
    names(output) <- rownames(summ)
  } else
  if (i == "Rhat_upper_CI"){
    output <- summ[,8]
    names(output) <- rownames(summ)
  } else
  if (i == "mcmc_chains") output <- methods::slot(x,"fit") else
  if (i == "model") output <- methods::slot(x,"original") else
  if (i == "initial values") output <- methods::slot(methods::slot(x,"original"),"init") else
  if (i == "mcmc setup") output <- methods::slot(x,"mcmc_setup") else
  if (i == "log_posterior") output <- as.data.frame(x)$log_Posterior else
  if (i == "area") output <- methods::slot(x,"area")

  output
})

#' The '$' method for the bayesPO_fit class.
#'
#' @param x A bayesPO_fit object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$","bayesPO_fit",function(x,name) x[[name]])

## as.array
#' as.array method for the bayesPO_fit class.
#'
#' @param x A bayesPO_fit object.
#' @param ... Ignored.
#' @return The MCMC chains organized in a way ready for the \code{bayesplot} package.
#' @export
#' @exportMethod as.array
methods::setMethod("as.array","bayesPO_fit",function(x,...) as.array.bayesPO_fit(x,...))

#' Put MCMC output in an \code{array}.
#'
#' Prepares the MCMC results to be used by plotting functions of the
#' bayesplot package.
#' @param x Output of the \code{\link{fit_bayesPO}} function.
#' @param ... Ignored in this version.
#' @return An \code{array} with dimensions I x C x P, where I stands for number of
#' iterations, C for number of chains and P for total number of parameters.
#' P is actually larger than the number of parameters in the model, as the
#' the generated sizes of the latent processes and the log-posterior are also
#' included.
#' @seealso \code{\link{fit_bayesPO}}
#' @method as.array bayesPO_fit
#' @export
as.array.bayesPO_fit <- function(x,...){
  nchains <- length(methods::slot(x,"fit"))
  chains <- do.call(rbind,methods::slot(x,"fit"))
  iterations <- nrow(methods::slot(x,"fit")[[1]])
  npar <- ncol(chains)

  ## Format to be used with bayesplot:: functions
  return(
    array(chains,
          dim = c(iterations,nchains,npar),
          dimnames = list(iterations=NULL,
                          chains=paste0("chain:",1:nchains),
                          parameters=methods::slot(x,"parnames")))
  )
}

## as.matrix
#' as.matrix method for the bayesPO_fit class.
#'
#' @param x A bayesPO_fit object.
#' @param ... Ignored.
#' @return The MCMC chains reformatted as a matrix.
#' @export
#' @exportMethod as.matrix
methods::setMethod("as.matrix","bayesPO_fit",function(x,...) as.matrix.bayesPO_fit(x,...))

#' Put MCMC output in a \code{matrix}.
#'
#' Create a matrix with all the MCMC results.
#' @param x Output of the \code{\link{fit_bayesPO}} function.
#' @param ... Ignored in this version.
#' @return A matrix where all the MCMC results are included.
#' @details The dimension of the output is I*C x P+2, where I stands for
#' number of iterations, C for number of chains and P for total number of
#' parameters. P is actually larger than the number of parameters in the model,
#' as the generated sizes of the latent processes and the log-posterior are
#' also included.
#'
#' Two extra columns are included to indicate to which chain and to which
#' iteration that draw belongs.
#' @seealso \code{\link{fit_bayesPO}}
#' @method as.matrix bayesPO_fit
#' @export
as.matrix.bayesPO_fit <- function(x,...){
  nchains <- length(methods::slot(x,"fit"))
  chains <- do.call(rbind,methods::slot(x,"fit"))
  iterations <- nrow(methods::slot(x,"fit")[[1]])
  parnames <- colnames(chains)
  chains <- cbind(chains,rep(factor(1:nchains),each=iterations))
  chains <- cbind(chains,rep(1:iterations,nchains))
  colnames(chains) <- c(parnames,"chain","iteration")

  return(chains)
}

## as.data.frame
#' as.data.frame method for the bayesPO_fit class.
#'
#' @param x A bayesPO_fit object.
#' @param row.names NULL or a character vector giving the row names for the
#' data frame. Missing values are not allowed.
#' @param optional logical. If TRUE, setting row names and converting column
#' names to syntactic names is optional. See help('as.data.frame') for more.
#' Leaving as \code{FALSE} is recommended.
#' @param ... Ignored.
#' @return The MCMC chains reformatted as a data.frame
#' @export
#' @exportMethod as.data.frame
methods::setMethod("as.data.frame","bayesPO_fit",function(x, row.names = NULL, optional = FALSE, ...) as.data.frame.bayesPO_fit(x, row.names = NULL, optional = FALSE, ...))

#' Put MCMC output in a \code{data.frame}.
#'
#' Create a \code{data.frame} with all the MCMC results.
#' @param x Output of the \code{\link{fit_bayesPO}} function.
#' @param row.names NULL or a character vector giving the row names for the
#' data frame. Missing values are not allowed.
#' @param optional logical. If TRUE, setting row names and converting column
#' names to syntactic names is optional. See help('as.data.frame') for more.
#' Leaving as \code{FALSE} is recommended.
#' @param ... Ignored in this version.
#' @return A data.frame where all the MCMC results are included.
#' @details The dimension of the output is I*C x P+2, where I stands for
#' number of iterations, C for number of chains and P for total number of
#' parameters. P is actually larger than the number of parameters in the model,
#' as the generated sizes of the latent processes and the log-posterior are
#' also included.
#'
#' Two extra columns are included to indicate to which chain and to which
#' iteration that draw belongs. This is to facilitate the use of plotting
#' results via the \code{ggplot2} package if desired.
#'
#' If \code{row.names} is left at \code{NULL} then row names are created as
#' CcIi where c is the chain and i is the iteration of that row.
#' @seealso \code{\link{fit_bayesPO}}
#' @method as.data.frame bayesPO_fit
#' @export
as.data.frame.bayesPO_fit = function(x, row.names = NULL, optional = FALSE, ...){
  nchains <- length(methods::slot(x,"fit"))
  chains <- do.call(rbind,methods::slot(x,"fit"))
  parnames <- colnames(chains)
  iterations <- nrow(methods::slot(x,"fit")[[1]])

  colsList <- list()
  for (pp in 1:length(parnames)) colsList[[parnames[pp]]] <- chains[,pp]

  if (!is.null(NULL)) row_names <- row.names else row_names <- paste0(rep(paste0("C",1:nchains),each=iterations),rep(paste0("I",1:iterations),nchains))
  output <- do.call(data.frame,c(colsList,list(check.names = TRUE, fix.empty.names = TRUE),list(row.names = row_names)))
  output$chain <- factor(rep(1:nchains,each=iterations))
  output$iteration <- rep(1:iterations,nchains)

  return(output)
}

