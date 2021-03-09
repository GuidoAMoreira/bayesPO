#### Individual priors ####
#' Retrieve the parameters in an adequate format for the C++ code.
#'
#' @param object The object from which the parameters will be extracted.
methods::setGeneric("retrievePars",function(object){standardGeneric("retrievePars")})

#### Beta-Delta priors ####
#' Generic class for the beta and delta parameters.
#'
#' @slot family The family of distributions of the prior.
#' @export
#' @exportClass BetaDeltaPrior
methods::setClass("BetaDeltaPrior",representation(family="character"))

#' Show method for the BetaDeltaPrior class.
#'
#' @param object the BetaDeltaPrior object.
#' @export
#' @exportMethod show
methods::setMethod("show","BetaDeltaPrior",function(object) cat("Prior for Beta and Delta from the",methods::slot(object,"family"),"family.\n"))

#' Print method for the BetaDeltaPrior class.
#'
#' @param x The BetaDeltaPrior object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print","BetaDeltaPrior",function(x,...) show(x))

#' @method print BetaDeltaPrior
#' @export
print.BetaDeltaPrior <- function(x,...) show(x)

#' Normal prior class for Beta and Delta parameters.
#'
#' @slot mu The mean vector for the prior.
#' @slot Sigma The covariance matrix for the prior.
#' @export
#' @exportClass NormalPrior
methods::setClass("NormalPrior",contains="BetaDeltaPrior",
                  representation = methods::representation(mu="numeric",Sigma="matrix"),
                  validity = function(object){
                    if (!isSymmetric(methods::slot(object,"Sigma"))) stop("Covariance matrix is not symmetric.")
                    if (any(eigen(methods::slot(object,"Sigma"),symmetric=TRUE)$values <= 0)) stop("Covariance matrix is not positive-definite.")
                    if (length(methods::slot(object,"mu")) != nrow(methods::slot(object,"Sigma"))) stop(paste0("Prior mean vector and covariance matrix have incompatible sizes. Mean length is ",length(methods::slot(object,"mu"))," and covariance matrix is ",nrow(methods::slot(object,"Sigma"))," x ",ncol(methods::slot(object,"Sigma")),"."))
                    TRUE
                  })

## Constructor
#' @export
NormalPrior <- function(mu,Sigma){
  mu <- as.numeric(mu)
  new("NormalPrior", mu = mu, Sigma = Sigma, family = "normal")
}

## Methods
#' Names method for the NormalPrior class.
#'
#' @param x The NormalPrior object.
#' @export
#' @exportMethod names
methods::setMethod("names","NormalPrior",function(x) c("mu","Sigma"))

#' The '$' method for the NormalPrior class.
#'
#' @param x The NormalPrior object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$","NormalPrior",function(x,name) methods::slot(x,name))

#' The '$<-' method for the NormalPrior class.
#'
#' @param x The NormalPrior object.
#' @param name The slot to be changed.
#' @param value New value.
#' @export
#' @exportMethod $<-
methods::setMethod("$<-","NormalPrior",function(x,name,value){
  methods::slot(x,name) <- value
  methods::validObject(x)
  x
})

#' Show methodfor the NormalPrior class.
#'
#' @param object The NormalPrior object.
#' @export
#' @exportMethod show
methods::setMethod("show","NormalPrior",function(object){
  cat("Normal prior\n\nMu mean vector:\n")
  print(methods::slot(object,"mu"))
  s = methods::slot(object,"Sigma")
  cat("Sigma covariance matrix:\n")
  if (max(s - diag(diag(s)))) # If Sigma is not a diagonal matrix
    print(methods::slot(object,"Sigma"))
  else if (min(diag(s)) != max(diag(s))) # If Sigma is a diagonal matrix of different values
    cat("diag(c(",paste(diag(s),collapse=","),"))\n",sep="")
  else
    cat(s[1],"* I, where I is an identity matrix.\n")
  invisible(object)
})

#' Print method for the NormalPrior class.
#'
#' @param x The NormalPrior object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print","NormalPrior",function(x,...) methods::show(x))

#' @export
print.NormalPrior <- function(x,...) methods::show(x)

#' Retrieve method for the NormalPrior class.
#'
#' @param object The NormalPrior object.
methods::setMethod("retrievePars","NormalPrior",function(object){
  list(mean = methods::slot(object,"mu"),
       covariance = solve(methods::slot(object,"Sigma")))
})

#### Lambda* priors ####
#' Generic class for the LambdaStar parameters.
#'
#' @slot mu The mean vector for the prior.
#' @export
#' @exportClass LambdaStarPrior
methods::setClass("LambdaStarPrior",methods::representation(family = "character"))

#' Gamma prior class for the LambdaStar parameter.
#'
#' @slot mu The mean vector for the prior.
#' @export
#' @exportClass GammaPrior
methods::setClass("GammaPrior",contains="LambdaStarPrior",
                  representation = methods::representation(shape = "numeric", rate = "numeric"),
                  validity = function(object){
                    for (par in c("shape","rate")){
                      if (length(methods::slot(object,par)) > 1) stop(paste0("Prior parameter ",par," for lambdaStar must have length 1."))
                      if (methods::slot(object,par) <= 0) stop("Prior parameters must be positive")
                    }
                    TRUE
                  })

## Constructor
#' @export
GammaPrior <- function(shape, rate) new("GammaPrior",shape = shape, rate = rate, family = "gamma")

## Methods
#' Names method for the GammaPrior class.
#'
#' @param x The GammaPrior object.
#' @export
#' @exportMethod names
methods::setMethod("names","GammaPrior",function(x) c("shape","rate"))

#' The '$' method for the GammaPrior class.
#'
#' @param x The GammaPrior object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$","GammaPrior",function(x,name) methods::slot(x,name))

#' The '$<-' method for the GammaPrior class.
#'
#' @param x The GammaPrior object.
#' @param name The slot to be changed.
#' @param value New value.
#' @export
#' @exportMethod $<-
methods::setMethod("$<-","GammaPrior",function(x,name,value){
  methods::slot(x,name) <- value
  methods::validObject(x)
  x
})

#' Show method for the GammaPrior class.
#'
#' @param object The GammaPrior object.
#' @export
#' @exportMethod show
methods::setMethod("show","GammaPrior",function(object){
  cat("Gamma prior\n\n")
  vec <- c(methods::slot(object,"shape"),methods::slot(object,"rate"))
  names(vec) <- names(object)
  print(vec)
  invisible(object)
})

#' Print method for the GammaPrior class.
#'
#' @param x The GammaPrior object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print","GammaPrior",function(x,...) methods::show(x))

#' @method print GammaPrior
#' @export
print.GammaPrior <- function(x,...) methods::show(x)

#' Retrieve method for the GammaPrior class.
#'
#' @param object The GammaPrior object.
methods::setMethod("retrievePars","GammaPrior",function(object){
  list(a = methods::slot(object,"shape"),
       b = methods::slot(object,"rate"))
})

#### Joint prior ####
#' Joint prior class for the bayesPO package parameters.
#'
#' @slot beta An object of BetaDeltaPrior S4 class with the appropriate Beta prior.
#' @slot delta An object of BetaDeltaPrior S4 class with the appropriate Delta prior.
#' @slot lambdaStar An object of LambdaStarPrior S4 class with the appropriate LambdaStar prior.
#' @export
#' @exportClass bayesPO_prior
methods::setClass("bayesPO_prior",
                  methods::representation(
                        beta="BetaDeltaPrior",
                        delta="BetaDeltaPrior",
                        lambdaStar="LambdaStarPrior"))

#' Names method for the bayesPO_prior class.
#'
#' @param x The bayesPO_prior object.
#' @export
#' @exportMethod names
methods::setMethod("names","bayesPO_prior",function(x) c("beta","delta","lambdaStar"))

#' The '$' method for the bayesPO_prior class.
#'
#' @param x The bayesPO_prior object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$","bayesPO_prior",function(x,name) methods::slot(x,name))

#' @export
prior <- function(beta,delta,lambdaStar){
  methods::new("bayesPO_prior",beta = beta,delta = delta,lambdaStar = lambdaStar)
}

#' Show method for the bayesPO_prior class.
#'
#' @param object The bayesPO_prior object.
#' @export
#' @exportMethod show
setMethod("show","bayesPO_prior",function(object){
  cat("Joint prior for a bayesPO model. Individual components:\n\nBeta:\n")
      show(methods::slot(object,"beta"))
      cat("\nDelta:\n")
      show(methods::slot(object,"delta"))
      cat("\nlambdaStar:\n")
      print(methods::slot(object,"lambdaStar"))
})

#' Print method for the bayesPO_prior class.
#'
#' @param x The bayesPO_prior object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
setMethod("print","bayesPO_prior",function(x,...) show(x))

#' @method print bayesPO_prior
#' @export
print.bayesPO_prior <- function(x,...) show(x)

#' The '$' method for the bayesPO_prior class.
#'
#' @param x The bayesPO_prior object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
setMethod("$","bayesPO_prior",function(x,name) methods::slot(x,name))

#' The '$<-' method for the bayesPO_prior class.
#'
#' @param x The bayesPO_prior object.
#' @param name The slot to be changed.
#' @param value New value.
#' @export
#' @exportMethod $<-
setMethod("$<-","bayesPO_prior",function(x,name,value){
  methods::slot(x,name) <- value
  methods::validObject(x)
  x
})
