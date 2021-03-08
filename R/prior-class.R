#' @exportClass bayesPO_prior
methods::setClass("bayesPO_prior",
                  methods::representation(
                        beta="BetaDeltaPrior",
                        delta="BetaDeltaPrior",
                        lambdaStar="LambdaStarPrior"))

#' @export
methods::setMethod("names","bayesPO_prior",function(x) c("beta","delta","lambdaStar"))

#' @export
methods::setMethod("$","bayesPO_prior",function(x,name) methods::slot(x,name))

#' @export
prior <- function(beta,delta,lambdaStar){
  methods::new("bayesPO_prior",beta = beta,delta = delta,lambdaStar = lambdaStar)
}

#' @export
setMethod("show","bayesPO_prior",function(object){
  cat("Joint prior for a bayesPO model. Individual components:\n\nBeta:\n")
      show(methods::slot(object,"beta"))
      cat("\nDelta:\n")
      show(methods::slot(object,"delta"))
      cat("\nlambdaStar:\n")
      print(methods::slot(object,"lambdaStar"))
})
#' @export
setMethod("print","bayesPO_prior",function(x,...) show(x))

#' @export
print.bayesPO_prior <- function(x,...) show(x)

#' @export
setMethod("$","bayesPO_prior",function(x,name) methods::slot(x,name))

#' @export
setMethod("$<-","bayesPO_prior",function(x,name,value){
  methods::slot(x,name) <- value
  methods::validObject(x)
  x
})

#### Individual priors ####
# Get parameters in a list to make it easier to build C++ code
methods::setGeneric("retrievePars",function(object){standardGeneric("retrievePars")})

#### Beta-Delta priors ####
#' @exportClass BetaDeltaPrior
methods::setClass("BetaDeltaPrior",representation(family="character"))

#' @export
methods::setMethod("show","BetaDeltaPrior",function(object) cat("Prior for Beta and Delta from the",methods::slot(object,"family"),"family.\n"))

#' @export
methods::setMethod("print","BetaDeltaPrior",function(x,...) show(x))

#' @export
print.BetaDeltaPrior <- function(x,...) show(x)

#' @export
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
#' @export
methods::setMethod("names","NormalPrior",function(x) c("mu","Sigma"))

#' @export
methods::setMethod("$","NormalPrior",function(x,name) methods::slot(x,name))

#' @export
methods::setMethod("$<-","NormalPrior",function(x,name,value){
  methods::slot(x,name) <- value
  methods::validObject(x)
  x
})

#' @export
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

#' @export
methods::setMethod("print","NormalPrior",function(x,...) methods::show(x))

#' @export
print.NormalPrior <- function(x,...) methods::show(x)

methods::setMethod("retrievePars","NormalPrior",function(object){
  list(mean = methods::slot(object,"mu"),
       covariance = solve(methods::slot(object,"Sigma")))
})

#### Lambda* priors ####
#' @exportClass LambdaStarPrior
methods::setClass("LambdaStarPrior",methods::representation(family = "character"))

#' @export
methods::setClass("GammaPrior",contains="LambdaStarPrior",
                  representation = methods::representation(shape = "numeric", rate = "numeric"),
         validity = function(object){
           for (par in c("shape","rate")){
             if (length(methods::slot(object,par)) > 1) stop(paste0("Prior parameter ",par," for lambdaStar must have length 1."))
             if (methods::slot(object,par) <= 0) stop("Prior parameters must be positive")
           TRUE
         })

## Constructor
#' @export
GammaPrior <- function(shape, rate) new("GammaPrior",shape = shape, rate = rate, family = "gamma")

## Methods
#' @export
methods::setMethod("names","GammaPrior",function(x) c("shape","rate"))

#' @export
methods::setMethod("$","GammaPrior",function(x,name) methods::slot(x,name))

#' @export
methods::setMethod("$<-","GammaPrior",function(x,name,value){
  methods::slot(x,name) <- value
  methods::validObject(x)
  x
  })

#' @export
methods::setMethod("show","GammaPrior",function(object){
  cat("Gamma prior\n\n")
  vec <- c(methods::slot(object,"shape"),methods::slot(object,"rate"))
  names(vec) <- names(object)
  print(vec)
  invisible(object)
})

#' @export
methods::setMethod("print","GammaPrior",function(x,...) methods::show(x))

#' @export
print.GammaPrior <- function(x,...) methods::show(x)

methods::setMethod("retrievePars","GammaPrior",function(object){
  list(a = methods::slot(object,"shape"),
       b = methods::slot(object,"rate"))
})
