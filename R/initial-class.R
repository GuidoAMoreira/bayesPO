#' Class for the initial values for the MCMC for the bayesPO package
#'
#' @slot beta Initial values for beta.
#' @slot delta Initial values for delta.
#' @slot lambdaStar Initial values for lambdaStar.
#' @slot tag Indicates the source of the initial values.
methods::setClass("bayesPO_initial",
                  methods::representation(beta = "numeric",
                        delta = "numeric",
                        lambdaStar = "numeric",
                        tag = "character"),
         validity = function(object){
           if (length(methods::slot(object,"beta")) < 1) stop("Beta must have at least 1 initial value")
           if (length(methods::slot(object,"delta")) < 1) stop("Delta must have at least 1 initial value")
           if (length(methods::slot(object,"lambdaStar")) != 1) stop("lambdaStar must have exactly 1 initial value")
           TRUE
         })

#' @export
methods::setMethod("initialize","bayesPO_initial",function(.Object,beta,delta,lambdaStar,tag){
  methods::slot(.Object,"beta") <- as.numeric(beta)
  methods::slot(.Object,"delta") <- as.numeric(delta)
  methods::slot(.Object,"lambdaStar") <- as.numeric(lambdaStar)
  methods::slot(.Object,"tag") <- tag
  .Object
})

## Methods
#' @export
methods::setMethod("names","bayesPO_initial",function(x) c("beta","delta","lambdaStar"))

#' @export
methods::setMethod("$","bayesPO_initial",function(x,name) methods::slot(x,name))

#' @export
methods::setMethod("+","bayesPO_initial",function(e1,e2) list(e1,e2))

#' @export
methods::setMethod("+",methods::signature(e1 = "bayesPO_initial", e2 = "list"),
                   function(e1,e2){
                     for (i in 1:length(e2))
                      if (!is(e2[[i]],"bayesPO_initial")) stop("Initial values can only be added to a list of other initial values.")
                     e2[[i+1]] = e1
                     e2
                   })
#' @export
methods::setMethod("+",methods::signature(e1 = "list", e2 = "bayesPO_initial"),
                   function(e1,e2){
                     n = length(e1)
                     l = list(e2)
                     for (i in 1:n){
                       if (!is(e1[[i]],"bayesPO_initial")) stop("Initial values can only be added to a list of other initial values.")
                       l[[i+1]] = e1[[i]]
                     }
                     l
                   })
#' @export
methods::setMethod("*",methods::signature(e1 = "bayesPO_initial",e2 = "numeric"),
                   function(e1,e2) {
                     if (e2 <= 0) stop("Can only multiply by a positive number.")
                     if (as.integer(e2) != e2) stop("Con only multiply by an integer.")
                     if (methods::slot(e1,"tag") != "random") {
                       message("Identical initial values are not recommended for independent Markov Chains.")
                       l = list()
                       for (i in 1:e2) l[[i]] = e1
                     } else {
                       l = list()
                       nb = length(methods::slot(e1,"beta"))
                       nd = length(methods::slot(e1,"delta"))
                       for (i in 1:e2) l[[i]] = initial(nb,nd,methods::slot(e1,"lambdaStar"),TRUE)
                     }
                     l
                   })
#' @export
methods::setMethod("*",methods::signature(e1 = "numeric", e2 = "bayesPO_initial"),function(e1,e2) e2*e1)

#' @export
methods::setMethod("show","bayesPO_initial",function(object){
  cat("Initial values for a bayesPO model.")
  if (methods::slot(object,"tag") == "supplied") cat(" Values were supplied by user.\n\n")
  else if (methods::slot(object,"tag") == "random") cat(" Values were randomly generated.\n\n")

  cat("beta:\n")
  print(methods::slot(object,"beta"))
  cat("\ndelta:\n")
  print(methods::slot(object,"delta"))
  cat("\nlambdaStar:\n")
  print(methods::slot(object,"lambdaStar"))

  invisible(object)
})

#' @export
methods::setMethod("print","bayesPO_initial",function(x,...) show(x))

#' @export
print.bayesPO_fit <- function(x,...) show(x)

## Initializing
## If random = FALSE, the other parameters set the initial values.
## If random = TRUE, beta and delta set the size of the vector and lambdaStar
## sets the random value average (a good guess is the size of the presence-only
## data set).
#' @export
initial <- function(beta=numeric(),delta=numeric(),lambdaStar=numeric(),random = FALSE){
  if (random) methods::new("bayesPO_initial",beta = stats::rnorm(beta), delta = stats::rnorm(delta), lambdaStar = stats::rgamma(1,lambdaStar,1),tag = "random")
  else methods::new('bayesPO_initial',beta=beta,delta=delta,lambdaStar=lambdaStar,tag="supplied")
}
