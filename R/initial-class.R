#' Class for the initial values for the MCMC for the bayesPO package
#'
#' @slot beta Initial values for beta.
#' @slot delta Initial values for delta.
#' @slot lambdaStar Initial values for lambdaStar.
#' @slot tag Indicates the source of the initial values.
#' @export
#' @exportClass bayesPO_initial
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
#' @exportMethod initialize
methods::setMethod("initialize", "bayesPO_initial",
                   function(.Object, beta, delta, lambdaStar, tag){
  methods::slot(.Object,"beta") <- as.numeric(beta)
  methods::slot(.Object,"delta") <- as.numeric(delta)
  methods::slot(.Object,"lambdaStar") <- as.numeric(lambdaStar)
  methods::slot(.Object,"tag") <- tag
  .Object
})

#' @name bayesPO_initial-class
#' @param x The bayesPO_initial object.
#' @export
#' @exportMethod names
methods::setMethod("names","bayesPO_initial", function(x) c("beta", "delta", "lambdaStar"))

#' @name bayesPO_initial-class
#' @param x The bayesPO_initial object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$","bayesPO_initial",function(x,name) methods::slot(x, name))

#' @name bayesPO_initial-class
#' @param e1 A bayesPO_initial object.
#' @param e2 Another bayesPO_initial object or a list with bayesPO_initial
#' objects for \strong{+} and a positive integer for \strong{*}.
#' @return \strong{\code{+}}: A list with the objects. Useful to start the
#' fit function.
#' @export
#' @exportMethod +
methods::setMethod("+", "bayesPO_initial", function(e1, e2) list(e1, e2))

#' @name bayesPO_initial-class
#' @param e1 A bayesPO_initial object.
#' @export
#' @exportMethod +
methods::setMethod("+",methods::signature(e1 = "bayesPO_initial", e2 = "list"),
                   function(e1, e2){
                     for (i in 1:length(e2))
                      if (!is(e2[[i]], "bayesPO_initial"))
                        stop("Initial values can only be added to a list of other initial values.")
                     e2[[i + 1]] = e1
                     e2
                   })

#' @export
#' @exportMethod +
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

#' @name bayesPO_initial-class
#' @param e1 A bayesPO_initial object.
#' @return \strong{\code{*}}: A list with \code{e2} random initial values.
#' @export
#' @exportMethod *
methods::setMethod("*", methods::signature(e1 = "bayesPO_initial", e2 = "numeric"),
                   function(e1, e2) {
                     if (e2 <= 0) stop("Can only multiply by a positive number.")
                     if (as.integer(e2) != e2) stop("Con only multiply by an integer.")
                     if (methods::slot(e1, "tag") != "random") {
                       message("Identical initial values are not recommended for independent Markov Chains.")
                       l = list()
                       for (i in 1:e2) l[[i]] = e1
                     } else {
                       l = list()
                       nb = length(methods::slot(e1, "beta"))
                       nd = length(methods::slot(e1, "delta"))
                       for (i in 1:e2) l[[i]] = initial(nb, nd, methods::slot(e1, "lambdaStar"), TRUE)
                     }
                     l
                   })

#' @export
#' @exportMethod *
methods::setMethod("*", methods::signature(e1 = "numeric", e2 = "bayesPO_initial"), function(e1, e2) e2 * e1)

#' @name bayesPO_initial-class
#' @param object A bayesPO_initial object.
#' @export
#' @exportMethod show
methods::setMethod("show", "bayesPO_initial", function(object){
  cat("Initial values for a bayesPO model.")
  if (methods::slot(object, "tag") == "supplied") cat(" Values were supplied by user.\n\n")
  else if (methods::slot(object, "tag") == "random") cat(" Values were randomly generated.\n\n")

  cat("beta:\n")
  print(methods::slot(object, "beta"))
  cat("\ndelta:\n")
  print(methods::slot(object, "delta"))
  cat("\nlambdaStar:\n")
  print(methods::slot(object, "lambdaStar"))

  invisible(object)
})

#' @export
#' @exportMethod print
methods::setMethod("print", "bayesPO_initial", function(x, ...) show(x))

#' @method print bayesPO_fit
#' @export
print.bayesPO_fit <- function(x, ...) show(x)

#' Initial values constructor for bayesPO modeling
#'
#' Helper function to create a valid set of initial values to be used with the
#' fit_bayesPO function.
#' @param beta Either a vector or a single integer. The vector is used if the
#' initial values are provided and the integer is used as the vector size to
#' be randomly generated.
#' @param delta Either a vector or a single integer. The vector is used if the
#' initial values are provided and the integer is used as the vector size to
#' be randomly generated.
#' @param lambdaStar A positive number.
#' @param random A logical value. If \code{TRUE}, then the initial values are
#' generated from standard normal distribution for \code{beta} and \code{delta}
#' and from a \code{Beta(lambdaStar, 1)} for \code{lambdaStar}. The latter is
#' generated as a low value due to potential explosive values resulting from
#' background area scaling.
#' @return A \code{bayesPO_initial} object.
#' @seealso \code{\link{bayesPO_initial-class}}.
#' @export
initial <- function(beta = numeric(), delta = numeric(), lambdaStar = numeric(),
                    random = FALSE){
  if (random) methods::new("bayesPO_initial", beta = stats::rnorm(beta),
                           delta = stats::rnorm(delta),
                           lambdaStar = stats::rbeta(1, lambdaStar, 1), # Simulated small for safety
                           tag = "random")
  else methods::new('bayesPO_initial', beta = beta, delta = delta,
                    lambdaStar = lambdaStar, tag = "supplied")
}
