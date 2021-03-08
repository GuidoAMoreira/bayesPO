#' Class that defines a model for the bayesPO package.
#'
#' @slot po The matrix containing the covariates values for the data.
#' @slot intensityLink A string informing about the chosen link for the intensity covariates. Current acceptable choice is only "logit"
#' @slot intensitySelection A vector containing the indexes of the selected intensity columns in the \code{po} matrix.
#' @slot observabilityLink A vector containing the indexes of the selected observability columns in the \code{po} matrix.
#' @slot intensitySelection A vector containing either the names of the columns chosen for observability covariates or the number of the column indexes.
#' @slot init A list with objects of class \code{bayesPO_initial} indicating the initial values for each chain. The number of objects will indicate how many chains will be run.
#' @slot prior An object of class \code{bayesPO_prior} which indicates the joint prior distribution for the model parameters.
#' @slot iSelectedColumns If the intensity covariates selection was made with the name of the columns, they are stored in this slot.
#' @slot iSelectedColumns If the observability covariates selection was made with the name of the columns, they are stored in this slot.
methods::setClass("bayesPO_model",
                  methods::representation(po="matrix",
                        intensityLink="character",
                        intensitySelection="numeric",
                        observabilityLink="character",
                        observabilitySelection="numeric",
                        init = "list",
                        prior = "bayesPO_prior",
                        iSelectedColumns = "character",
                        oSelectedColumns = "character"
                        ),
         validity = function(object){
           if (length(methods::slot(object,"intensityLink"))>1) stop("Argument intensityLink must have length 1.")
           if (length(methods::slot(object,"observabilityLink"))>1) stop("Argument observabilityLink must have length 1.")
           validLinks = c("logit")
           if (!(methods::slot(object,"intensityLink") %in% validLinks)) stop(paste0(methods::slot(object,"intensityLink")," is not a valid link. Accepted options for current version are\n",validLinks))
           if (!(methods::slot(object,"observabilityLink") %in% validLinks)) stop(paste0(methods::slot(object,"observabilityLink")," is not a valid link. Accepted options for current version are\n",validLinks))
           if (any(methods::slot(object,"intensitySelection") != as.integer(methods::slot(object,"intensitySelection")))) stop("intensitySelection must be integers.")
           if (any(methods::slot(object,"observabilitySelection") != as.integer(methods::slot(object,"observabilitySelection")))) stop("observabilitySelection must be integers.")
           nbSel = length(methods::slot(object,"intensitySelection"))+1; ndSel = length(methods::slot(object,"observabilitySelection"))+1
           for (i in methods::slot(object,"init")){
             if (!is(i,"bayesPO_initial")) stop("Initial values must be constructed with the initial function.")
             if (length(methods::slot(i,"beta")) != nbSel) stop(paste("\nInitial values for beta has the wrong size. Expected size:",nbSel))
             if (length(methods::slot(i,"delta")) != ndSel) stop(paste("\nInitial values for delta has the wrong size. Expected size:",ndSel))
           }
           if (length(methods::slot(methods::slot(methods::slot(object,"prior"),"beta"),"mu")) != nbSel) stop(paste("Prior for beta has wrong sized parameters. Expected size:",nbSel))
           if (length(methods::slot(methods::slot(methods::slot(object,"prior"),"delta"),"mu")) != ndSel) stop(paste("Prior for delta has wrong sized parameters. Expected size:",ndSel))
           TRUE
         })

#' @export
methods::setMethod("initialize","bayesPO_model",function(.Object,po,intensityLink,intensitySelection,
                                                observabilityLink,observabilitySelection,
                                                init,prior,iSelectedColumns,oSelectedColumns){
  cat("Loading data with",nrow(po),"observed points.\n")
  methods::slot(.Object,"po") <- po
  methods::slot(.Object,"intensityLink") <- intensityLink; methods::slot(.Object,"observabilityLink") = observabilityLink
  methods::slot(.Object,"intensitySelection") <- intensitySelection; methods::slot(.Object,"observabilitySelection") = observabilitySelection
  methods::slot(.Object,"init") <- init
  methods::slot(.Object,"prior") <- prior
  methods::slot(.Object,"iSelectedColumns") <- iSelectedColumns
  methods::slot(.Object,"oSelectedColumns") <- oSelectedColumns
  methods::validObject(.Object)
  cat("Data loaded successfully with ",length(intensitySelection)," intensity variables and ",
      length(observabilitySelection)," observability variables selected.\n",length(init)," chains ",
      ifelse(length(init)>1,"were","was")," initialized.\n",sep="")
  if (!length(iSelectedColumns)) cat("Intensity covariates selected with column indexes. Make sure the background covariates are in the same position.\n")
  if (!length(oSelectedColumns)) cat("Observability covariates selected with column indexes. Make sure the background covariates are in the same position.\n")
  .Object
})

#' @export
methods::setMethod("names","bayesPO_model",
          function(x) c("po","intensityLink","intensitySelection",
                        "observabilityLink","observabilitySelection",
                        "initial values","prior"))

#' @export
methods::setMethod("$","bayesPO_model",function(x,name){
  if (name %in% c("initial values","initial")) methods::slot(x,"init") else methods::slot(x,name)
})

#' @export
methods::setMethod("$<-","bayesPO_model",function(x,name,value){
  if (name %in% c("initial values","initial")) methods::slot(x,"init") = value else methods::slot(x,name) = value
  methods::validObject(x)
  x
})


