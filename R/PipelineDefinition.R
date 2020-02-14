.validatePipelineDef <- function(object){
  e <- c()
  if(!is.list(object@functions) || !all(sapply(object@functions, is.function))) 
    e <- c("`functions` should be a (named) list of functions!")
  if(!all(sapply(object@functions, FUN=function(x) "x" %in% names(formals(x))))) 
    e <- c(e, "Each function should at least take the argument `x`.")
  isf <- function(x) is.null(x) || is.function(x)
  if(!is.list(object@aggregation) || !all(sapply(object@aggregation, isf))) 
    stop("`aggregation` should be a list of functions and/or NULL slots!")
  if(!is.list(object@evaluation) || !all(sapply(object@evaluation, isf))) 
    stop("`evaluation` should be a list of functions and/or NULL slots!")
  if(!all(names(object@descriptions)==names(object@functions))) 
    e <- c(e, "descriptions do not match functions.")
  if(!all(names(object@evaluation)==names(object@functions))) 
    e <- c(e, "evaluation do not match functions.")
  if(!all(names(object@aggregation)==names(object@functions))) 
    e <- c(e, "aggregation do not match functions.")
  args <- unlist( lapply( object@functions, 
                          FUN=function(x){ setdiff(names(formals(x)), "x") }) )
  if(any(duplicated(args))) e <- c(e, paste("Some arguments (beside `x`) is",
    "used in more than one step, which is not currently supported."))
  if(length( wa <- setdiff(names(object@defaultArguments),args) )>0)
    e <- c(e, paste("The following default arguments are not in the pipeline's 
           functions:", paste(wa, collapse=", ")))
  if(length(e) == 0) TRUE else e
}

#' @import methods
#' @exportClass PipelineDefinition
setClass( "PipelineDefinition", 
          slots=representation( functions="list", descriptions="list", 
                                evaluation="list", aggregation="list", 
                                initiation="function",
                                defaultArguments="list", misc="list" ),
          prototype=prototype( functions=list(), descriptions=list(), 
                               evaluation=list(), aggregation=list(), 
                               initiation=identity,
                               defaultArguments=list(), misc=list() ),
          validity=.validatePipelineDef )



#' PipelineDefinition
#' 
#' Creates on object of class `PipelineDefinition` containing step functions,
#' as well as optionally step evaluation and aggregation functions.
#'
#' @param functions A list of functions for each step
#' @param descriptions A list of descriptions for each step
#' @param evaluation A list of optional evaluation functions for each step
#' @param aggregation A list of optional aggregation functions for each step
#' @param initiation A function ran when initiating a dataset
#' @param defaultArguments A lsit of optional default arguments
#' @param misc A list of whatever.
#' @param verbose Whether to output additional warnings (default TRUE).
#'
#' @return An object of class `PipelineDefinition`, with the slots functions,
#' descriptions, evaluation, aggregation, defaultArguments, and misc.
#' 
#' @aliases PipelineDefinition-class
#' @seealso \code{\link{PipelineDefinition-methods}}, \code{\link{addPipelineStep}}.
#' For an example pipeline, see \code{\link{scrna_pipeline}}.
#' @export
PipelineDefinition <- function( functions, descriptions=NULL, evaluation=NULL,
                                aggregation=NULL, initiation=identity, 
                                defaultArguments=list(), 
                                misc=list(), verbose=TRUE ){
  if(!is.list(functions) || !all(sapply(functions, is.function))) 
    stop("`functions` should be a (named) list of functions!")
  n <- names(functions)
	if(is.null(n)) n <- names(functions) <- paste0("step",1:length(functions))
  if(!is.null(descriptions)){
    if(length(descriptions)!=length(functions)) 
      stop("`descriptions` should have the same length as `functions`")
    if( !is.null(names(descriptions)) ){
      if(!all(names(descriptions)==names(functions)) )
        stop("The names of `descriptions` should match those of `functions`")
    }
  }else{
    descriptions <- lapply(functions,FUN=function(x) NULL)
  }
  if(!is.null(evaluation)){
    if(length(evaluation)!=length(functions)) 
      stop("`evaluation` should have the same length as `functions`")
    if( !is.null(names(evaluation)) ){
      if(!all(names(evaluation)==names(functions)) )
        stop("The names of `evaluation` should match those of `functions`")
    }
  }else{
    evaluation <- lapply(functions,FUN=function(x) NULL)
  }
  if(!is.null(aggregation)){
    if(length(aggregation)!=length(functions)) 
      stop("`aggregation` should have the same length as `functions`")
    if( !is.null(names(aggregation)) ){
      if(!all(names(aggregation)==names(functions)) )
        stop("The names of `aggregation` should match those of `functions`")
    }
  }else{
    aggregation <- lapply(functions,FUN=function(x) NULL)
  }
  names(aggregation)<-names(evaluation)<-names(descriptions)<-names(functions)
  if(is.null(misc)) misc <- list()
	x <- new("PipelineDefinition", functions=functions, descriptions=descriptions,
	       evaluation=evaluation, aggregation=aggregation, initiation=initiation,
	       defaultArguments=defaultArguments, misc=misc)

	w <- which( !sapply(x@aggregation,is.null) & 
	              sapply(x@evaluation,is.null) )
	if(verbose && length(w)>0){
	  warning(paste("An aggregation is defined for some steps that do not have",
	                "a defined evaluation function: ",
	                paste(names(x@functions)[w], collapse=", "),
	                "It is possible that evaluation is performed by the step's",
	                "function itself.") )
	}
	x
}

#' Methods for \code{\link{PipelineDefinition}} class
#' @name PipelineDefinition-methods
#' @rdname PipelineDefinition-methods
#' @aliases PipelineDefinition-method
#' @seealso \code{\link{PipelineDefinition}}, \code{\link{addPipelineStep}}
#' @param object An object of class \code{\link{PipelineDefinition}}
NULL

#' @rdname PipelineDefinition-methods
#' @importMethodsFrom methods show
#' @importFrom knitr opts_current
setMethod("show", signature("PipelineDefinition"), function(object){
  # colors and bold are going to trigger errors when rendered in a knit, so
  # we disable them when rendering
  isKnit <- tryCatch( isTRUE(getOption('knitr.in.progress')) || 
                        length(knitr::opts_current$get())>0,
                      error=function(e) FALSE)
  fns <- sapply(names(object@functions), FUN=function(x){ 
    x2 <- x
    if(!isKnit) x2 <- paste0("\033[1m",x,"\033[22m")
    y <- sapply( names(formals(object@functions[[x]])), FUN=function(n){
      if(!is.null(def <- object@defaultArguments[[n]]))
        n <- paste0(n,"=",deparse(def,100,FALSE))
      n
    })
    y <- paste0("  - ", x2, "(", paste(y, collapse=", "), ")")
    if(!is.null(object@evaluation[[x]]) || !is.null(object@aggregation[[x]])) 
      y <- paste0(y, ifelse(isKnit, " * ", " \033[34m*\033[39m "))
    if(!is.null(object@descriptions[[x]])){
      x2 <- object@descriptions[[x]]
      if(!isKnit) x2 <- paste0("\033[3m",x2,"\033[23m")
      y <- paste(y, x2, sep="\n")
    }
    y
  })
  cat("A PipelineDefinition object with the following steps:\n")
  cat(paste(fns,collapse="\n"))
  cat("\n")
})

#' @rdname PipelineDefinition-methods
#' @param x An object of class \code{\link{PipelineDefinition}}
setMethod("names", signature("PipelineDefinition"), function(x){
  names(x@functions)
})
#' @rdname PipelineDefinition-methods
setMethod("names<-", signature("PipelineDefinition"), function(x, value){
  if(any(duplicated(value))) stop("Some step names are duplicated!")
  names(x@functions) <- value
  names(x@evaluation) <- value
  names(x@aggregation) <- value
  names(x@descriptions) <- value
  validObject(x)
  x
})

#' @rdname PipelineDefinition-methods
setMethod("$", signature("PipelineDefinition"), function(x, name){
  x@functions[[name]]
})

#' @rdname PipelineDefinition-methods
setMethod("length", signature("PipelineDefinition"), function(x){
  length(x@functions)
})

#' @rdname PipelineDefinition-methods
setMethod("[",signature("PipelineDefinition"), function(x, i){
  new("PipelineDefinition", functions=x@functions[i], 
       descriptions=x@descriptions[i], evaluation=x@evaluation[i],
       aggregation=x@aggregation[i], misc=x@misc)
})

#' @rdname PipelineDefinition-methods
setMethod("as.list",signature("PipelineDefinition"), function(x){
  x@functions
})

#' @exportMethod arguments
setGeneric("arguments", function(object) args(object))
#' @rdname PipelineDefinition-methods
setMethod("arguments",signature("PipelineDefinition"), function(object){
  lapply(object@functions, FUN=function(x){ setdiff(names(formals(x)), "x") })
})

#' @exportMethod stepFn
setGeneric("stepFn", function(object, step, type) standardGeneric("stepFn"))
#' @param step The name of the step for which to set or get the function
#' @param type The type of function to set/get (either `functions`, `evaluation`,
#' `aggregation`, or `descriptions` - will parse partial matches)
#' @rdname PipelineDefinition-methods
setMethod("stepFn", signature("PipelineDefinition"), function(object, step, type){
  type <- match.arg(type, c("functions","evaluation","aggregation","descriptions"))
  step <- match.arg(step, names(object))
  slot(object, type)[[step]]
})
#' @exportMethod stepFn<-
setGeneric("stepFn<-", function(object, step, type, value) standardGeneric("stepFn<-"))
#' @rdname PipelineDefinition-methods
setMethod("stepFn<-", signature("PipelineDefinition"), function(object, step, type, value){
  type <- match.arg(type, c("functions","evaluation","aggregation","descriptions"))
  step <- match.arg(step, names(object))
  slot(object, type)[[step]] <- value
  validObject(object)
  object
})


#' addPipelineStep
#' 
#' Add a step to an existing \code{\link{PipelineDefinition}}
#'
#' @param object A \code{\link{PipelineDefinition}}
#' @param name The name of the step to add
#' @param after The name of the step after which to add the new step. If NULL, will
#' add the step at the beginning of the pipeline.
#' @param slots A optional named list with slots to fill for that step (i.e. `functions`,
#' `evaluation`, `aggregation`, `descriptions` - will be parsed)
#'
#' @return A \code{\link{PipelineDefinition}}
#' @seealso \code{\link{PipelineDefinition}}, \code{\link{PipelineDefinition-methods}}
#' @importFrom methods is slot
#' @export
#'
#' @examples
#' pd <- scrna_pipeline()
#' pd
#' pd <- addPipelineStep(pd, name="newstep", after="filtering", 
#'                       slots=list(description="Step that does nothing..."))
#' pd
addPipelineStep <- function(object, name, after=NULL, slots=list()){
  if(!is(object, "PipelineDefinition")) stop("object should be a PipelineDefinition")
  if(name %in% names(object)) stop("There is already a step with that name!")
  if(!is.null(after) && !(after %in% names(object))) 
    stop("`after` should either be null or the name of a step.")
  n <- c("functions","evaluation","aggregation","descriptions")
  if(length(slots)>0) names(slots) <- sapply(names(slots), choices=n, FUN=match.arg)
  if(!all(names(slots) %in% n)) stop( paste("fns should be a function or a list of", 
    "functions with one or more of the following names:\n", paste(n,collapse=", ")) )
  
  if(is.null(after)){
    i1 <- vector("integer")
    i2 <- seq_along(names(object))
  }else{
    w <- which(names(object)==after)
    i1 <- 1:w
    i2 <- (w+1):length(object)
    if(w==length(object)) i2 <- vector("integer")
  }
  ll <- list(NULL)
  names(ll) <- name
  for(f in n) slot(object,f) <- c(slot(object,f)[i1], ll, slot(object,f)[i2])
  for(f in names(slots)) stepFn(object, name, f) <- slots[[f]]
  if(is.null(stepFn(object, name, "functions"))) 
    stepFn(object, name, "functions") <- identity

  validObject(object)
  object
}