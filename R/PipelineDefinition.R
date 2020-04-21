.validatePipelineDef <- function(object){
  e <- c()
  if(!is.list(object@functions) || !all(sapply(object@functions, is.function))) 
    e <- c("`functions` should be a (named) list of functions!")
  if(!all(sapply(object@functions,FUN=function(x) "x" %in% names(formals(x))))) 
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
#' @seealso \code{\link{PipelineDefinition-methods}}, 
#' \code{\link{addPipelineStep}}. For an example pipeline, see 
#' \code{\link{scrna_pipeline}}.
#' @export
#' @examples
#' PipelineDefinition(
#'   list( step1=function(x, meth1){ get(meth1)(x) },
#'         step2=function(x, meth2){ get(meth2)(x) } )
#' )
PipelineDefinition <- function( functions, descriptions=NULL, evaluation=NULL,
                                aggregation=NULL, initiation=identity, 
                                defaultArguments=list(), 
                                misc=list(), verbose=TRUE ){
  if(!is.list(functions) || !all(sapply(functions, is.function))) 
    stop("`functions` should be a (named) list of functions!")
  n <- names(functions)
  if(is.null(n)) 
    n <- names(functions) <- paste0("step",seq_len(length(functions)))
  descriptions <- .checkInputList(descriptions, functions, FALSE)
  evaluation <- .checkInputList(evaluation, functions)
  aggregation2 <- .checkInputList(aggregation, functions)
  names(aggregation2)<-names(evaluation)<-names(descriptions)<-names(functions)
  for(f in names(aggregation2)){
    if(is.null(aggregation2[[f]]) && !is.null(evaluation[[f]]) &&
       !(f %in% names(aggregation)))
      aggregation2[[f]] <- defaultStepAggregation
  }
  if(is.null(misc)) misc <- list()
  x <- new("PipelineDefinition", functions=functions,descriptions=descriptions,
         evaluation=evaluation, aggregation=aggregation2, 
         initiation=initiation, defaultArguments=defaultArguments, misc=misc)
  
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

.checkInputList <- function( x, fns, containsFns=TRUE, 
                             name=deparse(substitute(x)) ){
  name <- paste0("`",name,"`")
  if(!is.null(x)){
    if(length(x)!=length(fns)){
      if(is.null(names(x)))
        stop("If ", name, " does not have the same length as the number of ",
             "steps, its slots should be named.")
      if(length(unknown <- setdiff(names(x),names(fns)))>0)
        stop("Some elements of ",name," (",paste(unknown,collapse=", "),")",
             "are unknown.")
      x <- lapply(names(fns), FUN=function(f){
        if(is.null(x[[f]])) return(NULL)
        x[[f]]
      })
      names(x) <- names(fns)
    } 
    if( !is.null(names(x)) ){
      if(!all(names(x)==names(fns)) )
        stop("The names of ",name," should match those of `functions`")
    }
  }else{
    x <- lapply(fns,FUN=function(x) NULL)
  }
  if( containsFns && 
      !all(sapply(x, FUN=function(x) is.null(x) || is.function(x))) )
    stop(name," should be a list of functions")
  x
}

#' Methods for \code{\link{PipelineDefinition}} class
#' @name PipelineDefinition-methods
#' @rdname PipelineDefinition-methods
#' @aliases PipelineDefinition-method
#' @seealso \code{\link{PipelineDefinition}}, \code{\link{addPipelineStep}}
#' @param object An object of class \code{\link{PipelineDefinition}}
#' @return Depends on the method.
#' @examples
#' pd <- mockPipeline()
#' length(pd)
#' names(pd)
#' pd$step1
#' pd[2:1]
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

#' get names of PipelineDefinition steps
#' @rdname PipelineDefinition-methods
#' @param x An object of class \code{\link{PipelineDefinition}}
setMethod("names", signature("PipelineDefinition"), function(x){
  names(x@functions)
})

#' set names of PipelineDefinition steps
#' @rdname PipelineDefinition-methods
#' @param value Replacement values
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
#' @param name The step name
setMethod("$", signature("PipelineDefinition"), function(x, name){
  x@functions[[name]]
})

#' @rdname PipelineDefinition-methods
setMethod("length", signature("PipelineDefinition"), function(x){
  length(x@functions)
})

#' @rdname PipelineDefinition-methods
#' @param i The index(es) of the steps
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
#' @rdname PipelineDefinition-methods
setGeneric("arguments", function(object) args(object))
#' @rdname PipelineDefinition-methods
setMethod("arguments",signature("PipelineDefinition"), function(object){
  lapply(object@functions, FUN=function(x){ setdiff(names(formals(x)), "x") })
})

#' @rdname PipelineDefinition-methods
#' @exportMethod defaultArguments
setGeneric("defaultArguments", function(object) NULL)
#' @exportMethod defaultArguments<-
#' @rdname PipelineDefinition-methods
setGeneric("defaultArguments<-", function(object, value) NULL)
#' @rdname PipelineDefinition-methods
setMethod("defaultArguments",signature("PipelineDefinition"), function(object){
  object@defaultArguments
})
#' @rdname PipelineDefinition-methods
setMethod( "defaultArguments<-",signature("PipelineDefinition"), 
           function(object, value){
  object@defaultArguments <- value
  validObject(object)
  object
})

#' @exportMethod stepFn
#' @rdname PipelineDefinition-methods
setGeneric("stepFn", function(object, step, type) standardGeneric("stepFn"))
#' @param step The name of the step for which to set or get the function
#' @param type The type of function to set/get, either `functions`, 
#' `evaluation`, `aggregation`, `descriptions`, or `initiation` (will parse 
#' partial matches)
#' @rdname PipelineDefinition-methods
setMethod("stepFn", signature("PipelineDefinition"), 
          function(object, step, type){
  ft <- c("functions","evaluation","aggregation","descriptions","initiation")
  type <- match.arg( type, ft )
  step <- match.arg(step, names(object))
  slot(object, type)[[step]]
})
#' @exportMethod stepFn<-
#' @rdname PipelineDefinition-methods
setGeneric( "stepFn<-", 
            function(object, step, type, value) standardGeneric("stepFn<-") )
#' @rdname PipelineDefinition-methods
setMethod( "stepFn<-", signature("PipelineDefinition"), 
           function(object, step, type, value){
  ft <- c("functions","evaluation","aggregation","descriptions","initiation")
  type <- match.arg(type, ft)
  if(type!="descriptions" &&  !is.function(value)) 
    stop("Replacement value should be a function.")
  if(type=="initiation"){
    slot(object, type) <- value
  }else{
    step <- match.arg(step, names(object))
    slot(object, type)[[step]] <- value
  }
  if(type=="evaluation" && !is.null(value)){
    # also add the default aggregation:
    if(is.null(slot(object, "aggregation")[[step]]))
      slot(object, "aggregation")[[step]] <- defaultStepAggregation
  }
  object
})


#' addPipelineStep
#' 
#' Add a step to an existing \code{\link{PipelineDefinition}}
#'
#' @param object A \code{\link{PipelineDefinition}}
#' @param name The name of the step to add
#' @param after The name of the step after which to add the new step. If NULL,
#' will add the step at the beginning of the pipeline.
#' @param slots A optional named list with slots to fill for that step (i.e. 
#' `functions`, `evaluation`, `aggregation`, `descriptions` - will be parsed)
#'
#' @return A \code{\link{PipelineDefinition}}
#' @seealso \code{\link{PipelineDefinition}}, 
#' \code{\link{PipelineDefinition-methods}}
#' @importFrom methods is slot
#' @export
#'
#' @examples
#' pd <- mockPipeline()
#' pd
#' pd <- addPipelineStep(pd, name="newstep", after="step1", 
#'                       slots=list(description="Step that does nothing..."))
#' pd
addPipelineStep <- function(object, name, after=NULL, slots=list()){
  if(!is(object, "PipelineDefinition")) 
    stop("object should be a PipelineDefinition")
  if(name %in% names(object)) stop("There is already a step with that name!")
  if(!is.null(after) && !(after %in% names(object))) 
    stop("`after` should either be null or the name of a step.")
  n <- c("functions","evaluation","aggregation","descriptions")
  if(length(slots)>0) 
    names(slots) <- sapply(names(slots), choices=n, FUN=match.arg)
  if(!all(names(slots) %in% n)) 
    stop(paste("fns should be a function or a list", 
    "with one or more of the following names:\n", paste(n,collapse=", ")))
  
  if(is.null(after)){
    i1 <- vector("integer")
    i2 <- seq_along(names(object))
  }else{
    w <- which(names(object)==after)
    i1 <- seq_len(w)
    i2 <- seq.int(from=w+1, to=length(object))
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

#' mockPipeline
#' 
#' A mock `PipelineDefinition` for use in examples.
#'
#' @return a `PipelineDefinition`
#' @export
#'
#' @examples
#' mockPipeline()
mockPipeline <- function(){
  PipelineDefinition(
    list( step1=function(x, meth1){ get(meth1)(x) },
          step2=function(x, meth2){ get(meth2)(x) } ),
    evaluation=list( step2=function(x) c(mean=mean(x), max=max(x)) ),
    descriptions=list( step1="This steps applies meth1 to x.",
                       step2="This steps applies meth2 to x."),
    defaultArguments=list(meth1=c("log","sqrt"), meth2="cumsum")
  )
}
