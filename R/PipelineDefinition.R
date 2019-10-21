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
#' @describeIn PipelineDefinition
#' @export
setClass( "PipelineDefinition", 
          slots=representation( functions="list", descriptions="list", 
                                evaluation="list", aggregation="list", 
                                defaultArguments="list", misc="list" ),
          prototype=prototype( functions=list(), descriptions=list(), 
                               evaluation=list(), aggregation=list(), 
                               defaultArguments=list(), misc=list() ),
          validity=.validatePipelineDef )



#' PipelineDefinition
#'
#' @param functions A list of functions for each step
#' @param descriptions A list of descriptions for each step
#' @param evaluation A list of optional evaluation functions for each step
#' @param aggregation A list of optional aggregation functions for each step
#' @param defaultArguments A lsit of optional default arguments
#' @param misc A list of whatever.
#' @param verbose Whether to output additional warnings (default TRUE).
#'
#' @return An object of class `PipelineDefinition`, with the slots functions,
#' descriptions, evaluation, aggregation, defaultArguments, and misc.
#' @export
PipelineDefinition <- function( functions, descriptions=NULL, evaluation=NULL,
                                aggregation=NULL, defaultArguments=list(), 
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
	         evaluation=evaluation, aggregation=aggregation, 
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

#' @export
setMethod("show", signature("PipelineDefinition"), function(object){
  blue("test")
  fns <- sapply(names(object@functions), FUN=function(x){ 
    y <- paste0("  - \033[1m",x,"\033[22m(",paste(names(formals(object@functions[[x]])),collapse=", "),")")
    if(!is.null(object@evaluation[[x]]) || !is.null(object@aggregation[[x]])) 
      y <- paste0(y," \033[34m*\033[39m ")
    if(!is.null(object@descriptions[[x]])) y <- paste(y,paste0("\033[3m",object@descriptions[[x]],"\033[23m"), sep="\n")
    y
  })
  cat("A PipelineDefinition object with the following steps:\n")
  cat(paste(fns,collapse="\n"))
  cat("\n")
})

#' @export
setMethod("names", signature("PipelineDefinition"), function(x){
  names(x@functions)
})

#' @export
setMethod("$", signature("PipelineDefinition"), function(x, name){
  x@functions[[name]]
})

#' @export
setMethod("length", signature("PipelineDefinition"), function(x){
  length(x@functions)
})

#' @export
setMethod("[",signature("PipelineDefinition"), function(x, i){
  new("PipelineDefinition", functions=x@functions[i], 
       descriptions=x@descriptions[i], evaluation=x@evaluation[i],
       aggregation=x@aggregation[i], misc=x@misc)
})

#' @export
setMethod("as.list",signature("PipelineDefinition"), function(x){
  x@functions
})

#' @importFrom base args
setGeneric("arguments", function(name) args(name))
#' @export
setMethod("arguments",signature("PipelineDefinition"), function(name){
  lapply(name@functions, FUN=function(x){ setdiff(names(formals(x)), "x") })
})

#' @export
setMethod("c",signature("PipelineDefinition"), function(x, ...){
  stop("Not yet implemented")
  ## TO DEBUG
  y <- unlist(list(...))
  ycl <- sapply(y, class)
  if( !all(ycl=="PipelineDefinition") ){
    if(!all( ycl %in% c("PipelineDefinition","function"))){
      stop("A PipelineDefinition can only be concatenated with other
           PipelineDefinition or functions.")
    }
    if(is.null(names(y)) || (any(names(y)=="" & ycl=="function")))
      stop("The function arguments should be named.")
    y <- lapply(names(y), FUN=function(x){
      if(is.function(y[[x]])) y[[x]] <- PipelineDefinition(y[x])
    })
  }
  y <- c(x,y)
  if(any(table(unlist(lapply(y,names)))>1))
    stop("Some step names are not unique.")
  new( "PipelineDefinition",
       functions=do.call(c, lapply(y, FUN=function(x) x@functions)),
       descriptions=do.call(c, lapply(y, FUN=function(x) x@functions)),
       aggregation=do.call(c, lapply(y, FUN=function(x) x@aggregation)),
       misc=do.call(c, lapply(y, FUN=function(x) x@misc))
    )
})
