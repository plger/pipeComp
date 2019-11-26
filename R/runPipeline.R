#' runPipeline
#' 
#' This function runs a pipeline with combinations of parameter variations on nested steps.
#' The pipeline has to be defined as a list of functions applied consecutively on their 
#' respective outputs. See `defaultPipelineDef()` for indications about how to build one.
#'
#' @param datasets A named vector of initial objects or paths to rds files.
#' @param alternatives The (named) list of alternative values for each parameters.
#' @param pipelineDef An object of class `PipelineDefinition`.
#' @param eg An optional matrix of indexes indicating the combination to run. 
#' Each column should correspond to an element of `alternatives`, and contain 
#' indexes relative to this element. If omitted, all combinations will be 
#' performed.
#' @param output.prefix An optional prefix for the output files.
#' @param nthreads Number of threads, default 6.
#' @param saveEndResults Logical; whether to save the output of the last step.
#' @param debug Logical (default FALSE). When enabled, disables multithreading 
#' and prints extra information.
#' @param ... passed to MulticoreParam. Can for instance be used to set `makeCluster` 
#' arguments, or set `threshold="TRACE"` when debugging in a multithreaded context.
#'
#' @return 
#' 
#' @import methods BiocParallel data.table
#' @export
runPipeline <- function( datasets, alternatives, pipelineDef, eg=NULL, output.prefix="",
                         nthreads=6, saveEndResults=TRUE, debug=FALSE, ...){
  mcall <- match.call()
  if(!is(pipelineDef,"PipelineDefinition")) 
    pipelineDef <- PipelineDefinition(pipelineDef)
  pipDef <- pipelineDef@functions
  .checkPipArgs(alternatives, pipDef)
  
  if(is.null(names(datasets)))
    names(datasets) <- paste0("dataset",seq_along(datasets))
  if(any(grepl("\\.",names(datasets)))) 
    warning("It is recommended not to use dots ('.') in dataset names to 
            facilitate browsing aggregated results.")
  
  # extract the arguments required by the pipeline
  args <- arguments(pipelineDef)
  
  # check that output folder exists, otherwise create it
  if(output.prefix!=""){
    x <- gsub("[^/]+$","",output.prefix)
    if(x!="" && !dir.exists(x)) dir.create(x, recursive=TRUE)
  }
  
  # prepare the combinations of parameters to use
  alt <- alternatives[unlist(args)]
  if(is.null(eg)){
    eg <- data.table(expand.grid(lapply(alt, FUN=function(x){ 1:length(x) })))
  }else{
    if(!all(names(alt) %in% colnames(eg))) stop("The columns of `eg` do not correspond to the arguments.")
    eg <- eg[,names(alt)]
  }
  eg <- as.matrix(setorder(eg))
  
  ## BEGIN .runPipelineF
  .runPipelineF <- function(dsi){
    if(is.character(datasets[[dsi]])){
      ds <- readRDS(datasets[[dsi]])
    }else{
      ds <- datasets[[dsi]]
    }
    dsname <- names(datasets)[dsi]

    if(debug) message(dsname)

    elapsed <- lapply(pipDef, FUN=function(x) list())
    elapsed.total <- list()
    
    objects <- c(list(BASE=ds), lapply(args[-length(args)],FUN=function(x) NULL))
    intermediate_return_objects <- lapply(args, FUN=function(x) list() )
    rm(ds)
    
    res <- sapply(1:nrow(eg), FUN=function(x) NULL)
    for(n in 1:nrow(eg)){
      newPar <- as.numeric(eg[n,])
      aa <- paste( mapply(an=names(alt), a=alt, i=newPar, 
                          FUN=function(an,a,i) paste0(an,"=",a[i]) )
                   ,collapse=", ")
      if(debug){
        message("
                ####################
                # Iteration ",n,"
                # ", aa)
      }
      if(n==1){
        oldPar <- rep(0,ncol(eg))
      }else{
        oldPar <- as.numeric(eg[n-1,])
      }
      # identify the first parameters that changes and the corresponding step
      chParam <- names(alt)[which(newPar!=oldPar)[1]]
      wStep <- which(sapply(args,FUN=function(x){ chParam %in% x }))
      # fetch the object from the previous step
      x <- objects[[which(names(objects)==names(args)[wStep])-1]]
      # proceed with the remaining steps of the pipeline
      for(step in names(args)[wStep:length(args)]){
        if(debug) message(step)
        # prepare the arguments
        a <- lapply(args[[step]], FUN=function(a) alt[[a]][newPar[which(colnames(eg)==a)]] )
        names(a) <- args[[step]]
        #a$x <- x
        #x <- do.call(pipDef[[step]], a)   ## unknown issue with do.call... we use a custom eval function:
        fcall <- .mycall(pipDef[[step]], a)
        if(debug) message(fcall)
        st <- Sys.time()
        x <- tryCatch( eval(fcall),
                       error=function(e){
  ## error report
   if(debug) save(x, step, pipDef, fcall, newPar, 
                  file=paste0(output.prefix,"runPipeline_error_TMPdump.RData"))
   msg <- paste0("Error in dataset `", dsi, "` with parameters:\n", aa, 
                "\nin step `", step, "`, evaluating command:\n`", fcall, "`\nError:\n", 
                e, "\n", ifelse(debug, paste("Current variables dumped in", 
                                             paste0(output.prefix,"runPipeline_error_TMPdump.RData")), ""))
   if(!debug || nthreads>1) print(msg)
   stop( msg )
  ## end error report
                       })

        # name the current results on the basis of the previous steps:
        ws <- 1:sum(sapply(args[1:which(names(args)==step)], length))
        ename <- .args2name(newPar[ws], alt[ws])
        # save elapsed time for this step
        elapsed[[step]][[ename]] <- as.numeric(Sys.time()-st)
        # save eventual intermediate objects (e.g. step evaluation)
        if(is.list(x) && all(c("x","intermediate_return") %in% names(x))){
          intermediate_return_objects[[step]][[ename]] <- x$intermediate_return
          x <- x$x
        }else{
          if(!is.null(pipelineDef@evaluation[[step]])){
            intermediate_return_objects[[step]][[ename]] <- pipelineDef@evaluation[[step]](x)
          }
        }
        objects[[step]] <- x
      }

      # compute total time for this iteration
      elapsed.total[[n]] <- sum(sapply(names(args),FUN=function(step){
        ws <- 1:sum(sapply(args[1:which(names(args)==step)], length))
        ename <- .args2name(newPar[ws], alt[ws])
        elapsed[[step]][[ename]]
      }))
      
      # return final results
      res[[n]] <- x
    }
    
    if(debug) message("
                      Completed running all variations.")
    
    # set names as the combination of arguments
    names(res) <- apply(eg,1,alt=alt,FUN=.args2name)
    names(elapsed.total) <- names(res)
    res <- lapply(res,FUN=function(x){
      ## refactor keeping names
      xn <- names(x)
      x <- factor(as.character(x))
      names(x) <- xn
      x
    })
    
    if(saveEndResults){
      resfile <- paste0(output.prefix,"res.",dsi,".endOutputs.rds")
      saveRDS(res, file=resfile)
    }else{
      resfile <- NULL
    }
    
    res <- list( res=resfile, elapsed=elapsed, elapsed.total=elapsed.total )
    ifile <- paste0(output.prefix,"res.",dsi,".stepIntermediateReturnObjects.rds")
    if(sum(sapply(intermediate_return_objects,length))>0){
      saveRDS(intermediate_return_objects, file=ifile)
      res[["intermediate"]] <- ifile
    }
    res
  }
  ## END .runPipelineF
  
  names(dsnames) <- dsnames <- names(datasets)
  if(!debug && nthreads>1 && length(datasets)>1){
    nthreads <- min(nthreads, length(datasets))
    message(paste("Using",nthreads,"threads"))
    res <- bplapply( dsnames, 
                     BPPARAM=MulticoreParam(nthreads, ...), 
                     FUN=.runPipelineF )
  }else{
    nthreads <- 1
    if(debug) message("Running in debug mode (single thread)")
    res <- lapply( dsnames, FUN=.runPipelineF)
  }
  
  names(res) <- names(datasets)
  elapsed <- lapply(res, FUN=function(x) x$elapsed)
  names(steps) <- steps <- names(elapsed[[1]])
  elapsed <- lapply(steps, FUN=function(x){
    x <- do.call(cbind, lapply(elapsed, FUN=function(y) unlist(y[[x]])))
    colnames(x) <- names(datasets)
    x
  })
  elapsed.total <- do.call(cbind, lapply(res, FUN=function(x) unlist(x$elapsed.total)))
  colnames(elapsed.total) <- names(datasets)
  
  res <- lapply(res, FUN=function(x) x$res)
  
  # save pipeline and resolved functions
  pipinfo <- list( pipDef=pipelineDef,
                   alts=lapply(alt, FUN=function(x){ 
                     if(is.numeric(x)) return(x)
                     lapply(x,FUN=function(x){
                       if(is.function(x)) return(x)
                       if(exists(x) && is.function(get(x))){
                         return(get(x))
                       }else{
                         return(x)
                       }
                     })
                   }),
                   sessionInfo=sessionInfo(),
                   call=mcall
  )
  saveRDS(pipinfo, file=paste0(output.prefix,"pipelineInfo.rds"))
  
  # running time
  elapsed <- list( stepwise=elapsed, total=elapsed.total )
  saveRDS(elapsed, file=paste0(output.prefix,"elapsed.rds"))

  message("
                  Finished running on all datasets, now aggregating results...")
    
  aggregateResults(output.prefix, pipelineDef)
}


# build function call (for a step of runPipeline) from list of arguments
.mycall <- function(fn, args){
  if(is.function(fn)) fn <- deparse(match.call()$fn)
  args <- paste(paste(names(args), sapply(args, FUN=function(x){
    if(is.numeric(x)) return(x)
    paste0("\"",x,"\"")
  }), sep="="), collapse=", ")
  parse(text=paste0(fn, "(x=x, ", args, ")"))
}

.checkPipArgs <- function(alternatives, pipDef=NULL){
  if(any(grepl(";|=",names(alternatives)))) 
    stop("Some of the pipeline arguments contain unaccepted characters (e.g. ';' or '=').")
  if(any(sapply(alternatives, FUN=function(x) any(grepl(";|=",x)))))
    stop("Some of the alternative argument values contain unaccepted characters (e.g. ';' or '=').")
  if(!is.null(pipDef)){
    args <- lapply(pipDef, FUN=function(x){ setdiff(names(formals(x)), "x") })
    if(!all(unlist(args) %in% names(alternatives))){
      missingParams <- setdiff(as.character(unlist(args)), names(alternatives))
      stop("`alternatives` should have the following slots defined in the pipeline:",
           paste(as.character(unlist(args)),collapse=", "),"
          The following are missing:
          ", paste(missingParams ,collapse=", "))
    }
    if(!all( sapply(alternatives, FUN=length)>0)){
      stop("All steps of `alternatives` should contain at least one option.")
    }
  }
}


.args2name <- function(x, alt){
  x2 <- mapply(a=alt,i=as.numeric(x),FUN=function(a,i) a[i])
  paste( paste0( names(alt), "=", x2), collapse=";" )
}