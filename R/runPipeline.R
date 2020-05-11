#' runPipeline
#' 
#' This function runs a pipeline with combinations of parameter variations on 
#' nested steps. The pipeline has to be defined as a list of functions applied 
#' consecutively on their respective outputs. See 'examples' for more details. 
#'
#' @param datasets A named vector of initial objects or paths to rds files.
#' @param alternatives The (named) list of alternative values for each 
#' parameter.
#' @param pipelineDef An object of class `PipelineDefinition`.
#' @param comb An optional matrix of indexes indicating the combination to run. 
#' Each column should correspond to an element of `alternatives`, and contain 
#' indexes relative to this element. If omitted, all combinations will be 
#' performed.
#' @param output.prefix An optional prefix for the output files.
#' @param nthreads Number of threads, default 1. If the memory requirements are
#' very high or the first steps very long to compute, consider setting this as
#' the number of datasets or below.
#' @param saveEndResults Logical; whether to save the output of the last step.
#' @param debug Logical (default FALSE). When enabled, disables multithreading 
#' and prints extra information.
#' @param ... passed to MulticoreParam. Can for instance be used to set 
#' `makeCluster` arguments, or set `threshold="TRACE"` when debugging in a 
#' multithreaded context.
#'
#' @examples 
#' pip <- mockPipeline()
#' datasets <- list( ds1=1:3, ds2=c(5,10,15) )
#' tmpdir1 <- paste0(tempdir(),"/")
#' res <- runPipeline(datasets, pipelineDef=pip, output.prefix=tmpdir1,
#'                    alternatives=list() )
#' # See the `pipeComp_scRNA` vignette for a more complex example
#' 
#' @return A SimpleList with elapsed time and the results of the evaluation 
#' functions defined by the given `pipelineDef`.
#' 
#' The results are also stored in the output folder with: 
#' \itemize{
#' \item The clustering results for each dataset (`endOutputs.rds` files),
#' \item A SimpletList of elapsed time and evaluations for each dataset 
#'  (`evaluation.rds` files),
#' \item A list of the `pipelineDef`, `alternatives`, `sessionInfo()` and
#'  function call used to produce the results (`runPipelineInfo.rds` file),
#' \item A copy of the SimpleList returned by the function 
#'  (`aggregated.rds`file). 
#' }
#' 
#' @importFrom utils sessionInfo
#' @import methods BiocParallel S4Vectors
#' @export
runPipeline <- function( datasets, alternatives, pipelineDef, comb=NULL, 
                         output.prefix="", nthreads=1, 
                         saveEndResults=TRUE, debug=FALSE, ...){
  mcall <- match.call()
  if(!is(pipelineDef,"PipelineDefinition")) 
    pipelineDef <- PipelineDefinition(pipelineDef)
  alternatives <- .checkPipArgs(alternatives, pipelineDef)
  pipDef <- pipelineDef@functions
  
  if(is.null(names(datasets)))
    names(datasets) <- paste0("dataset",seq_along(datasets))
  if(any(grepl(" ",names(datasets)))) 
    stop("Dataset names should not have spaces.")
  if(any(grepl("\\.",names(datasets)))) 
    warning("It is recommended not to use dots ('.') in dataset names to 
            facilitate browsing aggregated results.")
  
  # check that output folder exists, otherwise create it
  if(output.prefix!=""){
    x <- gsub("[^/]+$","",output.prefix)
    if(x!="" && !dir.exists(x)) dir.create(x, recursive=TRUE)
  }
  
  # prepare the combinations of parameters to use
  alt <- alternatives[unlist(arguments(pipelineDef))]
  if(is.null(comb)){
    eg <- buildCombMatrix(alt, TRUE)
  }else{
    eg <- .checkCombMatrix(comb, alt)
  }
  names(dsnames) <- dsnames <- names(datasets)
  
  if(!debug){
    if(any( class(nthreads) %in% 
            c("SnowParam","MulticoreParam","SerialParam") )){
      bpp <- nthreads
      nthreads <- bpnworkers(bpp)
    }else{
      bpp <- MulticoreParam(nthreads, ...)
    }
    message(paste("Running", nrow(eg), "pipeline settings on", length(datasets),
                  "datasets using", nthreads,"threads"))
  }else{
    nthreads <- 1
    bpp <- SerialParam()
  }
    
  if(!debug && nthreads>length(datasets)){
    # splitting downstream of datasets
    egs <- .splitEG(names(datasets), eg, nthreads)
    dsnames <- rep(names(datasets),each=length(egs))
    resfiles <- bpmapply( dsi=dsnames, eg=egs, BPPARAM=bpp,
                          FUN=function(dsi, eg) tryCatch(
                            .runPipelineF( dsi, datasets[[dsi]], pipelineDef, 
                                           alt, eg, output.prefix, noWrite=TRUE,
                                           saveEndResults=saveEndResults ),
                              error=function(e){
                                print(e)
                                stop(e)
                              })
                          )
    resfiles <- split(resfiles, dsnames)
    # reassembling and saving the results per dataset
    resfiles <- lapply(names(datasets), FUN=function(dsi){
      if(saveEndResults){
        res <- do.call(c, lapply(resfiles[[dsi]], FUN=function(x) x$res))
        saveRDS(res, file=paste0(output.prefix,"res.",dsi,".endOutputs.rds"))
      }
      if(length(resfiles[[dsi]])==1){
        res <- resfiles[[dsi]][[1]][c("evaluation", "elapsed")]
      }else{
        res <- resfiles[[dsi]][[1]]
        for(i in seq.int(from=2, to=length(resfiles[[dsi]]))){
          res <- .dsMergeResults(res, resfiles[[dsi]][[i]])
        }
      }
      metadata(res)$PipelineDefinition <- pipelineDef
      ifile <- paste0(output.prefix,"res.",dsi,".evaluation.rds")
      saveRDS(res, file=ifile)
      return(ifile)
    })
  }else if(!debug && nthreads>1 && length(datasets)>1){
    resfiles <- bplapply( dsnames, BPPARAM=bpp,
                          FUN=function(dsi){
                            tryCatch(
                    .runPipelineF(dsi, datasets[[dsi]], pipelineDef, alt, eg,
                                  output.prefix, saveEndResults=saveEndResults),
                                     error=function(e){
                                       print(e)
                                       print(traceback())
                                       stop(e)
                                     })
                          })
  }else{
    nthreads <- 1
    if(debug) message("Running in debug mode (single thread)")
    resfiles <- lapply( dsnames, FUN=function(dsi){
      .runPipelineF(dsi, datasets[[dsi]], pipelineDef, alt, eg, output.prefix, 
                    debug, saveEndResults=saveEndResults)
    })
  }

  message("
                Finished running on all datasets, now aggregating results...")
  
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
  saveRDS(pipinfo, file=paste0(output.prefix,"runPipelineInfo.rds"))
  
  names(resfiles) <- names(datasets)
  res <- lapply(resfiles, readRDS)
  res <- aggregatePipelineResults(res, pipelineDef)
  saveRDS(res, file=paste0(output.prefix,"aggregated.rds"))
  
  res
}


.runPipelineF <- function( dsi, ds, pipelineDef, alt, eg, output.prefix,
                           debug=FALSE, saveEndResults=FALSE, noWrite=FALSE ){
  
  if(debug) message(dsi)
  
  pipDef <- pipelineDef@functions
  args <- arguments(pipelineDef)
  
  ds <- tryCatch( pipelineDef@initiation(ds),
                  error=function(e){
                    stop("Error trying to initiate dataset ", dsi,"
",e)                      
                  })
  
  elapsed <- lapply(pipDef, FUN=function(x) list())
  elapsed.total <- list()
  
  objects <- c(list(BASE=ds),lapply(args[-length(args)],FUN=function(x)NULL))
  intermediate_return_objects <- lapply(args, FUN=function(x) list() )
  rm(ds)
  
  res <- vector(mode="list", length=nrow(eg))
  for(n in seq_len(nrow(eg))){
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
    wStep <- which(vapply(args,FUN=function(x){ chParam %in% x },logical(1)))
    # fetch the object from the previous step)
    while( is.null(x=objects[[which(names(objects)==names(args)[wStep])-1]]) &&
           wStep > 2 ) wStep <- wStep-1 # to handle steps without parameter
    if(is.null(x=objects[[which(names(objects)==names(args)[wStep])-1]])){
      # going back to original dataset
      x <- objects[[1]]
      wStep <- 1
    }else{
      x <- objects[[which(names(objects)==names(args)[wStep])-1]]
    }
    # proceed with the remaining steps of the pipeline
    for(step in names(args)[seq.int(from=wStep, to=length(args))]){
      if(debug) message(step)
      # prepare the arguments
      a <- lapply(args[[step]], FUN=function(a){
        alt[[a]][newPar[which(colnames(eg)==a)]]
      })
      names(a) <- args[[step]]
      #a$x <- x
      #x <- do.call(pipDef[[step]], a)   ## unknown issue with do.call...
      fcall <- .mycall(pipDef[[step]], a)
      if(debug) message(fcall)
      st <- Sys.time()
      x <- tryCatch( eval(fcall),
                     error=function(e){
                       .pipError(x, e, step, pipDef, fcall, newPar, 
                                 output.prefix, dsi, aa, debug)
                     })
      
      # name the current results on the basis of the previous steps:
      ws <- seq_len(sum( vapply(args[seq_len(which(names(args)==step))], 
                                length, integer(1)) ))
      ename <- .args2name(newPar[ws], alt[ws])
      # save elapsed time for this step
      elapsed[[step]][[ename]] <- as.numeric(Sys.time()-st)
      # save eventual intermediate objects (e.g. step evaluation)
      if(is.list(x) && all(c("x","intermediate_return") %in% names(x))){
        intermediate_return_objects[[step]][[ename]] <- x$intermediate_return
        x <- x$x
      }else{
        if(!is.null(pipelineDef@evaluation[[step]])){
          intermediate_return_objects[[step]][[ename]] <- tryCatch(
            pipelineDef@evaluation[[step]](x),
            error=function(e){
              fcall <- "pipelineDef@evaluation[[step]](x)"
              .pipError(x, e, step, pipDef, fcall, newPar, output.prefix, dsi,
                        aa, debug)
            })
        }
      }
      objects[[step]] <- x
    }
    
    # compute total time for this iteration
    elapsed.total[[n]] <- sum(vapply(names(args),FUN=function(step){
      ws <- seq_len(sum( vapply(args[seq_len(which(names(args)==step))], 
                                length, integer(1)) ))
      ename <- .args2name(newPar[ws], alt[ws])
      elapsed[[step]][[ename]]
    }, numeric(1) ))
    
    # return final results
    res[[n]] <- x
  }
  
  if(debug) message("
                      Completed running all variations.")
  
  # set names as the combination of arguments
  names(res) <- apply(eg,1,alt=alt,FUN=.args2name)
  names(elapsed.total) <- names(res)
  
  if(noWrite){
    out <- SimpleList( evaluation=intermediate_return_objects,
                       elapsed=list( stepwise=elapsed, total=elapsed.total ) )
    if(saveEndResults) out$res <- res
    return(out)
  }
  if(saveEndResults)
    saveRDS(res, file=paste0(output.prefix,"res.",dsi,".endOutputs.rds"))
  
  res <- SimpleList( evaluation=intermediate_return_objects,
                     elapsed=list( stepwise=elapsed, total=elapsed.total ) )
  metadata(res)$PipelineDefinition <- pipelineDef
  
  ifile <- paste0(output.prefix,"res.",dsi,".evaluation.rds")
  saveRDS(res, file=ifile)
  return(ifile)
}



# build function call (for a step of runPipeline) from list of arguments
.mycall <- function(fn, args){
  if(is.function(fn)) fn <- deparse(match.call()$fn)
  args <- paste(paste(names(args), vapply(args, FUN=function(x){
    if(is.numeric(x)) return(as.character(x))
    paste0("\"",x,"\"")
  }, character(1)), sep="="), collapse=", ")
  if(args!="") args <- paste(",",args)
  parse(text=paste0(fn, "(x=x", args, ")"))
}

.checkPipArgs <- function(alternatives, pipDef=NULL){
  if(any(grepl(";|=",names(alternatives)))) 
    stop("Some of the pipeline arguments contain unaccepted characters ",
      "(e.g. ';' or '=').")
  if(any(vapply(alternatives, function(x) any(grepl(";|=",x)), logical(1))))
    stop("Some of the alternative argument values contain unaccepted ",
      "characters (e.g. ';' or '=').")
  if(!is.null(pipDef)){
    def <- pipDef@defaultArguments
    for(f in names(alternatives)) def[[f]] <- alternatives[[f]]
    args <- arguments(pipDef)
    if(!all(unlist(args) %in% names(def))){
      missingParams <- setdiff(as.character(unlist(args)), names(def))
      stop("`alternatives` should have entries for the following slots defined",
        " in the pipeline: ", paste(missingParams ,collapse=", "))
    }
    if(!all( vapply(def, length, integer(1))>0)){
      stop("All steps of `alternatives` should contain at least one option.")
    }
    alternatives <- def
  }
  alternatives
}


.args2name <- function(x, alt){
  x2 <- mapply(a=alt,i=as.numeric(x),FUN=function(a,i) a[i])
  paste( paste0( names(alt), "=", x2), collapse=";" )
}

.pipError <- function(x, e, step, pipDef, fcall, newPar, output.prefix, dsi, aa, 
                      debug=TRUE){
  if(debug) save(x, step, pipDef, fcall, newPar, 
                 file=paste0(output.prefix,"runPipeline_error_TMPdump.RData"))
  msg <- paste0("Error in dataset `", dsi, "` with parameters:\n", aa, 
                "\nin step `", step, "`, evaluating command:\n`", fcall, "`
                Error:\n", e, "\n", 
                ifelse(debug, paste0("Current variables dumped in ", 
                                     output.prefix,
                                     "runPipeline_error_TMPdump.RData"), ""))
  if(!debug) print(msg)
  stop( msg )  
}

.splitEG <- function(datasets, eg, nthreads, tol=1){
  # we find at which step to split to maximize cluster use
  eg <- as.data.frame(eg)
  for(i in 1:ncol(eg)){
    if(length(unique(eg[,i]))*length(datasets) >= (nthreads-tol)){
      return(split(eg, eg[,seq_len(i)]))
    }
  }
  split(eg, seq_len(nrow(eg)))
}

