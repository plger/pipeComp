#' aggregateResults
#'
#' Aggregates the `stepIntermediateReturnObjects` (and running times) of a 
#' `runPipeline` results.
#'
#' @param path Path to the folder containing the output of runPipeline
#' @param pipDef A `pipelineDefinition` containing the aggregation
#' methods.
#'
#' @return A list with a slot for each step for which there is an aggregation 
#' method, or (if no aggregation method available) a list of the 
#' `stepIntermediateReturnObjects` of `runPipeline`
#' @export
aggregateResults <- function(path="./", pipDef=NULL){
  elapsed <- list.files(path, "elapsed\\.rds$",full.names=TRUE)
  if(length(elapsed)==0 || length(elapsed)>1 && 
     file.exists(paste0(path,"elapsed.rds"))){
      elapsed <- paste0(path,"elapsed.rds")
      res <- list.files( dirname(path), full.names=TRUE, 
                         pattern=paste0(gsub("\\.","\\\\.",basename(path)),
                                ".*\\.stepIntermediateReturnObjects\\.rds$") )
      pi <- list.files( dirname(path), full.names=TRUE, 
                       pattern=paste0(gsub("\\.","\\\\.",basename(path)),
                                     "pipelineInfo\\.rds$") )
  }else{
    res <- list.files( path, full.names=TRUE, 
                       pattern="\\.stepIntermediateReturnObjects\\.rds$" )
    pi <- list.files( path, full.names=TRUE, pattern="pipelineInfo\\.rds$")
  }
  if(length(elapsed)==0) stop("`path` does not appear to contain the results of
a `runPipeline`")
  if(length(elapsed)>1) stop("`path` appears to contain the results of more than
one `runPipeline` run; please append a prefix to the `path`.")
  elapsed <- readRDS(elapsed)

  names(res) <- sapply( strsplit(basename(res),".",fixed=T), 
                        FUN=function(x){ x[length(x)-2] } )
  res <- lapply(res, readRDS)
  if(is.null(pipDef)) pipDef <- readRDS(pi)$pipDef
  if(is.null(pipDef)) return(list(res=res, elapsed=elapsed))
  
  isn <- sapply(pipDef@aggregation, is.null)
  if(all(isn)){
    warning("No aggregation defined in the pipelineDefinition; 
returning raw results and running times.")
    return(list(res=res, elapsed=elapsed))
  }
  names(isn) <- isn <- names(isn)[!isn]
  fullnames <- parsePipNames(names(res[[1]][[length(res[[1]])]]))
  res <- lapply(isn, FUN=function(x){
    message(x)
    res2 <- pipDef@aggregation[[x]](lapply(res, FUN=function(y) y[[x]]))
    if(!is(res2, "list")) res2 <- list(res2)
    if(!is.null(elapsed$stepwise[[x]])){
      el <- as.data.frame(elapsed$stepwise[[x]], row.names=NULL)
      colnames(el) <- paste0("stepElapsed.",colnames(el))
      res2 <- lapply(res2, FUN=function(x){
        suppressWarnings(cbind(x, el))
      })
    }
    if(length(res2)==1) res2 <- res2[[1]]
    res2
  })
  if(!is.null(elapsed$total)){
    et <- as.data.frame(elapsed$total)
    row.names(et) <- NULL
    colnames(et) <- paste0("totalTime.",colnames(et))
    res$elapsedTotal <- cbind( fullnames,  et)
  }
  res
}
