readPipelineResults <- function(path=NULL, resfiles=NULL){
  if( (!is.null(path) && !is.null(resfiles)) ||
      (is.null(path) && is.null(resfiles)) )
    stop("Exactly one of `path` or `resfiles` should be given.")
  if(!is.null(path)){
    resfiles <- list.files(path, pattern="evaluation\\.rds", full.names=TRUE)
    if(length(resfiles)==0)
      resfiles <- list.files(dirname(path), full.names=TRUE,
                             pattern=paste0(gsub("\\.","\\\\.",basename(path)),
                                            ".*\\.evaluation\\.rds$") )
    if(length(resfiles)==0) stop("Could not find evaluation files.")
  }
  ds <- sapply(strsplit(basename(resfiles),"\\."), FUN=function(x) rev(x)[3])
  if(any(duplicated(ds)))
    stop("Some datasets appear to occur in more than one file.",
         "If these are the result of multiple `runPipeline` calls, please read",
         "them separately and use `mergePipelineResults()`.")
  names(resfiles) <- ds
  lapply(resfiles, readRDS)
}


#' aggregatePipelineResults
#'
#' Aggregates the evaluation and running times of `runPipeline` results. Results
#'  should be indicated either as a `path`` prefix or as a vector of 
#' paths to `evaluation\\.rds` files (`resfiles`).
#'
#' @param res A (named) list of results (per dataset), as produced by 
#' `readPipelineResults` (or `mergePipelineResults`).
#' @param pipDef An optional `pipelineDefinition` containing the aggregation 
#' methods. If omitted, that from the results will be used.
#'
#' @return A list with a slot for each step for which there is an aggregation 
#' method, or (if no aggregation method available) a list of the 
#' `stepIntermediateReturnObjects` of `runPipeline`
#'
#' @import S4Vectors
#' @export
aggregatePipelineResults <- function(res, pipDef=NULL){
  .checkRes(res, requirePDidentity=is.null(pipDef))
  if(is.null(pipDef)) pipDef <- metadata(res[[1]])$PipelineDefinition
  if(is.null(pipDef)) stop("No PipelineDefinition found!")
    
  # aggregating elapsed time
  names(steps) <- steps <- names(res[[1]]$elapsed$stepwise)
  reso <- SimpleList( evaluation=NULL, elapsed=list(
    stepwise=lapply(steps, FUN=function(step){
      .aggElapsed(lapply(res, FUN=function(x) x$elapsed$stepwise[[step]]))
    }),
    total=.aggElapsed(lapply(res, FUN=function(x) x$elapsed$total))
  ))
  metadata(reso)$PipelineDefinition <- pipDef
  
  isn <- sapply(pipDef@aggregation, is.null)
  if(all(isn)){
    warning("No aggregation defined in the pipelineDefinition; 
returning only running times.")
    return(reso)
  }
  
  names(isn) <- isn <- names(isn)[!isn]
  res <- lapply(res, FUN=function(x) x$evaluation)
  
  fullnames <- parsePipNames(names(res[[1]][[length(res[[1]])]]))
  reso$evaluation <- lapply(isn, FUN=function(x){
    message(x)
    pipDef@aggregation[[x]](lapply(res, FUN=function(y) y[[x]]))
  })
  reso
}

.aggElapsed <- function(res){
  el.tot <- parsePipNames(names(res[[1]]))
  el.tot <- el.tot[rep(1:nrow(el.tot), length(res)),,drop=FALSE]
  row.names(el.tot) <- NULL
  el.tot$dataset <- factor(rep(names(res), each=length(res[[1]])))
  el.tot$elapsed <- as.numeric(unlist(res))
  el.tot
}

.checkRes <- function(res1, res2=NULL, requirePDidentity=TRUE){
  res1 <- lapply(res1, FUN=function(x) x$evaluation)
  # check that all datasets have the same steps and number of analyses
  if(!all(table(unlist(lapply(res1, names)))==length(res1))) 
    stop("The different datasets were not produced with the same pipeline!")
  for(step in names(res1[[1]])){
    tt <- table(unlist(lapply(res1, FUN=function(x) names(x[[step]]) )))
    if(!all(tt==length(res1))) stop("The different datasets do not have the same",
        " runs, i.e. they include different sets of alternative parameters.")
  }
  if(!is.null(res2)){
    .checkRes(res2)
    res2 <- lapply(res2, FUN=function(x) x$evaluation)
    if(names(res1)[[1]] != names(res2)[[1]])
      stop("The two sets of results were not produced with the same pipeline.")
    if(is(res1,"SimpleList") && !is.null(metadata(res1)$PipelineDefinition) &&
       is(res2,"SimpleList") && !is.null(metadata(res2)$PipelineDefinition) ){
      # we require that the evaluation functions be the same
      pd1 <- metadata(res1)$PipelineDefinition
      pd2 <- metadata(res2)$PipelineDefinition
      if(!identical(pd1@evaluation,pd2@evaluation)){
        msg <- "The evaluation functions of the PipelineDefinitions are not identical"
        if(requirePDidentity) stop(msg)
        warning(msg)
      }
    }
  }
}

#' mergePipelineResults
#'
#' @param res1 A list of pipeline results, as produced by `readPipelineResults`
#' @param res2 A list of pipeline results, as produced by `readPipelineResults`
#'
#' @return A list of pipeline results.
#' @export
#'
#' @import S4Vectors
mergePipelineResults <- function(res1,res2){
  .checkRes(res1,res2)
  nn <- intersect(names(res1),names(res2))
  if(length(nn) != length(res1)){
    if(length(nn)==0) stop("The results do not share datasets.")
    warning("The results will be reduced to the subset of datasets they share.")
    res1 <- res1[nn]
    res2 <- res2[nn]
  }
  names(nn) <- nn
  names(steps) <- steps <- names(res1[[1]]$evaluation)
  lapply(nn, FUN=function(ds){
    edif <- setdiff( names(res2[[ds]]$elapsed$total), 
                     names(res1[[ds]]$elapsed$total))
    res <- SimpleList(
      evaluation=lapply(steps, FUN=function(step){
        r1 <- res1[[ds]]$evaluation[[step]]
        r2 <- res2[[ds]]$evaluation[[step]]
        dif <- setdiff(names(r1), names(r2))
        c(r1, r2[dif])
      }),
      elapsed=list(
        stepwise=lapply(steps, FUN=function(step){
          r1 <- res1[[ds]]$elapsed$stepwise[[step]]
          r2 <- res2[[ds]]$elapsed$stepwise[[step]]
          dif <- setdiff(names(r1), names(r2))
          c(r1, r2[dif])
        }),
        total=c(res1[[ds]]$elapsed$total,res1[[ds]]$elapsed$total[edif])
      )
    )
    if(!is.null(res1$PipelineDefinition)) 
      metadata(res)$PipelineDefinition <- res1$PipelineDefinition
    res
  })
}
