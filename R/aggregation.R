#' readPipelineResults
#'
#' @param path The path (e.g. folder or prefix) to the results. Either `path` 
#' or `resfiles` should be given.
#' @param resfiles A vector of paths to `*.evaluation.rds` files. Either `path` 
#' or `resfiles` should be given.
#'
#' @return A list of results.
#' @examples
#' # we produce mock pipeline results:
#' pip <- mockPipeline()
#' datasets <- list( ds1=1:3, ds2=c(5,10,15) )
#' tmpdir1 <- paste0(tempdir(),'/')
#' res <- runPipeline(datasets, pipelineDef=pip, output.prefix=tmpdir1,
#'                    alternatives=list() )
#' # we read the evaluation files:
#' res <- readPipelineResults(tmpdir1)
#' @export
readPipelineResults <- function(path = NULL, resfiles = NULL) {
  if ( (!is.null(path) && !is.null(resfiles)) || 
       (is.null(path) && is.null(resfiles)) ) 
    stop("Exactly one of `path` or `resfiles` should be given.")
  if (!is.null(path)) {
    resfiles <- list.files(path, pattern="evaluation\\.rds", full.names=TRUE)
    if (length(resfiles) == 0) 
      resfiles <- list.files(dirname(path), full.names=TRUE, 
                             pattern=paste0(gsub("\\.","\\\\.",basename(path)),
                                            ".*\\.evaluation\\.rds$"))
    if (length(resfiles) == 0) 
      stop("Could not find evaluation files.")
  }
  ds <- vapply(strsplit(basename(resfiles), "\\."), FUN=function(x) rev(x)[3], 
               FUN.VALUE=character(1))
  if (any(duplicated(ds))) 
    stop("Some datasets appear to occur in more than one file. If these are ",
          "the result of multiple `runPipeline` calls, please read them ",
          "separately and use `mergePipelineResults()`.")
  names(resfiles) <- ds
  lapply(resfiles, readRDS)
}


#' aggregatePipelineResults
#'
#' Aggregates the evaluation and running times of `runPipeline` results. 
#' Results should be indicated either as a `path`` prefix or as a vector of 
#' paths to `evaluation\\.rds` files (`resfiles`).
#'
#' @param res A (named) list of results (per dataset), as produced by 
#' \code{\link{readPipelineResults}} (or `mergePipelineResults`).
#' @param pipDef An optional \code{\link{PipelineDefinition}} containing the 
#' aggregation methods. If omitted, that from the results will be used.
#'
#' @return A list with a slot for each step for which there is an aggregation 
#' method, or (if no aggregation method available) a list of the 
#' `stepIntermediateReturnObjects` of `runPipeline`
#'
#' @import S4Vectors
#' @export
#' @examples
#' # we produce mock pipeline results:
#' pip <- mockPipeline()
#' datasets <- list( ds1=1:3, ds2=c(5,10,15) )
#' tmpdir1 <- paste0(tempdir(),'/')
#' res <- runPipeline(datasets, pipelineDef=pip, output.prefix=tmpdir1,
#'                    alternatives=list() )
#' # we read the evaluation files:
#' res <- readPipelineResults(tmpdir1)
#' # we aggregate the results (equivalent to the output of `runPipeline`):
#' res <- aggregatePipelineResults(res)
aggregatePipelineResults <- function(res, pipDef = NULL) {
  if (length(res) == 2 && all(names(res) == c("evaluation", "elapsed"))) 
    stop("`res` appears to be already aggregated, or the results of a single", 
         " dataset")
  .checkRes(res, requirePDidentity = is.null(pipDef))
  if (is.null(pipDef)) 
    pipDef <- metadata(res[[1]])$PipelineDefinition
  if (is.null(pipDef)) 
    stop("No PipelineDefinition found!")
  
  # aggregating elapsed time
  names(steps) <- steps <- names(res[[1]]$elapsed$stepwise)
  reso <- SimpleList(
    evaluation=NULL, 
    elapsed=list(stepwise=lapply(steps, FUN=function(step){
        .aggElapsed(lapply(res, FUN = function(x) x$elapsed$stepwise[[step]]))
      }),
      total=.aggElapsed(lapply(res, FUN=function(x) x$elapsed$total)))
  )
  metadata(reso)$PipelineDefinition <- pipDef
  
  isn <- vapply(stepFn(pipDef, type = "aggregation"), is.null, logical(1))
  if (all(isn)) {
    warning("No aggregation defined in the pipelineDefinition; returning only",
    " running times.")
    return(reso)
  }
  
  names(isn) <- isn <- names(isn)[!isn]
  res <- lapply(res, FUN = function(x) x$evaluation)
  
  fullnames <- parsePipNames(names(res[[1]][[length(res[[1]])]]))
  reso$evaluation <- lapply(isn, FUN = function(x) {
    message("Aggregating evaluation results for: ", x)
    stepFn(pipDef,type="aggregation")[[x]](lapply(res, FUN=function(y) y[[x]]))
  })
  reso
}

.aggElapsed <- function(res) {
  el.tot <- parsePipNames(names(res[[1]]))
  el.tot <- el.tot[rep(seq_len(nrow(el.tot)), length(res)), , drop = FALSE]
  row.names(el.tot) <- NULL
  el.tot$dataset <- factor(rep(names(res), each = length(res[[1]])))
  el.tot$elapsed <- as.numeric(unlist(res))
  el.tot
}

.checkRes <- function(res1, res2 = NULL, requirePDidentity = TRUE) {
  res1 <- lapply(res1, FUN = function(x) x$evaluation)
  # check that all datasets have the same steps and number of analyses
  if (!all(table(unlist(lapply(res1, names))) == length(res1))) 
    stop("The different datasets were not produced with the same pipeline!")
  for (step in names(res1[[1]])) {
    tt <- table(unlist(lapply(res1, FUN = function(x) names(x[[step]]))))
    if (!all(tt == length(res1))) 
      stop("The different datasets do not have the same runs, i.e. they ",
          "include different sets of alternative parameters.")
  }
  if (!is.null(res2)) {
    .checkRes(res2)
    res2 <- lapply(res2, FUN = function(x) x$evaluation)
    if (names(res1)[[1]] != names(res2)[[1]]) 
      stop("The two sets of results were not produced with the same pipeline.")
    if (is(res1,"SimpleList") && !is.null(metadata(res1)$PipelineDefinition) && 
        is(res2,"SimpleList") && !is.null(metadata(res2)$PipelineDefinition)) {
      # we require that the evaluation functions be the same
      pd1 <- metadata(res1)$PipelineDefinition
      pd2 <- metadata(res2)$PipelineDefinition
      if (!identical( stepFn(pd1, type = "evaluation"),
                      stepFn(pd2, type = "evaluation") )) {
        msg <- paste("The evaluation functions of the PipelineDefinitions are", 
                     "not identical")
        if (requirePDidentity) 
          stop(msg)
        warning(msg)
      }
    }
  }
}

#' mergePipelineResults
#' 
#' Merges the (non-aggregated) results of any number of runs of `runPipeline`
#' using the same \code{\link{PipelineDefinition}} (but on different datasets 
#' and/or using different parameters). First read the different sets of results 
#' using \code{\link{readPipelineResults}}, and pass them to this function.
#' 
#'
#' @param ... Any number of lists of pipeline results, each as produced by 
#' \code{\link{readPipelineResults}}
#' @param rr Alternatively, the pipeline results can be passed as a list (in 
#' which case `...` is ignored)
#' @param verbose Whether to print processing information
#'
#' @return A list of merged pipeline results.
#' @examples
#' # we produce 2 mock pipeline results:
#' pip <- mockPipeline()
#' datasets <- list( ds1=1:3, ds2=c(5,10,15) )
#' tmpdir1 <- paste0(tempdir(),'/')
#' res <- runPipeline(datasets, pipelineDef=pip, output.prefix=tmpdir1,
#'                    alternatives=list() )
#' alternatives <- list(meth1=c('log2','sqrt'), meth2='cumsum')
#' tmpdir2 <- paste0(tempdir(),'/')
#' res <- runPipeline(datasets, alternatives, pip, output.prefix=tmpdir2)
#' # we read the evaluation files:
#' res1 <- readPipelineResults(tmpdir1)
#' res2 <- readPipelineResults(tmpdir2)
#' # we merge them:
#' res <- mergePipelineResults(res1,res2)
#' # and we aggregate:
#' res <- aggregatePipelineResults(res)
#' @export
mergePipelineResults <- function(..., rr=NULL, verbose=TRUE){
  if(is.null(rr)) rr <- list(...)
  if(any(vapply(rr, FUN.VALUE=logical(1L), 
                FUN=function(x) is.null(x[[1]]$evaluation))))
    stop("Some of the objects appear not to be pipeline results.")
  if(is.null(names(rr))) names(rr) <- paste0("res",1:length(rr))
  for(ds in rr) .checkSameSteps(rr)
  .checkSameSteps(lapply(rr, FUN=function(x) x[[1]]))
  if(verbose) print(lapply(rr, FUN=function(res){
    vapply(res, FUN.VALUE=integer(1L), 
           FUN=function(x) length(rev(x$evaluation)[[1]]))}))
  
  # concatenate results
  names(ds) <- ds <- unique(unlist(lapply(rr, names)))
  hasdup <- FALSE
  res <- lapply(ds, FUN=function(dataset){
    r2 <- lapply(rr, FUN=function(x) x[[dataset]])
    r2 <- r2[!vapply(r2,is.null,logical(1L))]
    if(length(r2)==1) return(r2[[1]])
    if(!hasdup){
      rn <- unlist(lapply(r2, FUN=function(x) names(rev(x$evaluation)[[1]])))
      if(any(duplicated(rn))) hasdup <- TRUE
    }
    res <- r2[[1]]
    res$errors <- NULL
    for(i in seq(2,length(r2))){
      x <- setdiff(names(r2[[i]]$elapsed$total), names(res$elapsed$total))
      res$elapsed$total <- c(res$elapsed$total, r2[[i]]$elapsed$total[x])
      for(step in names(res$evaluation)){
        x <- setdiff(names(r2[[i]]$evaluation[[step]]), 
                     names(res$evaluation[[step]]))
        res$evaluation[[step]] <- c(res$evaluation[[step]], 
                                    r2[[i]]$evaluation[[step]][x])
        res$elapsed$stepwise[[step]] <- c(res$elapsed$stepwise[[step]], 
                                          r2[[i]]$elapsed$stepwise[[step]][x])
      }
    }
    res
  })
  rm(rr)
  if(hasdup && verbose) 
    message("Some analyses were redundant across results sets, ",
            "and only the first given was kept.")
  
  # retain only complete parameter sets
  hasmissing <- FALSE
  for(step in names(res[[1]]$evaluation)){
    a <- table(unlist(lapply(res,FUN=function(x) names(x$evaluation[[step]]))))
    if(length(a)>0){
      if(!all(w <- a==length(res))) hasmissing <- TRUE
      a <- sort(names(a)[w])
      if(length(a)==0) stop("No complete parameter set!")
      res <- lapply(res, FUN=function(x){
        x$evaluation[[step]] <- x$evaluation[[step]][a]
        x$elapsed$stepwise[[step]] <- x$elapsed$stepwise[[step]][a]
        if(step==rev(names(x$evaluation))[1])
          x$elapsed$total <- x$elapsed$total[a]
        x
      })
    }
  }
  if(hasmissing && verbose) 
    warning("Some parameter combinations were not available for all datasets ",
            "and were dropped.")
  
  if(verbose) print(list(merged=vapply(res, FUN.VALUE=integer(1L), 
           FUN=function(x) length(rev(x$evaluation)[[1]]))))
  res
}

.checkSameSteps <- function(x, stopifnot=TRUE){
  steps <- table(unlist(lapply(x, FUN=function(x) names(x$evaluation))))
  if(all(steps==length(x))) return(TRUE)
  if(stopifnot) stop("Some datasets have different steps - were then run using",
                      " the same PipelineDefinition?")
  FALSE
}

.dsMergeResults <- function(res1, res2) {
  names(steps) <- steps <- names(res1$evaluation)
  esteps <- vapply(res1$evaluation, FUN = function(x) {
    !is.null(x) && length(x) > 0
  }, logical(1))
  edif <- setdiff(names(res2$elapsed$total), names(res1$elapsed$total))
  res <- SimpleList(evaluation = lapply(steps[esteps], FUN = function(step) {
    r1 <- res1$evaluation[[step]]
    r2 <- res2$evaluation[[step]]
    dif <- setdiff(names(r2), names(r1))
    c(r1, r2[dif])
  }), elapsed = list(stepwise = lapply(steps, FUN = function(step) {
    r1 <- res1$elapsed$stepwise[[step]]
    r2 <- res2$elapsed$stepwise[[step]]
    dif <- setdiff(names(r2), names(r1))
    c(r1, r2[dif])
  }), total = c(res1$elapsed$total, res2$elapsed$total[edif])))
  if (!is.null(metadata(res1)$PipelineDefinition)) {
    metadata(res)$PipelineDefinition <- metadata(res1)$PipelineDefinition
  } else if (!is.null(metadata(res2)$PipelineDefinition)) {
    metadata(res)$PipelineDefinition < metadata(res2)$PipelineDefinition
  }
  res
}

#' defaultStepAggregation
#'
#' @param x A list of results per dataset, each containing a list (1 element 
#' per combination of parameters) of evaluation metrics (coercible to vectors 
#' or matrix).
#'
#' @return A data.frame.
#' @export
#' @importFrom dplyr bind_rows
defaultStepAggregation <- function(x) {
  y <- x[[1]][[1]]
  if (is(y, "list") || is(y, "SimpleList")) {
    res <- lapply(seq_along(y), FUN = function(i) {
      defaultStepAggregation(lapply(x, FUN = function(x) {
        lapply(x, FUN = function(x) x[[i]])
      }))
    })
    names(res) <- names(y)
    return(res)
  }
  dplyr::bind_rows(lapply(x, FUN = function(x) {
    x <- cbind(parsePipNames(names(x)), do.call(rbind, x))
    row.names(x) <- NULL
    x
  }), .id = "dataset")
}
