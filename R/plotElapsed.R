#' plotElapsed
#' 
#' Plot total elapsed time per run, split per step.
#'
#' @param res Aggregated pipeline results
#' @param steps The step(s) to plot (default all)
#' @param agg.by The parameters by which to aggregate (set to FALSE to disable 
#' aggregation)
#' @param agg.fn Aggregation function
#' @param width Width of the bar; default 0.9, use 1 to remove the gaps
#' @param split.datasets Logical; whether to split the datasets into facets
#' @param return.df Logical; whether to return the data.frame instead of plot
#'
#' @return A ggplot, or a data.frame if `return.df=TRUE`
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
plotElapsed <- function(res, steps=names(res$elapsed$stepwise), agg.by, 
                        agg.fn=mean, width=0.9, split.datasets=TRUE,
                        return.df=FALSE ){
  allsteps <- names(res$elapsed$stepwise)
  steps <- match.arg(steps, allsteps, several.ok = TRUE)
  args <- unlist(arguments(metadata(res)$PipelineDefinition)[steps])
  if(isFALSE(agg.by)){
    agg.by <- args
  }else{
    agg.by <- match.arg(agg.by, args, several.ok = TRUE)
  }
  for(f in steps){
    colnames(res$elapsed$stepwise[[f]]) <- 
      gsub("^elapsed$",f, colnames(res$elapsed$stepwise[[f]]))
  }
  m <- res$elapsed$stepwise[[1]]
  for(f in names(res$elapsed$stepwise)[-1]){
    fields <- intersect(colnames(m),colnames(res$elapsed$stepwise[[f]]))
    m <- merge(m, res$elapsed$stepwise[[f]], by=fields)
  }
  by2 <- agg.by
  if(split.datasets) by2 <- unique(c("dataset", by2))
  m <- aggregate( m[,steps,drop=FALSE], 
                  by=m[,by2,drop=FALSE], FUN=mean)
  m3 <- reshape2::melt(m, id.vars=by2, variable.name="step")
  m3 <- m3[which(m3$step %in% steps),]
  m3$id <- apply(m3[,agg.by,drop=FALSE], 1, collapse=" ", FUN=paste)
  if(return.df) return(m3)
  p <- ggplot(m3, aes_string('id', 'value', fill='step')) + 
    geom_bar( stat="identity", width=width,
              position=position_stack(reverse=TRUE) ) + 
    coord_flip() + xlab("") + ylab("Time (s)")
  if(split.datasets) p <- p + facet_wrap(~dataset, scales="free_x")
  p
}
