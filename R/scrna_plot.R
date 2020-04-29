#' scrna_evalPlot_silh
#' 
#' Plot a min/max/mean/median silhouette width heatmap from aggregated 
#' evaluation results of the `scrna_pipeline`.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param what What metric to plot, possible values are “minSilWidth”, 
#' “meanSilWidth” (default), “medianSilWidth”, or “maxSilWidth”.
#' @param dims If multiple sets of dimensions are available, which one to use
#' (defaults to the first).
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param filterExpr An optional filtering expression based on the columns of 
#' the target dataframe, (e.g. `filterExpr=param1=="value1"`).
#' @param value_format Format for displaying cells' values (no label by default)
#' @param reorder_rows Whether to sort rows (default FALSE). The row 
#' names themselves can also be passed to specify an order, or a 
#' `ComplexHeatmap`.
#' @param reorder_columns Whether to sort columns (default TRUE).
#' @param row_split Optional column (included in `agg.by`) by which to split
#' the rows.
#' @param show_heatmap_legend Passed to `Heatmap` (default FALSE)
#' @param show_column_names Passed to `Heatmap` (default FALSE)
#' @param col Colors for the heatmap
#' @param font_factor A scaling factor applied to fontsizes (default 1)
#' @param value_cols A vector of length 2 indicating the colors of the values
#' (above and below the mean), if printed
#' @param title Plot title
#' @param shortNames Logical; whether to use short row names (with only
#' the parameter values instead of the parameter name and value pairs), default
#' TRUE.
#' @param anno_legend Logical; whether to plot the legend for the datasets
#' @param ... Passed to `Heatmap`
#' 
#' @return A Heatmap
#' @export
#' @import ComplexHeatmap
#' @examples
#' data("exampleResults", package="pipeComp")
#' scrna_evalPlot_silh( exampleResults, agg.by=c("filt","norm"), 
#'                      row_split="norm" )
scrna_evalPlot_silh <- function( res, what=c("minSilWidth","meanSilWidth"), 
                                 dims=1, agg.by=NULL, agg.fn=mean, 
                                 filterExpr=NULL, value_format="", 
                                 reorder_rows=FALSE, reorder_columns=TRUE,
                                 show_heatmap_legend=TRUE, 
                                 show_column_names=FALSE, 
                                 col=rev(RColorBrewer::brewer.pal(n=11,"RdBu")),
                                 font_factor=1, row_split=NULL, shortNames=TRUE,
                                 value_cols=c("white","black"), title=NULL, 
                                 anno_legend=TRUE, ...){
  pd <- NULL
  if(is(res,"SimpleList")) pd <- metadata(res)$PipelineDefinition
  if(is.null(pd)) stop("Could not find the PipelineDefinition.")
  if(length(what)>1){
    fcall <- match.call()
    H <- NULL
    for(i in what){
      fcall$what <- i
      H <- H + eval(fcall)
    }
    return(H)
  }
  if(is(res,"SimpleList") && "evaluation" %in% names(res)) 
    res <- res$evaluation
  if(is(res,"list") && "dimreduction" %in% names(res))
    res <- res[["dimreduction"]]
  if(is(res,"list") && "silhouette" %in% names(res))
    res <- res[["silhouette"]]
  res <- res[[dims]]
  what_options <- setdiff( colnames(res), 
                           c("dataset","subpopulation",unlist(arguments(pd))) )
  what <- match.arg(what, choices=what_options)
  res <- .prepRes(res, what=what, agg.by=agg.by, pipDef=pd, returnParams=TRUE,
                  filt=substitute(filterExpr), shortNames=shortNames)
  pp <- res$pp
  res <- as.matrix(res$res)
  ro <- .dosort(res, reorder_rows)
  if(reorder_columns) res <- res[,order(colMeans(res), decreasing=TRUE)]
  
  cellfn <- .getCellFn(res, res, value_format, value_cols, font_factor)
  if(is.null(title)) title <- gsub("\\.","\n",what)
  if(!is.null(row_split)){
    if(row_split %in% colnames(pp)){
      row_split <- pp[,row_split]
    }else{
      warning("`row_split` wasn't found and will be ignored.")
      row_split <- NULL
    }
  }
  col <- .silScale(res, col)
  Heatmap( res, name=what, cluster_rows=FALSE, cluster_columns=FALSE, 
           show_heatmap_legend=show_heatmap_legend, row_order=ro,
           bottom_annotation=.ds_anno(colnames(res),anno_legend,font_factor), 
           show_column_names=show_column_names, cell_fun=cellfn, col=col,
           column_title=title, row_split=row_split,
           row_title_gp = gpar(fontsize = (14*font_factor)),
           column_title_gp = gpar(fontsize = (14*font_factor)),
           row_names_gp = gpar(fontsize = (12*font_factor)),
           column_names_gp = gpar(fontsize = (12*font_factor)), ...)
}

.getClustAggFields <- function(res, pipDef=NULL){
  if(is.null(pipDef)) pipDef <- metadata(res)$PipelineDefinition
  if("evaluation" %in% names(res)) res <- res$evaluation
  if("clustering" %in% names(res)) res <- res$clustering
  fields <- unlist(arguments(pipDef))
  fields <- setdiff(fields, c("k", "dims", "steps","resolution","min.size"))
  fields <- intersect(colnames(res), fields)
  fields[vapply( res[,fields,drop=FALSE], FUN=function(x) length(unique(x))>1,
                 logical(1) )]
}

.mergeFilterOut <- function(ll){
  if(length(ll)==1){
    ll[[1]]$N <- ll[[1]]$N.before
    return(ll[[1]])
  }
  # get back the initial number of cells
  mm <- ll[[1]]
  for(i in seq_len(length(ll)-1)){
    f <- setdiff(colnames(ll[[i]]), c("N.before","N.lost","pc.lost"))
    suf <- paste0(".",names(ll)[i+0:1])
    if(i>1) suf[1] <- ""
    mm <- merge(mm, ll[[i+1]], by=f, suffixes=suf)
  }
  mm$N <- mm[,grep("N\\.before",colnames(mm))[1]]
  mm$N.lost <- rowSums(mm[,grep("N\\.lost",colnames(mm))])
  mm$pc.lost <- 100*rowSums(mm[,grep("N\\.lost",colnames(mm))])/
    mm[,grep("N\\.before",colnames(mm))[1]]
  mm
}

#' scrna_evalPlot_filtering
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param steps Steps to include (default 'doublet' and 'filtering'); other 
#' steps will be averaged.
#' @param clustMetric Clustering accuracy metric to use (default `mean_F1``)
#' @param returnTable Logical; whether to return the data rather than plot.
#'
#' @return A ggplot, or a data.frame if `returnTable=TRUE`
#' @importFrom stats median
#' @export
#' @examples
#' data("exampleResults", package="pipeComp")
#' scrna_evalPlot_filtering(exampleResults)
scrna_evalPlot_filtering <- function(res, steps=c("doublet","filtering"), 
                                     clustMetric="mean_F1", returnTable=FALSE){
  param_fields <- unlist(arguments(metadata(res)$PipelineDefinition)[steps])
  res <- res$evaluation
  co <- .mergeFilterOut(res[steps])
  coI <- co[,c("dataset",param_fields)]
  ci <- split(seq_len(nrow(co)), coI, drop=TRUE)
  agfns <- list(min.lost=min, mean.lost=mean, median.lost=median, max.lost=max)
  x <- as.data.frame(do.call(cbind, lapply( agfns, FUN=function(agf){
    vapply( ci, function(x) agf(co[x,"pc.lost"]), numeric(1) )
  })))
  x$total.lost <- vapply(ci, FUN=function(x) sum(co[x,"N.lost"]), numeric(1))
  x$pc.lost <- vapply( ci, FUN.VALUE=numeric(1), 
                       FUN=function(x) 100*sum(co[x,"N.lost"])/sum(co[x,"N"]) )
  x <- cbind(coI[vapply(ci, FUN=function(x) x[1], integer(1)),], x)
  # get clustering data
  cl <- aggregate( res$clustering[,clustMetric,drop=FALSE], 
                   by=res$clustering[,c("dataset",param_fields)], FUN=mean )
  m <- merge(cl, x, by=c("dataset",param_fields))
  m$method <- apply( m[,param_fields,drop=FALSE], 1, collapse="+", FUN=paste )
  
  if(returnTable) return(m)
  if( length(param_fields)==2 && 
      all(sort(param_fields)==c("doubletmethod","filt")) ){
    return(ggplot(m, aes(max.lost, mean_F1, colour=filt, shape=doubletmethod))+ 
      geom_point(size=3) + facet_wrap(~dataset, scales="free") + 
      xlab("Max proportion of subpopulation excluded") +
      labs(colour="filterset", shape="doublet method"))
  }
  ggplot(m, aes(max.lost, mean_F1, colour=method, shape=method)) + 
    geom_point(size=3) + facet_wrap(~dataset, scales="free") + 
    xlab("Max proportion of subpopulation excluded")
}

#' scrna_describeDatasets
#' 
#' Plots descriptive information about the datasets
#'
#' @param sces A character vector of paths to SCE rds files, or a list of SCEs
#' @param pt.size Point size (for reduced dims)
#' @param ... Passed to geom_point()
#'
#' @return A plot_grid output
#'
#' @export
#' @import SingleCellExperiment scran scater ggplot2 scales
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom cowplot plot_grid
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
scrna_describeDatasets <- function(sces, pt.size=0.3, ...){
  if(is.null(names(sces))) names(sces) <- paste0("dataset",seq_along(sces))
  if(is.character(sces)){
    sces <- lapply(sces, FUN=function(x){
      sce <- readRDS(x)
      sce <- logNormCounts(sce)
      var.stats <- modelGeneVar(sce)
      sce <- denoisePCA(sce, technical=var.stats)
      reducedDim(sce, "tSNE") <- Rtsne(reducedDim(sce))$Y
      reducedDim(sce, "umap") <- umap(reducedDim(sce))
      sce
    })
  }
  tt <- lapply(sces, FUN=function(x) table(x$phenoid))
  cols <- lapply(tt, FUN=function(x){
    y <- getQualitativePalette(length(x))
    names(y) <- names(x)
    y
  })
  for(i in seq_along(sces)) 
    sces[[i]]$cluster <- cols[[i]][as.character(colData(sces[[i]])$phenoid)]
  cd <- lapply(sces, FUN=function(x) as.data.frame(colData(x)))
  names(cols2) <- cols2 <- unique(unlist(cols))
  noy <- theme( axis.line=element_blank(),axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(), axis.title.y=element_blank(),
                  strip.background = element_blank(), legend.position = "none",
                  aspect.ratio = 1, axis.text=element_text(size=10),
                  plot.margin=margin(l=0,r=0))
  cs <- scale_fill_manual(values=cols2)
  d <- data.frame( dataset=rep(names(sces), vapply(tt, length, integer(1))), 
                   cluster=unlist(cols), nb=unlist(tt) )
  p1 <- ggplot(d, aes(x=cluster, y=nb, fill=cluster)) +
    geom_bar(stat = "identity") + scale_y_log10() + 
    facet_wrap(~dataset, scales="free_y", ncol=1, strip.position = "left") + 
    coord_flip() + ylab("Number of cells") + xlab("") +  noy + cs
  for(x in names(sces))
    colData(sces[[x]])$cluster <- 
      unlist(cols)[paste(x, colData(sces[[x]])$phenoid,sep=".")]
  d <- suppressWarnings(dplyr::bind_rows(lapply(cd, FUN=function(x) 
    x[,c("total_counts","total_features","cluster")]))
  )
  d$dataset <- rep(names(cd), vapply(cd, nrow, integer(1)))
  pf <- function(d, x) ggplot(d, aes_string(x="cluster", y=x, fill="cluster"))+ 
    geom_violin() + xlab("") + coord_flip() + noy + cs +
    facet_wrap(~dataset, scales="free_y", ncol=1) + 
    theme( strip.text.x = element_blank(), 
           panel.grid.major.x=element_line(color="grey", linetype="dashed")) + 
    scale_y_continuous( trans = 'log10', breaks=10^pretty(log10(d[[x]]),n=3), 
                        labels=trans_format('log10', math_format(10^.x)) )
  rd <- function(dimred,...){
    d2 <- cbind(d, do.call(rbind, lapply(sces, FUN=function(x){
      reducedDim(x,dimred)
    })))
    colnames(d2)[5:6] <- c("x","y")
    ggplot(d2, aes(x,y,col=cluster)) + geom_point(...) + 
      facet_wrap(~dataset, scales="free", ncol=1) + 
      scale_color_manual(values=cols2) + xlab(paste0("\n",dimred)) +
      theme(axis.line=element_blank(), axis.text=element_blank(),
            axis.ticks=element_blank(), axis.title.y=element_blank(),
            strip.background = element_blank(), strip.text.x=element_blank(),
            legend.position = "none", aspect.ratio = 1 )
  }
  p5 <- rd("umap", size=pt.size, ...)
  plot_grid( p1, pf(d,"total_counts"), pf(d,"total_features"), 
             rd("tSNE", size=pt.size, ...), rd("umap", size=pt.size, ...),
             nrow=1 )
}

#' @importFrom circlize colorRamp2
.silScale <- function(x, cols=NULL){
  if(is.null(cols)) cols <- rev(brewer.pal(n=11,"RdBu"))
  if(is.function(cols)) cols <- cols(11)
  if(length(cols)!=11) stop("`cols` should contain 11 colors.")
  bb <- c( -seq(from=sqrt(abs(min(x, na.rm=TRUE))), to=0, length.out=6)^2, 
           seq(from=0, to=sqrt(max(x, na.rm=TRUE)), length.out=6)[-1]^2 )
  colorRamp2(bb, cols)
}

