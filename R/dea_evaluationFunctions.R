#' evaluateDEA
#' 
#' Evaluates a differential expression analysis (DEA).
#'
#' @param dea Expects a data.frame with logFC and FDR, as produced by `edgeR::topTags`, `limma::topTable` or 
#' `DESeq2::results`.
#' @param truth A data.frame containing the columns `expected.beta` (real logFC)
#' and `isDE` (logical indicating whether there is a difference or not; accepts
#' NA values)
#' @param th The significance thresholds for which to compute the metrics.
#'
#' @return A list with two slots: `logFC` (vector of metrics on logFC) and 
#' `significance` table of significance-related statistics.
#' 
#' @export
#' @examples
#' # fake DEA results
#' dea <- data.frame( row.names=paste0("gene",1:10), logFC=rnorm(10) )
#' dea$PValue <- dea$FDR <- c(2:8/100, 0.2, 0.5, 1)
#' truth <- data.frame( row.names=paste0("gene",1:10), expected.beta=rnorm(10),
#'                      isDE=rep(c(TRUE,FALSE,TRUE,FALSE), c(3,1,2,4)) )
#' evaluateDEA(dea, truth)
evaluateDEA <- function(dea, truth=NULL, th=c(0.01,0.05,0.1)){
  dea <- .homogenizeDEA(dea)
  if(is.null(truth)) truth <- metadata(dea)$truth
  dea <- cbind(dea, truth[row.names(dea),])
  dea <- dea[!is.na(dea$expected.beta),]
  # comparison of estimated and expected log2 folchanges
  res <- c(logFC.pearson=cor(dea$logFC, dea$expected.beta, 
                             use = "pairwise"),
           logFC.spearman=cor(dea$logFC, dea$expected.beta, 
                              use = "pairwise", method="spearman"),
           logFC.mad=median(abs(dea$logFC-dea$expected.beta),na.rm=TRUE),
           ntested=sum(!is.na(dea$PValue) & !is.na(dea$FDR)))
  # evaluation of singificance calls
  names(th) <- th
  res2 <- t(vapply(th, FUN.VALUE=vector(mode="numeric", length=6), FUN=function(x){
    called=sum(dea$FDR<x,na.rm=TRUE)
    P <- sum(dea$isDE)
    TP <- sum(dea$FDR<x & dea$isDE, na.rm=TRUE)
    c( TP=TP, FP=called-TP, TPR=TP/P, PPV=TP/called, FDR=1-TP/called, 
       FPR=(P-TP)/sum(!dea$isDE) )
  }))
  res2 <- cbind(threshold=as.numeric(row.names(res2)), as.data.frame(res2))
  row.names(res2) <- NULL
  list(logFC=res, significance=res2)
}



#' aggregateDEAeval
#'
#' Aggregates DEA results (as produced by `evaluateDEA`) across datasets.
#'
#' @param res A list (per-dataset) of lists of DEA evaluation results.
#'
#' @return A list of aggregated results.
#' @export
aggregateDEAeval <- function(res){
  lfc <- dplyr::bind_rows(lapply(res, FUN=function(x){
    x <- do.call(rbind, lapply(x, FUN=function(x) x$logFC))
    parsePipNames(as.data.frame(x))
  }), .id="dataset")
  sig <- dplyr::bind_rows(lapply(res, FUN=function(x){
    x <- dplyr::bind_rows( lapply(x, FUN=function(x) x$significance), 
                           .id="method" )
    cbind(parsePipNames(x$method), x[,-1])
  }), .id="dataset")
  lapply( list( logFC=lfc, significance=sig ), FUN=function(x){
    x$sva.method <- gsub("^sva\\.","",x$sva.method)
    x$dea.method <- gsub("^dea\\.","",x$dea.method)
    x
  })
}

#' dea_evalPlot_curve
#'
#' @param res 
#' @param scales Passed to `facet_grid`
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param xlim Optional vector of x limits
#' @param colourBy Name of column by which to colour
#' @param shapeBy Name of column determining the shape of the points. If 
#' omitted, the shape will indicate whether the nominal FDR is below or equal 
#' the real FDR.
#'
#' @return A ggplot.
#' @export
#' @examples
#' data("exampleDEAresults", package="pipeComp")
#' dea_evalPlot_curve(exampleDEAresults, agg.by=c("sva.method"))
dea_evalPlot_curve <- function(res, scales="free", agg.by=NULL, agg.fn=mean, 
                               xlim=c(NA,NA), colourBy="method", shapeBy=NULL){
  pd <- NULL
  if(is(res,"SimpleList")) pd <- metadata(res)$PipelineDefinition
  if(is(res,"SimpleList") && "evaluation" %in% names(res)) res <- res$evaluation
  if(is(res,"list") && "dea" %in% names(res)) res <- res$dea
  if(is(res,"list") && "significance" %in% names(res)) res <- res$significance
  if(is.null(agg.by)) agg.by <- intersect(colnames(res), unlist(arguments(pd)))
  d <- aggregate(res[,c("FDR","TPR")], by=res[,c("dataset","threshold",agg.by)], na.rm=TRUE, FUN=agg.fn)
  pp <- d[,agg.by,drop=FALSE]
  d$method <- apply(pp, 1, collapse=" > ", FUN=paste)
  if(all(c("filt","minCount") %in% colnames(d))) d$filter <- paste(d$filt, "(",d$minCount,")")
  p <- ggplot(d, aes_string("FDR", "TPR", group="method", colour=colourBy, 
                            shape=shapeBy)) + 
    geom_vline(xintercept=unique(d$threshold), linetype="dashed", 
               colour="darkgrey") + 
    geom_line(size=1) + geom_point(size=4)
  if(is.null(shapeBy)) p <- p + 
    geom_point(data=d[d$FDR>d$threshold,], size=3, colour="white")
  p + facet_wrap(~dataset, scales=scales) + 
    scale_x_sqrt( breaks=c(0.001,0.01,0.05,0.1,0.2,0.4), 
                  limits=xlim ) + 
    theme(axis.text.x=element_text(angle = 90, hjust = 1))
}

#' dea_evalPlot_sig
#'
#' @param res Aggregated results of the DEA pipeline
#' @param what What to plot
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param scale Logical; whether to scale columns (default FALSE)
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param reorder_rows Logical; whether to sort rows (default TRUE). The row 
#' names themselves can also be passed to specify an order, or a 
#' `ComplexHeatmap`.
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsize of column titles.
#' @param row_split Optional column by which to split rows.
#' @param value_cols A vector of length 2 indicating the colors of the values
#' (above and below the mean), if printed
#' @param title Plot title
#' @param anno_legend Logical; whether to plot the legend for the datasets
#' @param ... Passed to `Heatmap`
#' 
#' @return A Heatmap
#' @export
#' @examples
#' data("exampleDEAresults", package="pipeComp")
#' dea_evalPlot_sig( exampleDEAresults, agg.by=c("sva.method","dea.method"), 
#'                   row_split = "sva.method" )
dea_evalPlot_sig <- function( res, what=c("TPR","FDR"), threshold=0.05,
                              agg.by=NULL, agg.fn=mean, scale=TRUE, 
                              value_format="%.2f", reorder_rows=TRUE,
                              show_heatmap_legend=FALSE,
                              col=viridisLite::inferno(100), 
                              col_title_fontsize=12, show_column_names=FALSE,
                              value_cols=c("black","white"), title=NULL,
                              anno_legend=TRUE, row_split=NULL, ...){
  pd <- NULL
  if(is(res,"SimpleList")) pd <- metadata(res)$PipelineDefinition
  if(length(what)>1){
    ro <- H <- dea_evalPlot_sig(res, what[1], agg.by=agg.by, scale=scale, 
                            reorder_rows=reorder_rows, anno_legend=anno_legend,
                            row_split=row_split, 
                            show_column_names=show_column_names, ...)
    for(i in what[-1]) H <- H +
        dea_evalPlot_sig(res, i, agg.by=agg.by, scale=scale, reorder_rows=ro, 
                         anno_legend=anno_legend, row_split=row_split,
                         show_column_names=show_column_names, ...)
    return(H)
  }
  if(is(res,"SimpleList") && "evaluation" %in% names(res)) res <- res$evaluation
  if(is(res,"list") && "dea" %in% names(res)) res <- res$dea
  if(is(res,"list") && "significance" %in% names(res)) res <- res$significance
  res <- res[res$threshold==threshold,]
  what <- match.arg(what)
  res <- .prepRes(res, what=what, agg.by=agg.by, pipDef=pd, 
                          shortNames=TRUE, returnParams=TRUE)
  pp <- res$pp
  res2 <- res <- res$res
  if(scale) res2 <- .safescale(res)
  if(is(reorder_rows, "Heatmap")){
    ro <- row.names(reorder_rows@matrix)
  }else{
    if(length(reorder_rows)>1){
      ro <- reorder_rows
    }else{
      if(reorder_rows){
        ro <- order(rowMeans(res2, na.rm=TRUE), decreasing=TRUE)
      }else{
        ro <- seq_len(nrow(res2))
      }
    }
  }
  if(!is.null(row_split)){
    if(row_split %in% colnames(pp)){
      row_split <- pp[,row_split]
    }else{
      warning("`row_split` wasn't found and will be ignored.")
      row_split <- NULL
    }
  }
  if(!is.numeric(ro)) ro <- sapply(ro, FUN=function(x) which(row.names(res)==x))
  res2 <- as.matrix(res2)
  cellfn <- .getCellFn(res, res2, value_format, value_cols)
  if(is.null(title)) title <- what
  Heatmap( res2, name=what, cluster_rows=FALSE, cluster_columns=FALSE, 
           show_heatmap_legend=show_heatmap_legend, row_order=ro,
           bottom_annotation=.ds_anno(colnames(res),anno_legend), 
           show_column_names=show_column_names, cell_fun=cellfn, col=col,
           column_title=title, row_split=row_split,
           column_title_gp=gpar(fontisze=col_title_fontsize), ...)
}

#' dea_evalPlot_logFCs
#'
#' @param res Aggregated results of the DEA pipeline
#' @param what What to plot
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param scale Logical; whether to scale columns (default FALSE)
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param reorder_rows Logical; whether to sort rows (default TRUE). The row 
#' names themselves can also be passed to specify an order, or a 
#' `ComplexHeatmap`.
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsize of column titles.
#' @param row_split Optional column by which to split rows.
#' @param value_cols A vector of length 2 indicating the colors of the values
#' (above and below the mean), if printed
#' @param title Plot title
#' @param anno_legend Logical; whether to plot the legend for the datasets
#' @param ... Passed to `Heatmap`
#' 
#' @return A Heatmap
#' @export
#' @examples
#' data("exampleDEAresults", package="pipeComp")
#' dea_evalPlot_logFCs( exampleDEAresults, agg.by=c("sva.method","dea.method"), 
#'                   row_split = "sva.method" )
dea_evalPlot_logFCs <- function(res,
                                what=c("logFC.pearson","logFC.spearman","logFC.mad"), 
                                agg.by=NULL, agg.fn=mean, scale=TRUE, 
                                value_format="%.2f", reorder_rows=TRUE,
                                show_heatmap_legend=FALSE,
                                col=viridisLite::inferno(100), 
                                col_title_fontsize=12, row_split=NULL, 
                                value_cols=c("black","white"), title=NULL,
                                anno_legend=TRUE, ...){
  pd <- NULL
  if(is(res,"SimpleList")) pd <- metadata(res)$PipelineDefinition
  if(length(what)>1){
    ro <- H <- dea_evalPlot_logFCs(res, what[1], agg.by=agg.by, scale=scale, 
                             reorder_rows=reorder_rows, row_split=row_split, 
                             anno_legend=anno_legend, ...)
    for(i in what[-1]) H <- H +
        dea_evalPlot_logFCs(res, i, agg.by=agg.by, scale=scale, reorder_rows=ro, 
                            anno_legend=anno_legend, row_split=row_split, ...)
    return(H)
  }
  if(is(res,"SimpleList") && "evaluation" %in% names(res)) res <- res$evaluation
  if(is(res,"list") && "dea" %in% names(res)) res <- res$dea
  if(is(res,"list") && "logFC" %in% names(res)) res <- res$logFC
  what <- match.arg(what)
  res <- .prepRes(res, what=what, agg.by=agg.by, pipDef=pd, 
                          shortNames=TRUE, returnParams = TRUE)
  pp <- res$pp
  res2 <- res <- res$res
  if(scale) res2 <- .safescale(res)
  if(is(reorder_rows, "Heatmap")){
    ro <- row.names(reorder_rows@matrix)
  }else{
    if(length(reorder_rows)>1){
      ro <- reorder_rows
    }else{
      if(reorder_rows){
        ro <- order(rowMeans(res2, na.rm=TRUE), decreasing=TRUE)
      }else{
        ro <- seq_len(nrow(res2))
      }
    }
  }
  if(!is.numeric(ro)) ro <- sapply(ro, FUN=function(x) which(row.names(res)==x))
  res2 <- as.matrix(res2)
  cellfn <- .getCellFn(res, res2, value_format, value_cols)
  if(is.null(title)) title <- gsub("\\.","\n",what)
  if(!is.null(row_split)){
    if(row_split %in% colnames(pp)){
      row_split <- pp[,row_split]
    }else{
      warning("`row_split` wasn't found and will be ignored.")
      row_split <- NULL
    }
  }
  Heatmap( res2, name=what, cluster_rows=FALSE, cluster_columns=FALSE, 
           show_heatmap_legend=show_heatmap_legend, row_order=ro,
           bottom_annotation=.ds_anno(colnames(res),anno_legend), 
           show_column_names = FALSE, cell_fun=cellfn, col=col,
           column_title=title, row_split=row_split,
           column_title_gp=gpar(fontisze=col_title_fontsize), ...)
}