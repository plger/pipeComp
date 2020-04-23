#' evalHeatmap
#' 
#' General heatmap representation of aggregated evaluation results.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param step Name of the step for which to plot the evaluation results
#' @param what What metric to plot
#' @param what2 If the step has more than one benchmark data.frame, which one
#' to use.
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param scale Logical; whether to scale columns (default FALSE)
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param reorder_rows Logical; whether to sort rows (default FALSE). The row 
#' names themselves can also be passed to specify an order, or a 
#' `ComplexHeatmap`.
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
#' data("exampleDEAresults", package="pipeComp")
#' evalHeatmap( exampleDEAresults, what=c("logFC.pearson", "TPR","FDR"), 
#'              agg.by=c("sva.method","dea.method"), row_split = "sva.method" )
evalHeatmap <- function( res, step=NULL, what, what2=NULL, agg.by=NULL, 
                         agg.fn=mean, scale=TRUE, value_format="%.2f", 
                         reorder_rows=FALSE, show_heatmap_legend=FALSE, 
                         show_column_names=FALSE, col=viridisLite::inferno(100),
                         font_factor=1, row_split=NULL, shortNames=TRUE,
                         value_cols=c("black","white"), title=NULL, 
                         anno_legend=TRUE, ...){
  pd <- NULL
  if(is(res,"SimpleList")) pd <- metadata(res)$PipelineDefinition
  if(is.null(pd)) stop("Could not find the PipelineDefinition.")
  if(is.null(step)) step <- rev(names(res$evaluation))[1]
  step <- match.arg(step, names(res$evaluation))
  if(length(what)>1){
    ro <- H <- evalHeatmap(res, step=step, what[1], agg.by=agg.by, scale=scale, 
                           reorder_rows=reorder_rows, row_split=row_split, 
                           show_column_names=show_column_names,
                           font_factor=font_factor, shortNames=shortNames,
                           show_heatmap_legend=show_heatmap_legend,
                           anno_legend=anno_legend, ... )
    for(i in what[-1]) H <- H +
        evalHeatmap( res, step=step, what=i, agg.by=agg.by, scale=scale, 
                     reorder_rows=ro, anno_legend=anno_legend, 
                     show_column_names=show_column_names,
                     show_heatmap_legend=show_heatmap_legend,
                     font_factor=font_factor, shortNames=shortNames,
                     row_split=row_split, ... )
    return(H)
  }
  if(is(res,"SimpleList") && "evaluation" %in% names(res)) res <- res$evaluation
  if(is(res,"list") && step %in% names(res)) res <- res[[step]]
  if(is(res,"list")){
    if(is.null(what2)){
      # try to guess the table
      test <- which(vapply( res, FUN=function(x) what %in% colnames(x), 
                            logical(1) ))
      if(length(test)==0) stop("`what2` is undefined and could not be guessed.")
      if(length(test)>1) message("`what2` is undefined, and `what` is present ",
        "in more than one results data.frame. '", test[1], "' will be used.")
      what2 <- names(res)[test[1]]
    }
    what2 <- match.arg(what2, names(res))
    res <- res[[what2]]
  }
  what_options <- setdiff( colnames(res), 
                           c("dataset",unlist(arguments(pd))) )
  what <- match.arg(what, choices=what_options)
  res <- .prepRes(res, what=what, agg.by=agg.by, pipDef=pd, 
                  shortNames=shortNames, returnParams = TRUE)
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
  if(!is.numeric(ro))
    ro <- vapply(ro, FUN=function(x) which(row.names(res)==x), integer(1))
  res2 <- as.matrix(res2)
  cellfn <- .getCellFn(res, res2, value_format, value_cols, font_factor)
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
           bottom_annotation=.ds_anno(colnames(res),anno_legend,font_factor), 
           show_column_names=show_column_names, cell_fun=cellfn, col=col,
           column_title=title, row_split=row_split,
           row_title_gp = gpar(fontsize = (14*font_factor)),
           column_title_gp = gpar(fontsize = (14*font_factor)),
           row_names_gp = gpar(fontsize = (12*font_factor)),
           column_names_gp = gpar(fontsize = (12*font_factor)), ...)
}

# NA-robust scale
.safescale <- function(x){
  if(!is.null(dim(x))){
    y <- apply(x,2,.safescale)
    row.names(y) <- row.names(x)
    return(y)
  }
  if(all(is.na(x))) return(x)
  if(sd(x,na.rm=TRUE)>0) return(base::scale(x))
  if(sum(!is.na(x))==0) return(base::scale(as.numeric(!is.na(x))))
  rep(0,length(x))
}

#' @importFrom ComplexHeatmap HeatmapAnnotation
.ds_anno <- function(x, legend=TRUE, font_factor=1){
  y <- vapply( strsplit(gsub("stepElapsed\\.","",x)," "),
               FUN=function(x) x[[1]], FUN.VALUE=character(1))
  cols <- getQualitativePalette(length(unique(y)))
  names(cols) <- sort(unique(y))
  HeatmapAnnotation(dataset=y, col=list(dataset=cols), 
                    annotation_legend_param=list(
                      labels_gp=gpar(fontsize=10*font_factor),
                      title_gp=gpar(fontsize=10*font_factor, fontface = "bold")
                    ),
                    show_annotation_name = FALSE, show_legend=legend)
}

.getReducedNames <- function(res){
  if(is.character(res)) res <- parsePipNames(res)
  res <- res[,grep("^stepElapsed\\.",colnames(res),invert=TRUE)]
  pp <- pp[,apply(pp,2,FUN=function(x) length(unique(x))>1),drop=FALSE]
  if(ncol(pp)>1){
    y <- apply(pp,1,FUN=function(x){
      x <- paste0(colnames(pp),"=",x)
      paste(x, collapse("; "))
    })
  }else{
    y <- apply(pp,1,collapse=" ",FUN=paste)
  }
  y
}

.prepRes <- function(res, what=NULL, agg.by=NULL, agg.fn=mean, pipDef=NULL, 
                     long=FALSE, shortNames=FALSE, returnParams=FALSE){
  if(is.null(what) && !long) stop("`what` should be defined.")
  if(!is.null(agg.by)){
    agg.by2 <- c(agg.by, intersect(colnames(res),c("dataset","subpopulation")))
    if(is.null(what)){
      what <- colnames(res)[which(
        vapply(res, function(x) is.integer(x) | is.numeric(x), logical(1))
      )]
      what <- setdiff(what, agg.by)
    }
    res <- aggregate(res[,what,drop=FALSE], by=res[,agg.by2,drop=FALSE], 
                     FUN=agg.fn)
    pp <- res[,agg.by,drop=FALSE]
  }else{
    if(is.null(pipDef)) stop("Either `agg.by` or `pipDef` should be given.")
    pp <- res[,intersect(colnames(res), unlist(arguments(pipDef))),drop=FALSE]
  }
  pp$method <- res$method <- .getReducedNames(pp, short=shortNames)
  if(long) return(res)
  
  if("subpopulation" %in% colnames(res))
    res$dataset <- paste(res$dataset, res$subpopulation)
  r2 <- reshape2::dcast( res, method~dataset, value.var=what,  
                         drop=FALSE, fun.aggregate=agg.fn)
  pp <- pp[!duplicated(pp$method),,drop=FALSE]
  row.names(pp) <- pp$method
  row.names(r2) <- r2[,1]
  res <- as.matrix(r2[,-1,drop=FALSE])
  if(returnParams) return(list(res=res, pp=pp[row.names(res),,drop=FALSE]))
  res
}

.getReducedNames <- function(res, short=FALSE){
  if(is.character(res)) res <- parsePipNames(res)
  pp <- vapply(res, FUN=function(x) length(unique(x))>1, logical(1))
  pp <- res[,pp,drop=FALSE]
  if(!short && ncol(pp)>1){
    y <- apply(pp,1,FUN=function(x){
      x <- paste0(colnames(pp),"=",x)
      paste(x, collapse="; ")
    })
  }else{
    y <- apply(pp,1,collapse=" > ",FUN=paste)
  }
  y
}
.getCellFn  <- function( res, res2, value_format, cols=c("black","white"), 
                         font_factor=1 ){
  resmid <- range(res2, na.rm=TRUE)
  resmid <- resmid[1]+(resmid[2]-resmid[1])/2
  function(j, i, x, y, width, height, fill){
    if(value_format=="" || is.null(value_format) || is.na(value_format))
      return(NULL)
    lab <- sprintf(value_format, res[i,j])
    if(!any(abs(res)>1, na.rm=TRUE)){
      lab <- gsub("^0\\.",".",lab)
      lab <- gsub("^-0\\.","-.",lab)
    } 
    lab <- gsub("^1.00$","1",lab)
    lab <- gsub("^.00$","0",lab)
    cols <- ifelse(res2[i,j]>resmid,cols[1],cols[2])
    grid.text(lab, x, y, gp = gpar(fontsize = 10*font_factor, col=cols))
  }
}
