#' evalHeatmap
#' 
#' General heatmap representation of aggregated evaluation results. By default,
#' the actual metric values are printed in the cells, and while the coloring is
#' determined by \code{\link{colCenterScale}} (number of matrix-median absolute
#'  deviations from the column means). Unless the total number of analyses is
#'  small, it is strongly recommended to use the `agg.by` argument to limit the
#'  size and improve the readability of the heatmap.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param step Name of the step for which to plot the evaluation results. If 
#' unspecified, will use the latest step that has evaluation results.
#' @param what What metric to plot.
#' @param what2 If the step has more than one benchmark data.frame, which one
#' to use. The function will attempt to guess that automatically based on 
#' `what`, and will notify in case of ambiguity.
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param filterExpr An optional filtering expression based on the columns of 
#' the target dataframe, (e.g. `filterExpr=param1=="value1"`).
#' @param scale Controls the scaling of the columns for the color mapping. Can 
#' either be a logical (TRUE will use NA-safe column z-scores, FALSE will not 
#' scale) or a function performing the scaling. The default uses the
#'  `colCenterScale` function (per-column centering, but per-matrix variance
#'  scaling).
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param reorder_rows Logical; whether to sort rows (default FALSE). The row 
#' names themselves can also be passed to specify an order, or a 
#' `ComplexHeatmap`.
#' @param row_split Optional column (included in `agg.by`) by which to split
#' the rows. Alternatively, an expression using the columns (retained after
#' aggregation) can be passed.
#' @param show_heatmap_legend Passed to `Heatmap` (default FALSE)
#' @param show_column_names Passed to `Heatmap` (default FALSE)
#' @param col Colors for the heatmap. By default, will apply linear mapping (if
#' the data is not scaled) or signed sqrt mapping (if scaled) on the 
#' `viridisLite::inferno` palette. To disable the signed sqrt-transformation,
#' simply pass `col=viridisLite::inferno(100)` or your own palette.
#' @param font_factor A scaling factor applied to fontsizes (default 1)
#' @param value_cols A vector of length 2 indicating the colors of the values
#' (above and below the mean), if printed
#' @param title Plot title
#' @param shortNames Logical; whether to use short row names (with only
#' the parameter values instead of the parameter name and value pairs), default
#' TRUE.
#' @param name Heatmap name (e.g. used for the legend)
#' @param anno_legend Logical; whether to plot the legend for the datasets
#' @param ... Passed to `Heatmap`
#' 
#' @return A Heatmap
#' @export
#' @importFrom viridisLite inferno
#' @importFrom grid gpar
#' @import ComplexHeatmap
#' @examples
#' data("exampleResults", package="pipeComp")
#' evalHeatmap( exampleResults, step="clustering", what=c("ARI","MI","min_pr"), 
#'              agg.by=c("filt", "norm"), row_split = "norm" ) +
#' evalHeatmap( exampleResults, step="clustering", what="ARI", 
#'              agg.by=c("filt", "norm"), filterExpr=n_clus==true.nbClusts, 
#'              name="ARI at true k", title="ARI at
#' true K" )
evalHeatmap <- function( res, step=NULL, what, what2=NULL, agg.by=NULL, 
                         agg.fn=mean, filterExpr=NULL, scale="colCenterScale", 
                         value_format="%.2f", reorder_rows=FALSE, 
                         show_heatmap_legend=FALSE, show_column_names=FALSE, 
                         col=NULL, font_factor=0.9, 
                         row_split=NULL, shortNames=TRUE,
                         value_cols=c("black","white"), title=NULL, 
                         name=NULL, anno_legend=TRUE, ...){
  pd <- NULL
  if(is(res,"SimpleList")) pd <- metadata(res)$PipelineDefinition
  #if(is.null(pd)) stop("Could not find the PipelineDefinition.")
  if(is.null(step)){
    step <- rev(names(res$evaluation))[1]
    if(length(res$evaluation)>1) message("Using step '",step,"'.")
  }
  step <- match.arg(step, names(res$evaluation))
  if(length(what)>1){
    title <- lapply(seq_along(what), FUN=function(i){
      tryCatch( title[[i]], error=function(e) what[[i]] )
    })
    fcall <- match.call()
    H <- NULL
    for(i in seq_along(what)){
      fcall$what <- what[[i]]
      fcall$title <- title[[i]]
      H <- H + eval(fcall)
    }
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
        "in more than one results data.frame (",
        paste(paste0("`",names(res)[test],"`"), collapse=", "), "). ",
        names(res)[test[1]], "' will be used.")
      what2 <- names(res)[test[1]]
    }
    what2 <- match.arg(what2, names(res))
    res <- res[[what2]]
  }
  if(!is.null(pd)){
    what_options <- setdiff( colnames(res), 
                           c("dataset",unlist(arguments(pd))) )
    what <- match.arg(what, choices=what_options)
  }
  res <- .prepRes(res, what=what, agg.by=agg.by, pipDef=pd, 
                  filt=substitute(filterExpr), shortNames=shortNames, 
                  returnParams=TRUE)
  pp <- res$pp
  res <- res$res
  res2 <- .doscale(res, param=scale)
  ro <- .dosort(res2, reorder_rows)
  res2 <- as.matrix(res2)
  cellfn <- .getCellFn(res, res2, value_format, value_cols, font_factor)
  if(is.null(title)) title <- gsub("\\.","\n",what)
  suppressWarnings({
    if(!tryCatch(is.null(row_split), error=function(e) FALSE)){
      if(tryCatch(is.character(row_split), error=function(e) FALSE)){
        if(row_split %in% colnames(pp)){
          row_split <- pp[,row_split]
        }else{
          warning("`row_split` wasn't found and will be ignored.")
          row_split <- NULL
        }
      }else{
        row_split <- eval(substitute(row_split), envir=pp)
      }
    }
  })
  if(is.null(name)) name <- what
  if(is.null(col))
    col <- .defaultColorMapping(res2, center=!(is.logical(scale) && !scale))
  Heatmap( res2, name=name, cluster_rows=FALSE, cluster_columns=FALSE, 
           show_heatmap_legend=show_heatmap_legend, row_order=ro,
           bottom_annotation=.ds_anno(colnames(res),anno_legend,font_factor), 
           show_column_names=show_column_names, cell_fun=cellfn, col=col,
           column_title=title, row_split=row_split,
           row_title_gp = gpar(fontsize = (14*font_factor)),
           column_title_gp = gpar(fontsize = (14*font_factor)),
           row_names_gp = gpar(fontsize = (12*font_factor)),
           column_names_gp = gpar(fontsize = (12*font_factor)), ...)
}

.doscale <- function(res, param){
  if(is.null(param)) param <- colCenterScale
  if(is.logical(param) && !param) return(res)
  if(is.logical(param) && param) return(.safescale(res))
  if(is.function(param)) return(param(res))
  if(is.character(param) && is.function(get(param))) return(get(param)(res))
  warning("Unknown scaling parameter, will use column z-scores.")
  .safescale(res)
}

#' @import ComplexHeatmap
.dosort <- function(res, reorder_rows){
  if(is(reorder_rows, "Heatmap")){
    ro <- row_order(reorder_rows)
  }else{
    if(length(reorder_rows)>1){
      ro <- reorder_rows
    }else{
      if(reorder_rows){
        ro <- order(rowMeans(res, na.rm=TRUE), decreasing=TRUE)
      }else{
        ro <- seq_len(nrow(res))
      }
    }
  }
  if(!is.numeric(ro))
    ro <- vapply(ro, FUN=function(x) which(row.names(res)==x), integer(1))
  ro
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

.defaultColorMapping <- function(x, center=TRUE){
  if(!center) return(viridisLite::inferno(101))
  q <- max(abs(range(x, na.rm=TRUE)))
  b <- c( -seq(from=sqrt(q), to=0, length.out=51)^2,
          seq(from=0, to=sqrt(q), length.out=51)[-1]^2 )
  colorRamp2(b, viridisLite::inferno(101))
}

#' colCenterScale
#'
#' Matrix scaling by centering columns separately and then performing variance 
#' scaling on the whole matrix, in a NA-robust fashion. With the default 
#' arguments, the output will be the number of (matrix-)median absolute 
#' deviations from the column-median.
#'
#' @param x A numeric matrix.
#' @param centerFn The function for calculating centers. Should accept the 
#' `na.rm` argument. E.g. `centerFn=mean` or `centerFn=median`.
#' @param scaleFn The function for calculating the (matrix-wise) scaling
#' factor. Should accept the `na.rm` argument. Default `median(abs(x))`.
#'
#' @return A scaled matrix of the same dimensions as `x`.
#' @export
#'
#' @examples
#' # random data with column mean differences
#' d <- cbind(A=rnorm(5, 10, 2), B=rnorm(5, 20, 2), C=rnorm(5,30, 2))
#' colCenterScale(d)
colCenterScale <- function(x, centerFn=median, 
                           scaleFn=function(x,na.rm) median(abs(x),na.rm=na.rm)
                           ){
  if(is.null(dim(x))) stop("`x` should be a numeric matrix or data.frame.")
  centers <- apply(x, MARGIN=2, na.rm=TRUE, FUN=centerFn)
  x2 <- t(t(x)-centers)
  x2 <- x2/scaleFn(x2, na.rm=TRUE)
  if(all(is.na(x2))) x2[is.na(x2)] <- 0
  x2
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

.prepRes <- function(res, what=NULL, agg.by=NULL, agg.fn=mean, pipDef=NULL, 
                     filt=NULL, long=FALSE, shortNames=FALSE, 
                     returnParams=FALSE){
  if(is.null(what) && !long) stop("`what` should be defined.")
  if(!is.null(filt) && filt!=TRUE){
    if(is.null(pipDef)) stop("Filtering requires providing a `pipDef`.")
    res <- res[which(eval(filt, res)),]
  }
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

#' @importFrom grid grid.text
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

.scaledLegend <- function(){
  ComplexHeatmap::Legend(
    col_fun=circlize::colorRamp2(1:100, viridisLite::inferno(100)), 
    at=c(1,50,100), title="MADs", labels=c("worst","median","best"), 
    direction="horizontal" )
}