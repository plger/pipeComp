#' scrna_evalPlot_DR
#' 
#' Plotting aggregated evaluation results at the level of dimensionality 
#' reduction for the scRNA pipelines.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param what What to plot (default plots main metrics)
#' @param covar Covariate of interest (used only if `what` is  'covar' or 
#' 'covarRes'
#' @param reorder_rows Logical; whether to sort rows (default TRUE)
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param scale Logical; whether to scale columns (default FALSE)
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsite of the column titles
#' @param ... Passed to `Heatmap`
#'
#' @return One or several `Heatmap` object.
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom viridisLite inferno
scrna_evalPlot_DR <- function(res, what=c("auto","silhouette", "covar", "covarRes", "varExpl", "elapsed"),
                              covar=c("log10_total_counts","log10_total_features",
                                      "total_features"), reorder_rows=TRUE,
                              agg.by=NULL, agg.fn=mean, scale=FALSE, 
                              show_heatmap_legend=FALSE, value_format="%.2f", 
                              col=NULL, col_title_fontsize=11, ...){
  what <- match.arg(what)
  covar <- covar[1]
  if("dimreduction" %in% names(res)) res <- res$dimreduction
  if(what=="auto"){
    # v <- lapply( c("PCtop5.R2", paste0("PC1_covar.", c("log10_total_counts", "total_features"))), FUN=function(x){
    #   x <- res[[x]][,grep("^stepElapsed", colnames(res[[x]]), invert=TRUE)]
    #   as.numeric(as.matrix(x))
    # })
    # corcol <- circlize::colorRamp2(range(v, na.rm=TRUE), c("darkblue","orange"))
    return( 
      scrna_evalPlot_DR(res, "silhouette", scale=scale, reorder_rows=FALSE, show_heatmap_legend=FALSE, col_title_fontsize=col_title_fontsize, ...) + 
      scrna_evalPlot_DR(res, "varExpl", scale=scale, reorder_rows=FALSE, show_heatmap_legend=FALSE, col_title_fontsize=col_title_fontsize, ...) + 
      scrna_evalPlot_DR(res, "covarRes", scale=scale, reorder_rows=FALSE, show_heatmap_legend=FALSE, col_title_fontsize=col_title_fontsize, ...) +
      scrna_evalPlot_DR(res, "covarRes", covar="total_features", scale=scale, reorder_rows=FALSE, show_heatmap_legend=FALSE, col_title_fontsize=col_title_fontsize, ...)
    )
  }
  el <- grep("^stepElapsed\\.", colnames(res[[1]]), value=TRUE)
  res <- switch( what,
                 silhouette=res$clust.avg.silwidth.top_10_dims,
                 covar=res[[paste0("PC1_covar.",covar)]],
                 covarRes=res[[paste0("PC1_covarR.",covar)]],
                 varExpl=res$PCtop5.R2,
                 elapsed=res$PCtop5.R2[,el],
                 stop("Unknown plot type requested")
                 )
  res2 <- res <- .prepRes(res, agg.by, agg.fn, elapsed=what=="elapsed")
  if(scale) res2 <- base::scale(res)
  if(reorder_rows){
    ro <- order(rowMeans(res2), decreasing=TRUE)
  }else{
    ro <- seq_len(nrow(res2))
  }
  co <- order(colMeans(res), decreasing=TRUE)
  res <- res[ro,co]
  res2 <- res2[ro,co]
  name <- switch( what,
                 silhouette="silhouette\nwidth",
                 covar=covar,
                 covarRes=covar,
                 varExpl="Variance\nexplained",
                 elapsed="Running\ntime (s)",
                 NULL
  )
  title <- switch( what,
                 silhouette="average silhouette width\nper subpopluation",
                 covar=paste0("correlation with\n",covar),
                 covarRes=paste0("residual corr with\n",covar),
                 varExpl="var explained by\nsubpopulations",
                 elapsed="Running time (s)",
                 NULL
  )
  if(what=="silhouette" && (is.null(col) || length(col)==11) && !scale) return(
    Heatmap( res2, name=name, cluster_rows=FALSE, col=.silScale(res2, col), 
             cluster_columns=FALSE, bottom_annotation=.ds_anno(colnames(res)), 
             show_column_names = FALSE, 
             cell_fun=.getCellFn(res,res2,value_format,c("white","black")), 
             show_heatmap_legend=show_heatmap_legend, column_title=title, 
             column_title_gp=gpar(fontsize=col_title_fontsize), ...))

  if(is.null(col)) col <- viridisLite::inferno(100)
  Heatmap( res2, name=name, cluster_rows=FALSE,
           cluster_columns=FALSE, bottom_annotation=.ds_anno(colnames(res)), 
           show_column_names = FALSE, cell_fun=.getCellFn(res,res2,value_format),
           col=col, show_heatmap_legend=show_heatmap_legend, column_title=title, 
           column_title_gp=gpar(fontsize=col_title_fontsize), ...)
}

#' scrna_evalPlot_clust
#' 
#' Plotting aggregated evaluation results at the level of clustering for the 
#' scRNA pipelines.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param what What to plot (default plots main metrics)
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param scale Logical; whether to scale columns (default FALSE)
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param reorder_rows Logical; whether to sort rows (default TRUE)
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsize of column titles.
#' @param ... Passed to `Heatmap`
#'
#' @return One or several `Heatmap` object.
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom viridisLite inferno
scrna_evalPlot_clust <- function(res, what="auto", agg.by=NULL, agg.fn=mean, 
                                 scale=FALSE, value_format="%.2f", reorder_rows=TRUE, 
                                 show_heatmap_legend=FALSE,
                                 col=viridisLite::inferno(100), 
                                 col_title_fontsize=12, ...){
  if("clustering" %in% names(res)) res <- res$clustering
  what_options <- sapply(strsplit(colnames(res)[grep(" ",colnames(res))]," "),FUN=function(x) x[[2]])
  what <- match.arg(what, c(what_options, "elapsed", "auto"))
  if(what=="auto"){
    return( 
      scrna_evalPlot_clust(res, "ARI", scale=scale, reorder_rows=FALSE, ...) + 
      scrna_evalPlot_clust(res, "NMI", scale=scale, reorder_rows=FALSE, ...) + 
      scrna_evalPlot_clust(res, "mean_re", scale=scale, reorder_rows=FALSE, ...) + 
      scrna_evalPlot_clust(res, "min_re", scale=scale, reorder_rows=FALSE, ...) +
      scrna_evalPlot_clust(res, "mean_pr", scale=scale, reorder_rows=FALSE, ...) + 
      scrna_evalPlot_clust(res, "min_pr", scale=scale, reorder_rows=FALSE, ...)
    )
  }
  res <- res[,grep(ifelse(what=="elapsed","stepElapsed",paste0(" ",what)), colnames(res))]
  res2 <- res <- .prepRes(res, agg.by, agg.fn, elapsed=what=="elapsed")
  if(scale) res2 <- base::scale(res)
  if(reorder_rows){
    ro <- order(rowMeans(res2), decreasing=TRUE)
  }else{
    ro <- seq_len(nrow(res2))
  }
  co <- order(colMeans(res), decreasing=TRUE)
  res <- res[ro,co]
  res2 <- res2[ro,co]
  cellfn <- .getCellFn(res,res2,value_format)
  title <- switch( what,
                   elapsed="Running time (s)",
                   gsub("_re$","\nrecall",gsub("_pr$","\nprecision",what))
  )
  row.names(res2) <- gsub("resolution=", "res=", gsub("norm=norm.","",row.names(res2),fixed=TRUE))
  Heatmap( res2, name=what, cluster_rows=FALSE, show_heatmap_legend=show_heatmap_legend, 
           cluster_columns=FALSE, bottom_annotation=.ds_anno(colnames(res)), 
           show_column_names = FALSE, cell_fun=cellfn, col=col,
           column_title=title, column_title_gp=gpar(fontisze=col_title_fontsize), ...)
}

#' scrna_evalPlot_clustAtTrueK
#' 
#' Plotting aggregated evaluation results at the level of clustering for the 
#' scRNA pipelines, restricting to those analysis that have the right number of
#' clusters.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param what Metric to plot
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param scale Logical; whether to scale columns (default FALSE)
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param reorder_rows Logical; whether to sort rows (default TRUE)
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsize of the column titles
#' @param ... Passed to `Heatmap`
#'
#' @return A `Heatmap` object.
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom viridisLite inferno
scrna_evalPlot_clustAtTrueK <- function(res, what="ARI", agg.by=NULL, 
                                        agg.fn=mean, scale=FALSE, 
                                        value_format="%.2f", reorder_rows=TRUE, 
                                        show_heatmap_legend=FALSE,
                                        col=viridisLite::inferno(100), 
                                        col_title_fontsize=11, ...){
  if("clustering" %in% names(res)) res <- res$clustering
  pp <- parsePipNames(row.names(res))
  if(is.null(agg.by)){
    pp <- pp[,setdiff(colnames(pp), c("dims","k","steps","resolution","min.size"))]
    pp <- pp[,apply(pp,2,FUN=function(x) length(unique(x)))>1,drop=FALSE]
    agg.by <- colnames(pp)
  }
  if(length(agg.by)>0){
    ds <- gsub(paste0(" ",what), "", grep(paste0(" ",what), colnames(res), value=TRUE, fixed=TRUE), fixed=TRUE)
    pp <- pp[,agg.by,drop=FALSE]
    ps <- split(seq_len(nrow(pp)), .getReducedNames(pp))
    res <- t(sapply(ps, FUN=function(x){
      x <- vapply(ds, FUN.VALUE=double(1), FUN=function(y){
        w <- which(res[,paste(y,"n_clus")]==res[,paste(y,"true.nbClusts")])
        suppressWarnings(agg.fn( res[w,paste(y,what)] ))
      })
      names(x) <- ds
      x
    }))
  }else{
    res <- res[,grep(paste0(" ",what), colnames(res))]
  }
  res2 <- res
  if(scale) res2 <- base::scale(res)
  if(reorder_rows){
    ro <- order(rowMeans(res2, na.rm=TRUE), decreasing=TRUE)
  }else{
    ro <- seq_len(nrow(res2))
  }
  co <- order(colMeans(res), decreasing=TRUE)
  res <- res[ro,co]
  res2 <- res2[ro,co]
  cellfn <- .getCellFn(res,res2,value_format)
  title <- paste(what, "at true\n# of clusters")
  row.names(res2) <- gsub("resolution=", "res=", gsub("norm=norm.","",row.names(res2),fixed=TRUE))
  Heatmap( res2, name=what, cluster_rows=FALSE, show_heatmap_legend=show_heatmap_legend, 
           cluster_columns=FALSE, bottom_annotation=.ds_anno(colnames(res)), 
           show_column_names = FALSE, cell_fun=cellfn, col=col,
           column_title=title, column_title_gp=gpar(fontisze=col_title_fontsize), ...)
}


.ds_anno <- function(x){
  y <- sapply(strsplit(gsub("stepElapsed\\.","",x)," "),FUN=function(x) x[[1]])
  cols <- GTscripts::getQualitativePalette(length(unique(y)))
  names(cols) <- sort(unique(y))
  ComplexHeatmap::HeatmapAnnotation(dataset=y, col=list(dataset=cols), show_annotation_name = FALSE)
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

.prepRes <- function(res, agg.by, agg.fn, elapsed=FALSE){
  if(elapsed){
    res <- as.matrix(res[,grep("stepElapsed",colnames(res))])
  }else{
    res <- as.matrix(res[,grep("stepElapsed",colnames(res),invert=TRUE)])
  }
  pp <- parsePipNames(row.names(res))
  if(!is.null(agg.by)){
    ag <- aggregate(res, by=pp[,agg.by,drop=FALSE], FUN=agg.fn)
    pp <- ag[,agg.by,drop=FALSE]
    res <- ag[,setdiff(colnames(ag),agg.by)]
  }
  row.names(res) <- .getReducedNames(pp)
  res
}

.ds_anno <- function(x){
  y <- sapply(strsplit(gsub("stepElapsed\\.","",x)," "),FUN=function(x) x[[1]])
  cols <- GTscripts::getQualitativePalette(length(unique(y)))
  names(cols) <- sort(unique(y))
  ComplexHeatmap::HeatmapAnnotation(dataset=y, col=list(dataset=cols), show_annotation_name = FALSE)
}

.getReducedNames <- function(res){
  if(is.character(res)) res <- parsePipNames(res)
  pp <- res[,grep("^stepElapsed\\.",colnames(res),invert=TRUE),drop=FALSE]
  pp <- pp[,apply(pp,2,FUN=function(x) length(unique(x))>1),drop=FALSE]
  if(ncol(pp)>1){
    y <- apply(pp,1,FUN=function(x){
      x <- paste0(colnames(pp),"=",x)
      paste(x, collapse="; ")
    })
  }else{
    y <- apply(pp,1,collapse=" ",FUN=paste)
  }
  y
}
.getCellFn  <- function(res,res2,value_format,cols=c("black","white")){
  resmid <- range(res2, na.rm=TRUE)
  resmid <- resmid[1]+(resmid[2]-resmid[1])/2
  function(j, i, x, y, width, height, fill){
    if(value_format=="" || is.null(value_format) || is.na(value_format)) return(NULL)
    lab <- sprintf(value_format, res[i,j])
    if(!any(abs(res)>1)){
      lab <- gsub("^0\\.",".",lab)
      lab <- gsub("^-0\\.","-.",lab)
    } 
    grid.text(lab, x, y, gp = gpar(fontsize = 10, col=ifelse(res2[i,j]>resmid,cols[1],cols[2])))
  }
}

#' describeDatasets
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
#' @import SummarizedExperiment SingleCellExperiment scran ggplot2 scales
#' @importFrom scater logNormCounts
#' @importFrom cowplot plot_grid
describeDatasets <- function(sces, pt.size=0.3, ...){
  if(is.null(names(sces))) names(sces) <- paste0("dataset",seq_along(sces))
  if(is.character(sces)){
    sces <- lapply(sces, FUN=function(x){
      sce <- readRDS(x)
      sce <- scater::logNormCounts(sce)
      var.stats <- modelGeneVar(sce)
      sce <- denoisePCA(sce, technical=var.stats)
      reducedDim(sce, "tSNE") <- Rtsne::Rtsne(reducedDim(sce))$Y
      reducedDim(sce, "umap") <- uwot::umap(reducedDim(sce))
      sce
    })
  }
  tt <- lapply(sces, FUN=function(x) table(x$phenoid))
  cols <- lapply(tt, FUN=function(x){
    y <- GTscripts::getQualitativePalette(length(x))
    names(y) <- names(x)
    y
  })
  for(i in seq_along(sces)) sces[[i]]$cluster <- cols[[i]][as.character(sces[[i]]$phenoid)]
  cd <- lapply(sces, FUN=function(x) as.data.frame(colData(x)))
  names(cols2) <- cols2 <- unique(unlist(cols))
  noy <- theme( axis.line=element_blank(),axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(), axis.title.y=element_blank(),
                  strip.background = element_blank(), legend.position = "none",
                  aspect.ratio = 1, axis.text=element_text(size=10),
                  plot.margin=margin(l=0,r=0))
  cs <- scale_fill_manual(values=cols2)
  d <- data.frame(dataset=rep(names(sces),sapply(tt,length)), cluster=unlist(cols), nb=unlist(tt))
  p1 <- ggplot(d, aes(x=cluster, y=nb, fill=cluster)) + geom_bar(stat = "identity") + 
    scale_y_log10() + facet_wrap(~dataset, scales="free_y", ncol=1, strip.position = "left") + 
    coord_flip() + ylab("Number of cells") + xlab("") +  noy + cs
  for(x in names(sces)) sces[[x]]$cluster <- unlist(cols)[paste(x, sces[[x]]$phenoid,sep=".")]
  d <- suppressWarnings(dplyr::bind_rows(lapply(cd, FUN=function(x) 
    x[,c("total_counts","total_features","cluster")]))
  )
  d$dataset <- rep(names(cd),sapply(cd,nrow))
  pf <- function(d, x) ggplot(d, aes_string(x="cluster", y=x, fill="cluster")) + 
    geom_violin() + xlab("") + coord_flip() + noy + cs +
    facet_wrap(~dataset, scales="free_y", ncol=1) + 
    theme( strip.text.x = element_blank(), 
           panel.grid.major.x=element_line(color="grey", linetype="dashed")) + 
    scale_y_continuous( trans = 'log10', breaks=10^pretty(log10(d[[x]]),n=3), 
                        labels=trans_format('log10', math_format(10^.x)) )
  rd <- function(dimred,...){
    d2 <- cbind(d, do.call(rbind, lapply(sces, FUN=function(x) reducedDim(x,dimred) )))
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

.silScale <- function(x, cols=NULL){
  if(is.null(cols)) cols <- rev(RColorBrewer::brewer.pal(n=11,"RdBu"))
  if(is.function(cols)) cols <- cols(11)
  if(length(cols)!=11) stop("`cols` should contain 11 colors.")
  bb <- c( -seq(from=sqrt(abs(min(x))), to=0, length.out=6)^2, 
           seq(from=0, to=sqrt(max(x)), length.out=6)[-1]^2 )
  circlize::colorRamp2(bb, cols)
}