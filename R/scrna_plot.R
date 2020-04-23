#' scrna_evalPlot_DR
#' 
#' Plotting aggregated evaluation results at the level of dimensionality 
#' reduction for the scRNA pipelines.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param what What to plot (default plots main metrics)
#' @param covar.type The measure of association to covariates to use, if 
#' applicable.
#' @param reorder_rows Logical; whether to sort rows (default TRUE)
#' @param reorder_columns Logical; whether to sort columns
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param scale Logical; whether to scale columns (default FALSE)
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsite of the column titles
#' @param value_cols A vector of length 2 indicating the colors of the values
#' (above and below the mean), if printed
#' @param title Plot title
#' @param anno_legend Logical; whether to include the legend for the datasets
#' @param ... Passed to `Heatmap`
#'
#' @return One or several `Heatmap` object.
#' @export
#' @examples
#' data("exampleResults", package="pipeComp")
#' scrna_evalPlot_DR(exampleResults)
#'
#' @import ComplexHeatmap grid S4Vectors
#' @importFrom viridisLite inferno
#' @importFrom stats quantile
scrna_evalPlot_DR <- function(res, what=c("auto"), 
                              covar.type=c("PC1.covar.adjR2", "meanAbsCorr.covariate2"),
                              reorder_rows=TRUE, reorder_columns=NULL,
                              agg.by=NULL, agg.fn=mean, scale=FALSE, 
                              show_heatmap_legend=FALSE, value_format="%.2f", 
                              col=NULL, col_title_fontsize=11, 
                              value_cols=c("black","white"), title=NULL,
                              anno_legend=TRUE, ...){
  pipDef <- metadata(res)$PipelineDefinition
  if(length(what)==1 && what=="auto"){
    H <- scrna_evalPlot_DR(res, "meanSilWidth", scale=FALSE, 
                           reorder_rows=reorder_rows, value_format=value_format,
                           show_heatmap_legend=TRUE, 
                           col_title_fontsize=col_title_fontsize, agg.by=agg.by, 
                           agg.fn=agg.fn, ...)
    ro <- row.names(H@matrix)
    return( H +
      scrna_evalPlot_DR(res, "log10_total_counts", scale=scale, reorder_rows=ro,
                        value_format=value_format, show_heatmap_legend=FALSE, 
                        col_title_fontsize=col_title_fontsize, agg.by=agg.by, 
                        agg.fn=agg.fn, ...) +
      scrna_evalPlot_DR(res, "total_features", scale=scale, reorder_rows=ro, 
                        value_format=value_format, show_heatmap_legend=FALSE, 
                        col_title_fontsize=col_title_fontsize, agg.by=agg.by, 
                        agg.fn=agg.fn, ...)
    )
  }
  if(length(what)>1){
    H <- scrna_evalPlot_DR(res, what[1], scale=scale, reorder_rows=reorder_rows,
                           value_format=value_format, show_heatmap_legend=FALSE,
                           col_title_fontsize=col_title_fontsize, agg.by=agg.by, 
                           agg.fn=agg.fn, ...)
    ro <- row.names(H@matrix)
    for(i in what[-1]){
      H <- H + scrna_evalPlot_DR(res, i, scale=scale, reorder_rows=reorder_rows,
                       value_format=value_format, show_heatmap_legend=FALSE,
                       col_title_fontsize=col_title_fontsize, agg.by=agg.by, 
                       agg.fn=agg.fn, ...)
    }
    return(H)
  }
  if("evaluation" %in% names(res)) res <- res$evaluation
  if("dimreduction" %in% names(res)) res <- res$dimreduction
  isSil <- FALSE
  if( what %in% 
      c("minSilWidth", "meanSilWidth", "medianSilWidth", "maxSilWidth")){
    res <- res$silhouette[[1]]
    if(is.null(title)) title <- "Subpopulation silhouette"
    isSil <- TRUE
    sname <- "silhouette\nwidth"
  }else if(what %in% 
           c("log10_total_counts","log10_total_features","total_features")){
    covar.type <- match.arg(covar.type)
    if(is.null(title))
      title <- ifelse( covar.type=="meanAbsCorr.covariate2",
                       paste0("mean(abs(corr)) with\n", what),
                       paste0("var explained by\n", what) )
    res <- res[[covar.type]]
    sname <- what
  }else if(what=="varExpl" || what=="varExpl.subpops"){
    if(is.null(title)) title <- "var explained by\nsubpopulations"
    res <- res$varExpl.subpops
    sname <- "variance\nexplained"
  }else{
    stop("Unknown metric!")
  }
  res2 <- res <- .prepRes(res, what=what, agg.by, agg.fn, pipDef=pipDef)
  row.names(res2) <- row.names(res) <- 
    gsub("resolution=", "res=", gsub("norm=norm\\.","",row.names(res2)))
  if(scale) res2 <- base::scale(res)
  res2 <- as.matrix(res2)
  if(is(reorder_rows, "Heatmap")){
    ro <- row.names(reorder_rows@matrix)
  }else{
    if(length(reorder_rows)>1){
      ro <- reorder_rows
    }else{
      if(reorder_rows){
        if(isSil){
          ro <- order( colSums(apply(res2,1,prob=c(0.05,0.5),FUN=quantile)) +
                         rowMeans(res2), decreasing=TRUE)
        }else{
          ro <- order(rowMeans(res2), decreasing=TRUE)
        }
      }else{
        ro <- seq_len(nrow(res2))
      }
    }
  }
  if(is.null(reorder_columns)) reorder_columns <- isSil
  if(reorder_columns){
    co <- order(colMeans(res), decreasing=TRUE)
  }else{
    co <- seq_len(ncol(res))
  }
  res <- res[ro,co]
  res2 <- res2[ro,co]
  if(isSil && (is.null(col) || length(col)==11) && !scale) return(
    Heatmap( res2, name=sname, cluster_rows=FALSE, col=.silScale(res2, col), 
             cluster_columns=FALSE, show_column_names = FALSE, 
             bottom_annotation=.ds_anno(colnames(res),anno_legend), 
             cell_fun=.getCellFn(res,res2,value_format,c("white","black")), 
             show_heatmap_legend=show_heatmap_legend, column_title=title, 
             column_title_gp=gpar(fontsize=col_title_fontsize), ...))

  if(is.null(col)) col <- viridisLite::inferno(100)
  Heatmap( res2, name=sname, cluster_rows=FALSE, cluster_columns=FALSE, 
           bottom_annotation=.ds_anno(colnames(res), anno_legend), 
           show_column_names=FALSE, 
           cell_fun=.getCellFn(res,res2,value_format,value_cols),
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
#' @param atTrueK Logical; whether to restrict analyses to those giving the 
#' right number of clusters
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param scale Logical; whether to scale columns (default FALSE)
#' @param value_format Format for displaying cells' values (use 
#' `value_format=""` to disable)
#' @param reorder_rows Logical; whether to sort rows (default TRUE). The row 
#' names themselves can also be passed to specify an order, or a 
#' `ComplexHeatmap`.
#' @param reorder_columns Logical; whether to sort columns
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsize of column titles.
#' @param value_cols A vector of length 2 indicating the colors of the values
#' (above and below the mean), if printed
#' @param ... Passed to `Heatmap`
#' @param title Plot title
#' @param anno_legend Logical; whether to plot the legend for the datasets
#'
#' @return One or several `Heatmap` object.
#' @export
#' @examples
#' data("exampleResults", package="pipeComp")
#' scrna_evalPlot_clust(exampleResults)
#'
#' @import ComplexHeatmap grid
#' @importFrom viridisLite inferno
#' @importFrom circlize colorRamp2
scrna_evalPlot_clust <- function(res, what="auto", atTrueK=FALSE,
                                 agg.by=NULL, agg.fn=mean, 
                                 scale=FALSE, value_format="%.2f", 
                                 reorder_rows=TRUE, reorder_columns=FALSE,
                                 show_heatmap_legend=FALSE,
                                 col=viridisLite::inferno(100), 
                                 col_title_fontsize=12,
                                 value_cols=c("black","white"), title=NULL,
                                 anno_legend=TRUE, ...){
  pipDef <- tryCatch(metadata(res)$PipelineDefinition, error=function(e) NULL)
  if(any(c("auto","delta.nbClust") %in% what))
    res$evaluation$clustering$delta.nbClust <- 
      res$evaluation$clustering$n_clus - res$evaluation$clustering$true.nbClusts
  if(is.null(agg.by)) agg.by <- .getClustAggFields(res)
  if(length(what)==1 && what=="auto"){
    H <- scrna_evalPlot_clust(res, "MI", agg.by=agg.by, atTrueK=FALSE, 
                              scale=scale, reorder_rows=TRUE, ...)
    H2 <- tryCatch(scrna_evalPlot_clust(res, "ARI", agg.by=agg.by, atTrueK=TRUE, 
                                        scale=scale, reorder_rows=H, ...),
                   error=function(e) NULL)
    if(is.null(H2) || nrow(H)!=nrow(H2)){
      H2 <- NULL
      warning("ARI at the true number of clusters is omitted, most likely ",
              "because some of the methods never yielded a clustering ",
              "with the true number of clusters...")
    }
    return( H +
      scrna_evalPlot_clust(res, "ARI", agg.by=agg.by, atTrueK=FALSE, 
                           scale=scale, reorder_rows=H, ...) + 
      scrna_evalPlot_clust(res, "min_F1", agg.by=agg.by, atTrueK=FALSE, 
                           scale=scale, reorder_rows=H, ...) + H2 + 
      scrna_evalPlot_clust(res, "delta.nbClust", agg.by=agg.by, atTrueK=FALSE,
                           col=circlize::colorRamp2(c(-5, 0, 5), 
                                                    c("red", "white", "blue")),
                           value_cols = c("black","black"), scale=FALSE, 
                           reorder_rows=H, value_format = "%.1f", ...)
    )
  }
  if(length(what)>1){
    ro <- H <- scrna_evalPlot_clust(res, what[1], agg.by=agg.by, 
                                    atTrueK=atTrueK, scale=scale, 
                                    reorder_rows=TRUE, ...)
    for(i in what[-1]) H <- H +
        scrna_evalPlot_clust(res, i, agg.by=agg.by, atTrueK=atTrueK, 
                                   scale=scale, reorder_rows=ro, ...)
    return(H)
  }
  if("evaluation" %in% names(res)) res <- res$evaluation
  if("clustering" %in% names(res)) res <- res$clustering
  if(atTrueK){
    res <- res[which(res$n_clus==res$true.nbClusts),,drop=FALSE]
    if(any(vapply(res[,agg.by], function(x) any(table(x))==0, logical(1))))
      warning("Some methods never yielded the correct number of clusters...")
    res2 <- res <- .prepRes(res, what=what, agg.by, agg.fn, pipDef=pipDef)
  }else{
    res2 <- res <- .prepRes(res, what=what, agg.by, agg.fn, pipDef=pipDef)
  }
  if(scale) res2 <- .safescale(res)
  row.names(res2) <- row.names(res) <- 
    gsub("resolution=", "res=", gsub("norm=norm\\.","",row.names(res2)))
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
  if(reorder_columns){
    co <- order(colMeans(res), decreasing=TRUE)
  }else{
    co <- seq_len(ncol(res))
  }
  res <- res[ro,co]
  res2 <- res2[ro,co]
  res2 <- as.matrix(res2)
  if(sum(!is.na(res2))<2) stop("Too few non-NA values to plot!")
  cellfn <- .getCellFn(res,res2,value_format, value_cols)
  if(is.null(title)){
    title <- gsub("_re$","\nrecall",gsub("_pr$","\nprecision",what))
    if(title=="elapsed") title <- "Computing time (s)"
    if(atTrueK) title <- paste(title, "at true \n# of clusters")
  }
  Heatmap( res2, name=paste0(what,ifelse(atTrueK,"at\ntrue K","")), 
           cluster_rows=FALSE, show_heatmap_legend=show_heatmap_legend, 
           cluster_columns=FALSE, 
           bottom_annotation=.ds_anno(colnames(res),anno_legend), 
           show_column_names = FALSE, cell_fun=cellfn, col=col,
           column_title=title, 
           column_title_gp=gpar(fontisze=col_title_fontsize), ...)
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

.renameHrows <- function(h, f=function(x) gsub("norm\\.","",x)){
  if(is(h,"HeatmapList")){
    h@ht_list <- lapply(h@ht_list, f=f, FUN=.renameHrows)
    return(h)
  }
  h@row_names_param$anno@var_env$value <- 
    f(h@row_names_param$anno@var_env$value)
  h
}
