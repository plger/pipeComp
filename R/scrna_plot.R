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
#' @param reorder_columns Logical; whether to sort columns
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
#' @import ComplexHeatmap matrixStats
#' @importFrom viridisLite inferno
scrna_evalPlot_DR <- function(res, what=c("auto","silhouette", "covar", "covarRes", "varExpl", "elapsed"),
                              covar=c("log10_total_counts","log10_total_features",
                                      "total_features"), reorder_rows=TRUE,
                              reorder_columns=what=="silhouette",
                              agg.by=NULL, agg.fn=mean, scale=FALSE, 
                              show_heatmap_legend=FALSE, value_format="%.2f", 
                              col=NULL, col_title_fontsize=11, 
                              anno_legend=TRUE, ...){
  what <- match.arg(what)
  covar <- covar[1]
  if("dimreduction" %in% names(res)) res <- res$dimreduction
  if(what=="auto"){
    H <- scrna_evalPlot_DR(res, "silhouette", scale=FALSE, reorder_rows=reorder_rows, value_format=value_format, show_heatmap_legend=TRUE, col_title_fontsize=col_title_fontsize, agg.by=agg.by, agg.fn=agg.fn, ...)
    ro <- row.names(H@matrix)
    return( H +
      scrna_evalPlot_DR(res, "varExpl", scale=scale, reorder_rows=ro, value_format=value_format, show_heatmap_legend=TRUE, col_title_fontsize=col_title_fontsize, agg.by=agg.by, agg.fn=agg.fn, ...) + 
      scrna_evalPlot_DR(res, "covarRes", scale=scale, reorder_rows=ro, value_format=value_format, show_heatmap_legend=FALSE, col_title_fontsize=col_title_fontsize, agg.by=agg.by, agg.fn=agg.fn, ...) +
      scrna_evalPlot_DR(res, "covarRes", covar="total_features", scale=scale, reorder_rows=ro, value_format=value_format, show_heatmap_legend=FALSE, col_title_fontsize=col_title_fontsize, agg.by=agg.by, agg.fn=agg.fn, ...)
    )
  }
  el <- grep("^stepElapsed\\.", colnames(res[[1]]), value=TRUE)
  res <- switch( what,
                 silhouette=res$clust.avg.silwidth,
                 covar=res[[paste0("PC1_covar.",covar)]],
                 covarRes=res[[paste0("PC1_covarR.",covar)]],
                 varExpl=res$PCtop5.R2,
                 elapsed=res$PCtop5.R2[,el],
                 stop("Unknown plot type requested")
                 )
  res2 <- res <- .prepRes(res, agg.by, agg.fn, elapsed=what=="elapsed")
  if(scale) res2 <- base::scale(res)
  res2 <- as.matrix(res2)
  if(is(reorder_rows, "Heatmap")){
    ro <- row.names(reorder_rows@matrix)
  }else{
    if(length(reorder_rows)>1){
      ro <- reorder_rows
    }else{
      if(reorder_rows){
        if(what=="silhouette"){
          ro <- order(rowMedians(res2)+rowMeans(res2)+rowMins(res2), decreasing=TRUE)
        }else{
          ro <- order(rowMeans(res2), decreasing=TRUE)
        }
      }else{
        ro <- seq_len(nrow(res2))
      }
    }
  }
  if(reorder_columns){
    co <- order(colMeans(res), decreasing=TRUE)
  }else{
    co <- 1:ncol(res)
  }
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
             cluster_columns=FALSE, show_column_names = FALSE, 
             bottom_annotation=.ds_anno(colnames(res),anno_legend), 
             cell_fun=.getCellFn(res,res2,value_format,c("white","black")), 
             show_heatmap_legend=show_heatmap_legend, column_title=title, 
             column_title_gp=gpar(fontsize=col_title_fontsize), ...))

  if(is.null(col)) col <- viridisLite::inferno(100)
  Heatmap( res2, name=name, cluster_rows=FALSE, cluster_columns=FALSE, 
           bottom_annotation=.ds_anno(colnames(res), anno_legend), 
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
#' @param reorder_columns Logical; whether to sort columns
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsize of column titles.
#' @param ... Passed to `Heatmap`
#'
#' @return One or several `Heatmap` object.
#' @export
#'
#' @import ComplexHeatmap matrixStats
#' @importFrom viridisLite inferno
#' @importFrom grid gpar
scrna_evalPlot_clust <- function(res, what="auto", agg.by=NULL, agg.fn=mean, 
                                 scale=FALSE, value_format="%.2f", 
                                 reorder_rows=TRUE, reorder_columns=FALSE,
                                 show_heatmap_legend=FALSE,
                                 col=viridisLite::inferno(100), 
                                 col_title_fontsize=12, anno_legend=TRUE, ...){
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
    co <- 1:ncol(res)
  }
  res <- res[ro,co]
  res2 <- res2[ro,co]
  res2 <- as.matrix(res2)
  cellfn <- .getCellFn(res,res2,value_format)
  title <- switch( what,
                   elapsed="Running time (s)",
                   gsub("_re$","\nrecall",gsub("_pr$","\nprecision",what))
  )
  row.names(res2) <- gsub("resolution=", "res=", gsub("norm=norm.","",row.names(res2),fixed=TRUE))
  Heatmap( res2, name=what, cluster_rows=FALSE, show_heatmap_legend=show_heatmap_legend, 
           cluster_columns=FALSE, bottom_annotation=.ds_anno(colnames(res),anno_legend), 
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
#' @param reorder_columns Logical; whether to sort columns
#' @param show_heatmap_legend Passed to `Heatmap`
#' @param col Colors for the heatmap
#' @param col_title_fontsize Fontsize of the column titles
#' @param ... Passed to `Heatmap`
#'
#' @return A `Heatmap` object.
#' @export
#'
#' @import ComplexHeatmap matrixStats
#' @importFrom viridisLite inferno
#' @importFrom grid gpar
scrna_evalPlot_clustAtTrueK <- function(res, what="ARI", agg.by=NULL, 
                                        agg.fn=mean, scale=FALSE, 
                                        value_format="%.2f", reorder_rows=TRUE, 
                                        reorder_columns=FALSE,
                                        show_heatmap_legend=FALSE,
                                        col=viridisLite::inferno(100), 
                                        col_title_fontsize=11, title=NULL,
                                        name=NULL, anno_legend=TRUE, ...){
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
        res <- res[x,]
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
  if(reorder_columns){
    co <- order(colMeans(res, na.rm=TRUE), decreasing=TRUE)
  }else{
    co <- 1:ncol(res)
  }
  res <- res[ro,co]
  res2 <- res2[ro,co]
  cellfn <- .getCellFn(res,res2,value_format)
  if(is.null(title)) title <- paste(what, "at true\n# of clusters")
  row.names(res2) <- gsub("resolution=", "res=", gsub("norm=norm.","",row.names(res2),fixed=TRUE))
  Heatmap( res2, name=ifelse(is.null(name),what,name), cluster_rows=FALSE, 
           show_heatmap_legend=show_heatmap_legend, cluster_columns=FALSE, 
           bottom_annotation=.ds_anno(colnames(res),anno_legend), 
           show_column_names = FALSE, cell_fun=cellfn, col=col,
           column_title=title, column_title_gp=gpar(fontisze=col_title_fontsize), ...)
}

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


.ds_anno <- function(x, legend=TRUE){
  y <- sapply(strsplit(gsub("stepElapsed\\.","",x)," "),FUN=function(x) x[[1]])
  cols <- getQualitativePalette(length(unique(y)))
  names(cols) <- sort(unique(y))
  ComplexHeatmap::HeatmapAnnotation(dataset=y, col=list(dataset=cols), show_annotation_name = FALSE, show_legend=legend)
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
    if(!any(abs(res)>1, na.rm=TRUE)){
      lab <- gsub("^0\\.",".",lab)
      lab <- gsub("^-0\\.","-.",lab)
    } 
    grid.text(lab, x, y, gp = gpar(fontsize = 10, col=ifelse(res2[i,j]>resmid,cols[1],cols[2])))
  }
}

.mergeFilterOut <- function(ll){
  if(length(ll)==1) return(ll[[1]])
  # get back the initial number of cells
  N <- attr(ll[[1]],"initialNumbers")
  if(is.null(N)){
    x <- abs(as.matrix(ll[[1]]))
    pc <- x[,grep("^pcOut",colnames(x))]
    n <- x[,grep("^nOut",colnames(x))]
    N <- matrixStats::colMedians(round(100*n/pc), na.rm = TRUE)
    N[is.na(N)] <- 100
  }
  ll2 <- lapply(ll, FUN=function(x){
    x <- x[,grep("^nOut",colnames(x))]
    attr(x, "initialNumbers") <- NULL
    cbind(parsePipNames(row.names(x)), x)
  })
  x <- rev(ll2)[[1]]
  for(i in seq.int(length(ll)-1,1)){
    y <- ll2[[i]]
    of <- intersect(colnames(x),colnames(y))
    of <- of[grep("^nOut", of, invert=TRUE)]
    p1 <- apply(x[,of,drop=FALSE],1,FUN=function(x) 
      paste(paste0(of,"=",x), collapse=";"))
    y <- y[,grep("^nOut", colnames(y))]
    row.names(y) <- apply(ll2[[i]][,of,drop=FALSE],1,FUN=function(x) 
      paste(paste0(of,"=",x), collapse=";"))
    colnames(y) <- NULL
    x <- x[,grep("^nOut", colnames(x))]
    x+y[p1,]
  }
  pc <- round(100*t(t(x)/N),3)
  colnames(pc) <- gsub("^nOut","pcOut", colnames(pc))
  cbind(x, pc, deparse.level=0)
}

scrna_evalPlot_filtering <- function(res, steps=c("doublet","filtering"), returnTable=FALSE){
  co <- .mergeFilterOut(res[steps])
  co <- abs(as.matrix(co[,grep("^pcOut",colnames(co))]))
  ds <- sapply(strsplit(colnames(co),"\\."), FUN=function(x) x[2])
  ds <- split(seq_len(ncol(co)), ds)
  df <- data.frame(
    method=rep(row.names(co), length(ds)),
    dataset=rep(names(ds), each=nrow(co)),
    maxPCout=unlist(lapply(ds, FUN=function(x) rowMaxs(co[,x]))),
    medianPCout=unlist(lapply(ds, FUN=function(x) rowMedians(co[,x])))
  )
  ds <- sapply(strsplit(colnames(res$clustering)," "), FUN=function(x) x[1])
  ds <- split(seq_len(ncol(res$clustering)), ds)
  ds <- ds[grep("stepElapsed",names(ds),invert=TRUE)]
  ds <- lapply(ds, FUN=function(x){
    cl <- res$clustering[,x]
    by <- sapply(strsplit(row.names(cl),";"), FUN=function(x) paste(x[1:2],collapse=";"))
    meanF1 <- aggregate(cl[[grep("mean_F1",colnames(cl))]], by=list(by), na.rm=TRUE, FUN=mean)
    topF1 <- aggregate(cl[[grep("mean_F1",colnames(cl))]], by=list(by), na.rm=TRUE, FUN=max)
    w <- which(cl[[grep("n_clus",colnames(cl))]]==cl[[grep("nbClusts",colnames(cl))]])
    cl <- cl[w,]
    tmp <- aggregate(cl[[grep("mean_F1",colnames(cl))]], by=list(by[w]), na.rm=TRUE, FUN=mean)
    F1atK <- tmp[,2]
    names(F1atK) <- tmp[,1]
    d <- data.frame(method=meanF1[,1], meanF1=as.numeric(meanF1[,2]), 
                    topF1=as.numeric(topF1[,2]), F1atK=F1atK[as.character(meanF1[,1])])
    d
  })
  ds2 <- bind_rows(ds,.id="dataset")
  m <- merge(df, ds2, by=c("dataset","method"))
  if(returnTable) return(m)
  ggplot(m, aes(maxPCout, meanF1, group=method, colour=method)) + geom_point(size=3) + 
    facet_wrap(~dataset, scales="free") + xlab("Proportion") + 
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
#' @import SummarizedExperiment SingleCellExperiment scran scater ggplot2 scales
#' @importFrom cowplot plot_grid
scrna_describeDatasets <- function(sces, pt.size=0.3, ...){
  if(is.null(names(sces))) names(sces) <- paste0("dataset",seq_along(sces))
  if(is.character(sces)){
    sces <- lapply(sces, FUN=function(x){
      sce <- readRDS(x)
      sce <- logNormCounts(sce)
      var.stats <- modelGeneVar(sce)
      sce <- denoisePCA(sce, technical=var.stats)
      reducedDim(sce, "tSNE") <- Rtsne::Rtsne(reducedDim(sce))$Y
      reducedDim(sce, "umap") <- uwot::umap(reducedDim(sce))
      sce
    })
  }
  tt <- lapply(sces, FUN=function(x) table(x$phenoid))
  cols <- lapply(tt, FUN=function(x){
    y <- getQualitativePalette(length(x))
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
  bb <- c( -seq(from=sqrt(abs(min(x, na.rm=TRUE))), to=0, length.out=6)^2, 
           seq(from=0, to=sqrt(max(x, na.rm=TRUE)), length.out=6)[-1]^2 )
  circlize::colorRamp2(bb, cols)
}

.renameHrows <- function(h, f=function(x) gsub("norm\\.","",x)){
  if(is(h,"HeatmapList")){
    h@ht_list <- lapply(h@ht_list, f=f, FUN=.renameHrows)
    return(h)
  }
  h@row_names_param$anno@var_env$value <- f(h@row_names_param$anno@var_env$value)
  h
}
