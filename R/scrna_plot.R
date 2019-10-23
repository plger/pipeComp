scrna_evalPlot_DR <- function(res, what=c("auto","silhouette", "covar", "covarRes", "varExpl", "elapsed"),
                              covar=c("log10_total_counts","log10_total_features",
                                      "total_features"), reorder_rows=TRUE,
                              agg.by=NULL, agg.fn=mean, scale=FALSE, 
                              show_heatmap_legend=FALSE, value_format="%.2f", 
                              col=viridisLite::inferno(100), ...){
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
      scrna_evalPlot_DR(res, "silhouette", scale=scale, reorder_rows=FALSE, show_heatmap_legend=FALSE) + 
      scrna_evalPlot_DR(res, "varExpl", scale=scale, reorder_rows=FALSE, show_heatmap_legend=FALSE) + 
      scrna_evalPlot_DR(res, "covarRes", scale=scale, reorder_rows=FALSE, show_heatmap_legend=FALSE) +
      scrna_evalPlot_DR(res, "covarRes", covar="total_features", scale=scale, reorder_rows=FALSE, show_heatmap_legend=FALSE)
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
  cellfn <- .getCellFn(res,res2,value_format)
  name <- switch( what,
                 silhouette="silhouette\nwidth",
                 covar=covar,
                 covarRes=covar,
                 varExpl="Variance\nexplained",
                 elapsed="Running\ntime (s)",
                 NULL
  )
  title <- switch( what,
                 silhouette="Average silhouette width\nper subpopluation",
                 covar=paste0("Correlation with\n",covar),
                 covarRes=paste0("Residual corr with\n",covar),
                 varExpl="Variance explained\nby subpopulations",
                 elapsed="Running time (s)",
                 NULL
  )
  Heatmap( res2, name=name, cluster_rows=FALSE,
           cluster_columns=FALSE, bottom_annotation=.ds_anno(colnames(res)), 
           show_column_names = FALSE, cell_fun=cellfn, col=col,
           show_heatmap_legend=show_heatmap_legend,
           column_title=title, column_title_gp=gpar(fontisze=10), ...)
}


scrna_evalPlot_clust <- function(res, what="auto", agg.by=NULL, agg.fn=mean, 
                                 scale=FALSE, value_format="%.2f", reorder_rows=TRUE, 
                                 show_heatmap_legend=FALSE, 
                                  col=viridisLite::inferno(100), ...){
  if("clustering" %in% names(res)) res <- res$clustering
  what_options <- sapply(strsplit(colnames(res)[grep(" ",colnames(res))]," "),FUN=function(x) x[[2]])
  what <- match.arg(what, c(what_options, "elapsed", "auto"))
  if(what=="auto"){
    return( 
      scrna_evalPlot_clust(res, "ARI", scale=scale, reorder_rows=FALSE, ...) + 
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
  Heatmap( res2, name=what, cluster_rows=FALSE, show_heatmap_legend=show_heatmap_legend, 
           cluster_columns=FALSE, bottom_annotation=.ds_anno(colnames(res)), 
           show_column_names = FALSE, cell_fun=cellfn, col=col,
           column_title=title, column_title_gp=gpar(fontisze=10), ...)
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
  pp <- res[,grep("^stepElapsed\\.",colnames(res),invert=TRUE)]
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
.getCellFn  <- function(res,res2,value_format){
  resmid <- range(res2, na.rm=TRUE)
  resmid <- resmid[1]+(resmid[2]-resmid[1])/2
  function(j, i, x, y, width, height, fill){
    if(value_format=="" || is.null(value_format) || is.na(value_format)) return(NULL)
    lab <- sprintf(value_format, res[i,j])
    if(!any(abs(res)>1)){
      lab <- gsub("^0\\.",".",lab)
      lab <- gsub("^-0\\.","-.",lab)
    } 
    grid.text(lab, x, y, gp = gpar(fontsize = 10, col=ifelse(res2[i,j]>resmid,"black","white")))
  }
}

#' describeDatasets
#' 
#' Plots descriptive information about the datasets
#'
#' @param sces A character vector of paths to SCE rds files, or a list of SCEs
#' @param ridges Logical; whether to plot ridges
#'
#' @return A plot_grid output
#'
#' @export
#' @import SingleCellExperiment scran ggplot2 ggridges
#' @importFrom scater plotReducedDim logNormCounts
#' @importFrom cowplot plot_grid
describeDatasets <- function(sces, ridges=TRUE){
  options(scipen=10000)
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
  cd <- lapply(sces, FUN=function(x) as.data.frame(colData(x)))
  cols <- lapply(cd, FUN=function(x){ 
    x <- as.character(unique(x$phenoid))
    y <- GTscripts::getQualitativePalette(length(x))
    names(y) <- x
    y
  })
  noaxes <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank())
  names(cols) <- names(cd)
  if(ridges){
    pf <- function(d, x) ggplot(cd[[d]], aes_string(x=x, y="phenoid", fill="phenoid")) + geom_density_ridges() + scale_x_log10() + ylab("") + scale_fill_manual(values=cols[[d]]) + theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_blank())
  }else{
    pf <- function(d, x) ggplot(cd[[d]], aes_string(x="phenoid", y=x, fill="phenoid")) + geom_violin() + scale_y_log10() + xlab("") + scale_fill_manual(values=cols[[d]]) + coord_flip() + theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_blank())
  }
  nplots <- 5
  pl <- suppressWarnings(suppressMessages({ c( 
    lapply(names(cd), FUN=function(x){
      tmp <- as.data.frame(table(cd[[x]]$phenoid))
      tmp$label <- paste0(tmp$Var1, " (", tmp$Freq, ")")
      ggplot(cd[[x]], aes(phenoid, fill=phenoid)) + geom_histogram(stat="count", show.legend=FALSE) + 
        scale_y_log10() + coord_flip() + scale_fill_manual(values=cols[[x]]) + ylab("Number of cells") + 
        xlab(x) + geom_text(data=tmp, mapping=aes(Var1, y=1, label=label), inherit.aes = FALSE, nudge_y = 1) + 
        theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
              axis.line.x = element_blank(), axis.title.y=element_text(size=12,face="bold"))
    }),
    lapply(names(cd), FUN=function(x) pf(x, "total_counts")),
    lapply(names(cd), FUN=function(x) pf(x, "total_features")),
    lapply(names(sces), FUN=function(x) plotReducedDim(sces[[x]], "tSNE", colour_by ="phenoid") + scale_fill_manual(values=cols[[x]]) + ggtitle(paste(x,"tSNE")) + noaxes),
    lapply(names(sces), FUN=function(x) plotReducedDim(sces[[x]], "umap", colour_by ="phenoid") + scale_fill_manual(values=cols[[x]]) + ggtitle(paste(x,"UMAP")) + noaxes)
  )}))
  pl <- pl[as.numeric(sapply(seq_along(cd), FUN=function(x) 1:nplots*length(cd)-length(cd)+x))]
  pl <- lapply(pl, FUN=function(x) x + theme(legend.position = "none", aspect.ratio=1, axis.text=element_text(size=10), axis.title=element_text(size=11)))
  suppressMessages(plot_grid(plotlist=pl, nrow = length(cd), ncol=nplots))
}