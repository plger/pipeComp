#' scrna_evalPlot_silh
#' 
#' Plot a min/max/mean/median silhouette width heatmap from aggregated 
#' evaluation results of the `scrna_pipeline`.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param what What metric to plot, possible values are “minSilWidth”, 
#' “meanSilWidth” (default), “medianSilWidth”, or “maxSilWidth”.
#' @param step Name of the step for which to plot the evaluation results. 
#' Defaults to "dimreduction".
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
#' the rows. Alternatively, an expression using the columns (retained after
#' aggregation) can be passed.
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
                                 step="dimreduction", dims=1, agg.by=NULL, 
                                 agg.fn=mean, filterExpr=NULL, value_format="", 
                                 reorder_rows=FALSE, reorder_columns=TRUE,
                                 show_heatmap_legend=TRUE, 
                                 show_column_names=FALSE, 
                                 col=rev(RColorBrewer::brewer.pal(n=11,"RdBu")),
                                 font_factor=0.9, row_split=NULL, 
                                 shortNames=TRUE, value_cols=c("white","black"),
                                 title=NULL, anno_legend=TRUE, ...){
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
  if(is(res,"list") && step %in% names(res))
    res <- res[[step]]
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
  ll <- lapply(ll, FUN=function(x){
    if(any(is.na(x$N.lost))){
      w <- which(is.na(x$N.lost))
      x$N.lost[w] <- x$N.before[w]
      x$pc.lost[w] <- 100
    }
    x
  })
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
  mm$pc.lost <- 100*mm$N.lost/mm$N
  mm
}

#' scrna_evalPlot_filtering
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param steps Steps to include (default 'doublet' and 'filtering'); other 
#' steps will be averaged.
#' @param clustMetric Clustering accuracy metric to use (default `mean_F1``)
#' @param filterExpr An optional filtering expression based on the columns of 
#' the clustering evaluation (e.g. `filterExpr=param1=="value1"` or 
#' `filterExpr=n_clus==true.nbClusts`).
#' @param atNearestK Logical; whether to restrict analyses to those giving the 
#' smallest deviation from the real number of clusters (default FALSE).
#' @param returnTable Logical; whether to return the data rather than plot.
#' @param point.size Size of the points
#' @param ... passed to `geom_point`
#'
#' @return A ggplot, or a data.frame if `returnTable=TRUE`
#' @importFrom stats median
#' @import ggplot2
#' @export
#' @examples
#' data("exampleResults", package="pipeComp")
#' scrna_evalPlot_filtering(exampleResults)
scrna_evalPlot_filtering <- function(res, steps=c("doublet","filtering"), 
                                     clustMetric="mean_F1", filterExpr=TRUE,
                                     atNearestK=FALSE, returnTable=FALSE,
                                     point.size=2.2, ...){
  param_fields <- tryCatch(
    unlist(arguments(metadata(res)$PipelineDefinition)[steps]),
    error=function(e) unlist(arguments(scrna_pipeline())[steps]) )
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
  cl <- res$clustering[eval(substitute(filterExpr), res$clustering),]
  
  if(atNearestK){
    cl$absdiff <- abs(cl$n_clus-cl$true.nbClusts)
    cl <- do.call(rbind, lapply(split(cl, cl$dataset), FUN=function(x){
      cl[cl$absdiff==min(cl$absdiff),]
    }))
  }
  
  cl <- aggregate( cl[,clustMetric,drop=FALSE], 
                   by=cl[,c("dataset",param_fields)], FUN=mean )
  m <- merge(cl, x, by=c("dataset",param_fields))
  m$method <- apply( m[,param_fields,drop=FALSE], 1, collapse="+", FUN=paste )
    
  if(returnTable) return(m)
  if( length(param_fields)==2 && 
      all(sort(param_fields)==c("doubletmethod","filt")) ){
    return(ggplot(m, aes(max.lost, mean_F1, colour=filt, shape=doubletmethod))+ 
      geom_point(size=point.size, ...) + facet_wrap(~dataset, scales="free") + 
      xlab("Max proportion of subpopulation excluded") +
      labs(colour="filter set", shape="doublet removal"))
  }
  ggplot(m, aes(max.lost, mean_F1, colour=method, shape=method)) + 
    geom_point(size=point.size, ...) + facet_wrap(~dataset, scales="free") + 
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


#' scrna_evalPlot_overall
#' 
#' Plots a multi-level summary heatmap of many analyses of the `scrna_pipeline`.
#'
#' @param res Aggregated pipeline results (i.e. the output of `runPipeline` or
#' `aggregateResults`)
#' @param agg.by The paramters by which to aggregate.
#' @param width The width of individual heatmap bodies.
#' @param datasets_as_columnNames Logical; whether dataset names should be 
#' printed below the columns (except for silhouette) rather than using a
#' legend.
#' @param rowAnnoColors Optional list of colors for the row annotation variables
#' (passed to `HeatmapAnnotation(col=...)`)
#' @param column_names_gp Passed to each calls to `Heatmap`
#' @param column_title_gp Passed to each calls to `Heatmap`
#' @param heatmap_legend_param Passed to each calls to `Heatmap`
#' @param ... Passed to each calls to `Heatmap`
#'
#' @return A HeatmapList
#' @importFrom stats hclust
#' @export
#'
#' @examples
#' data("exampleResults")
#' h <- scrna_evalPlot_overall(exampleResults)
#' draw(h, heatmap_legend_side="bottom")
scrna_evalPlot_overall <- function(res, agg.by=NULL, width=NULL, 
                              datasets_as_columnNames=TRUE, 
                              rowAnnoColors=NULL,
                              column_names_gp=gpar(fontsize=10),
                              column_title_gp=gpar(fontsize=12),
                              heatmap_legend_param=list( by_row=TRUE,
                                                direction="horizontal", nrow=1), 
                              ... ){
  a <- arguments(metadata(res)$PipelineDefinition)
  if(is.null(agg.by)){
    agg.by <- c(unlist(a[-length(a)]),c("clustmethod", "dims"))
    agg.by <- agg.by[sapply(agg.by, FUN=function(x) 
                              length(unique(res$evaluation$clustering[[x]]))>1)]
  }
  agg.by <- as.character(agg.by)
  if(!all(agg.by %in% unlist(a)))
    stop("`agg.by` should be a vector of pipeline parameters.")
  
  # dimred
  sil <- res$evaluation$dimreduction$silhouette
  if(is(sil,"list")) sil <- sil[[1]]
  ll <- lapply(c("minSilWidth", "meanSilWidth"), FUN=function(x){
    .prepRes( sil, what=x, agg.by=intersect(agg.by, colnames(sil)), 
              returnParams=TRUE, shortNames=TRUE, 
              pipDef=metadata(res)$PipelineDefinition )
  })
  pp1 <- ll[[1]]$pp
  pp1$method <- NULL
  ll1 <- lapply(ll, FUN=function(x) x$res)

  # clustering
  ll <- lapply(c("ARI", "MI"), FUN=function(x){
    .prepRes( res$evaluation$clustering, what=x, returnParams=TRUE, 
              agg.by=agg.by, shortNames=TRUE, 
              pipDef=metadata(res)$PipelineDefinition )
  })
  ll[[3]] <- .prepRes(res$evaluation$clustering, what="ARI", returnParams=TRUE,
                      agg.by=agg.by, filt=expr(true.nbClusts==n_clus), 
                      shortNames=TRUE, pipDef=metadata(res)$PipelineDefinition)
  pp <- ll[[1]]$pp
  pp$method <- NULL
  ll2 <- lapply(ll, FUN=function(x){
    x <- as.data.frame(x$res)
    x <- as.matrix(colCenterScale(x[row.names(pp),]))
    row.names(x) <- row.names(pp)
    x
  })
  
  # merge dimred and clust results
  tmp <- as.character(apply( pp[,colnames(pp1),drop=FALSE], 1, 
                             collapse=" > ",FUN=paste ))
  ll1 <- lapply(ll1, FUN=function(x){
    x <- x[tmp,]
    row.names(x) <- row.names(ll2[[1]])
    x
  })
  ll <- c(ll1,ll2)

  # get max % lost
  pclost <- scrna_evalPlot_filtering(res, returnTable=TRUE)
  filt.agg.by <- intersect(agg.by,unlist(a[1:2]))
  pclost <- aggregate( pclost[,"max.lost",drop=FALSE], 
                        pclost[,c(filt.agg.by,"dataset"),drop=FALSE], FUN=mean)
  if(length(filt.agg.by)>0){
    f <- as.formula(paste(paste(filt.agg.by,collapse="+"),"~dataset"))
    pclost <- reshape2::dcast(pclost, f, value.var="max.lost")
    row.names(pclost) <- apply( pclost[,seq_along(filt.agg.by),drop=FALSE], 1, 
                                collapse=" > ",FUN=paste )
    tmp <- as.character(apply( pp[,filt.agg.by,drop=FALSE], 1, collapse=" > ",
                               FUN=paste ))
    pclost <- pclost[tmp,setdiff(colnames(pclost), filt.agg.by),]
    row.names(pclost) <- row.names(ll2[[1]])
  }else{
    pclost <- matrix(pclost$max.lost, nrow=nrow(ll2[[1]]), ncol=nrow(pclost), 
                     byrow=TRUE, 
                     dimnames=list(row.names(ll2[[1]])), pclost$dataset)
  }
  pclost <- apply(pclost,1,FUN=max)
  
  ll2 <- list( list(mat=ll[[1]], title="min silh.\nwidth", 
                    cluster_columns=TRUE, name="silhouette width"),
               list(mat=ll[[2]], title="mean silh.\nwidth", 
                    cluster_columns=TRUE, show_heatmap_legend=FALSE),
               list(mat=ll[[3]], title="mean ARI", name="ARI (MADs)", 
                    cluster_columns=FALSE),
               list(mat=ll[[4]], title="mean MI", name="MI (MADs)", 
                    cluster_columns=FALSE),
               list(mat=ll[[5]], title="mean ARI\nat true k", 
                    name="ARI at true k (MADs)", cluster_columns=FALSE)
               )

  if("doubletmethod" %in% colnames(pp))
    pp$doubletmethod <- gsub("^doublet\\.","",pp$doubletmethod)
  if("clustmethod" %in% colnames(pp))
    pp$clustmethod <- gsub("^clust\\.","",pp$doubletmethod)
  for(f in c("filt","sel","norm")){
    if(f %in% colnames(pp)) pp[[f]] <- gsub(paste0("^",f,"\\."),"",pp[[f]])
  }
  if(is.null(rowAnnoColors)){
    ha <- HeatmapAnnotation( which="row",
                             "max\n% lost"=anno_barplot( pclost, bar_width=0.85, 
                                                         border=FALSE,
                                                         width=unit(1.5,"cm"),
                                                         gp=gpar(fill="#282828",
                                                                 col="#282828") ),
                             df=pp, annotation_legend_param=list("side"="right") )
  } else {
    ha <- HeatmapAnnotation( which="row", col=rowAnnoColors,
                             "max\n% lost"=anno_barplot( pclost, bar_width=0.85, 
                                                         border=FALSE,
                                                         width=unit(1.5,"cm"),
                                                         gp=gpar(fill="#282828",
                                                                 col="#282828") ),
                             df=pp, annotation_legend_param=list("side"="right") )
  }
  
  h <- hclust(dist(do.call(cbind, ll)))
  silhscale <- .silScale(cbind(ll2[[1]]$mat, ll2[[2]]$mat))
  H <- NULL
  for(i in seq_along(ll2)){
    if(i==1){
      hi <- h
    }else{
      hi <- FALSE
    }
    if(grepl("silh", ll2[[i]]$title)){
      col <- silhscale
      scn <- FALSE
    }else{
      col <- .defaultColorMapping(ll2[[i]]$mat, center=TRUE)
      scn <- datasets_as_columnNames
    }
    la <- ra <- NULL
    if(i==length(ll2)) ra <- ha
    ba <- .ds_anno(colnames(ll2[[i]]$mat), 
                   legend=(!datasets_as_columnNames && i==1))
    name <- ifelse(is.null(ll2[[i]]$name),ll2[[i]]$title,ll2[[i]]$name)
    wi <- ifelse(is.null(width), unit(ifelse(i<=2,3.2,2.5), "cm"), width)
    H <- H + Heatmap( ll2[[i]]$mat, name=name, cluster_rows=hi, 
                      show_row_names=FALSE, show_column_names=scn, 
                      heatmap_legend_param=heatmap_legend_param, 
                      column_title=ll2[[i]]$title, col=col, width=wi, 
                      cluster_columns=ll2[[i]]$cluster_columns, 
                      show_column_dend=FALSE, bottom_annotation=ba, 
                      left_annotation=la, right_annotation=ra, 
                      column_names_gp=column_names_gp,
                      column_title_gp=column_title_gp,
                      show_heatmap_legend=ifelse(
                        is.null(ll2[[i]]$show_heatmap_legend),TRUE,
                        ll2[[i]]$show_heatmap_legend), ... )
  }
  H
}