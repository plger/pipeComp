#' evaluateDEA
#' 
#' Evaluates a differential expression analysis (DEA).
#'
#' @param dea Expects a data.frame with logFC and FDR, as produced by 
#' `edgeR::topTags`, `limma::topTable` or `DESeq2::results`.
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
  ## we make sure that the column names of `dea` are standard:
  dea <- pipeComp:::.homogenizeDEA(dea)
  ## within Pipecomp, the truth should be passed along with the `dea` object, 
  ## so we retrieve it here:
  if(is.null(truth)) truth <- metadata(dea)$truth
  dea <- cbind(dea, truth[row.names(dea),])
  ## we get rid of genes for which the truth is unknown:
  dea <- dea[!is.na(dea$expected.beta),]
  ## comparison of estimated and expected log2 folchanges:
  res <- c(logFC.pearson=cor(dea$logFC, dea$expected.beta, use = "pairwise"),
           logFC.spearman=cor(dea$logFC, dea$expected.beta, 
                              use = "pairwise", method="spearman"),
           logFC.mad=median(abs(dea$logFC-dea$expected.beta),na.rm=TRUE),
           ntested=sum(!is.na(dea$PValue) & !is.na(dea$FDR)))
  ## evaluation of singificance calls
  names(th) <- th
  res2 <- t(vapply( th, FUN.VALUE=vector(mode="numeric", length=6), 
                    FUN=function(x){
            ## for each significance threshold, calculate the various metrics
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

#' dea_evalPlot_curve
#'
#' @param res Aggregated results of the DEA pipeline
#' @param scales Passed to `facet_grid`
#' @param agg.by Aggregate results by these columns (default no aggregation)
#' @param agg.fn Function for aggregation (default mean)
#' @param xlim Optional vector of x limits
#' @param colourBy Name of column by which to colour
#' @param shapeBy Name of column determining the shape of the points. If 
#' omitted, the shape will indicate whether the nominal FDR is below or equal 
#' the real FDR.
#' @param pointsize Size of the points
#'
#' @return A ggplot.
#' @export
#' @examples
#' data("exampleDEAresults", package="pipeComp")
#' dea_evalPlot_curve(exampleDEAresults, agg.by=c("sva.method"))
dea_evalPlot_curve <- function(res, scales="free", agg.by=NULL, agg.fn=mean, 
                               xlim=c(NA,NA), colourBy="method", shapeBy=NULL,
                               pointsize=4){
  pd <- NULL
  if(is(res,"SimpleList")) pd <- metadata(res)$PipelineDefinition
  if(is(res,"SimpleList") && "evaluation" %in% names(res)) res <- res$evaluation
  if(is(res,"list") && "dea" %in% names(res)) res <- res$dea
  if(is(res,"list") && "significance" %in% names(res)) res <- res$significance
  if(is.null(agg.by)) agg.by <- intersect(colnames(res), unlist(arguments(pd)))
  d <- aggregate(res[,c("FDR","TPR")], by=res[,c("dataset","threshold",agg.by)], na.rm=TRUE, FUN=agg.fn)
  pp <- d[,agg.by,drop=FALSE]
  d$method <- apply(pp, 1, collapse=" > ", FUN=paste)
  if(all(c("filt","minCount") %in% colnames(d))){
    d$filter <- paste(d$filt, "(",d$minCount,")")
    d$filter[which(d$filt=="none")] <- "none"
  } 
  p <- ggplot(d, aes_string("FDR", "TPR", group="method", colour=colourBy, 
                            shape=shapeBy)) + 
    geom_vline(xintercept=unique(d$threshold), linetype="dashed", 
               colour="darkgrey") + 
    geom_line(size=1) + geom_point(size=pointsize)
  if(is.null(shapeBy)) p <- p + 
    geom_point(data=d[d$FDR>d$threshold,], size=3, colour="white")
  p + facet_wrap(~dataset, scales=scales) + 
    scale_x_sqrt( breaks=c(0.001,0.01,0.05,0.1,0.2,0.4), 
                  limits=xlim ) + 
    theme(axis.text.x=element_text(angle = 90, hjust = 1))
}
