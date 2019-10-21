scrna_evalPlot_DR <- function(res, what=c("silhouette", "covar", "covarRes", "varExpl"), elapsed=TRUE){
  what <- match.arg(what)
  if("dimreduction" %in% names(res)){
    res <- res$dimreduction
  }
  if(what=="silhouette"){
    res <- res$clust.avg.silwidth.top_10_dims
    w <- grep("^stepElapsed\\.",colnames(res))
    el <- res[,w]
    res <- res[,-w]
    w <- grep(" ",colnames(res))
    params <- res[,-w]
    res <- res[,w]
    return( Heatmap(el, name="elapsed", col=c("lightgrey","black")) + 
              Heatmap(scale(res), name="Avg.cluster\nsilhouette\nwidth") )
  }
}