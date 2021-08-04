#' knnPrecision
#' 
#' Calculates the precision-at-k (PRK) curve, i.e. the proportion of correct k 
#' nearest neighbors across values of k, as a measure of the accuracy of a 
#' kNN network judged against a 'true' one.
#'
#' @param pred,truth The matrices of predicted and true values, respectively, 
#' with variables as columns. Alternatively, if `isKNN=TRUE`, matrices of 
#' nearest neighbor indices (the first column being the first neighbor or each 
#' item), as produced by \code{\link[BiocNeighbors]{findKNN}$index}. Either 
#' way, the rows are assumed to correspond between the two objects.
#' @param kmax The maximum k nearest neighbor to consider.
#' @param isKNN Logical; whether `pred` and `truth` are already computed kNNs.
#' If FALSE (default), both kNN are first computed.
#'
#' @return A list containing 1) a matrix of kNN precisions (as %) for each 
#' item, at each value of k, and 2) the area under the PRK curve for each item.
#' @export
#' @importFrom BiocNeighbors findKNN AnnoyParam
#'
#' @examples
#' x <- runif(500)*10
#' e1 <- matrix(rpois(length(x),x),ncol=10)
#' e2 <- matrix(rpois(length(x),x),ncol=10)
#' ac <- knnPrecision(e1,e2)
#' hist(ac$AUC, xlab="Area under the PRK curve")
#' boxplot(as.data.frame(ac$precisions), xlab="Number of nearest neighbors", 
#'         ylab="Proportion correct", outline=FALSE)
knnPrecision <- function(pred, truth, kmax=NULL, isKNN=FALSE){
  stopifnot(!is.null(dim(pred)) && !is.null(dim(truth)))
  stopifnot(nrow(pred)==nrow(truth))
  if(!isKNN){
    if(is.null(kmax)){
      kmax <- min(150,ceiling(nrow(pred)/2))
      message("Using maximum k of ", kmax)
    }
    stopifnot(nrow(pred)>kmax)
    rn <- row.names(pred)
    pred <- findKNN(pred, k=kmax, BNPARAM=AnnoyParam(), 
                               get.distance=FALSE)$index
    truth <- findKNN(truth, k=kmax, BNPARAM=AnnoyParam(), 
                               get.distance=FALSE)$index
    row.names(truth) <- row.names(pred) <- rn
  }else{
    if(is.null(kmax)) kmax <- min(150,ncol(pred))
  }
  stopifnot(dim(pred)==dim(truth))
  stopifnot(ncol(pred)>=kmax)
  
  names(o) <- o <- seq_len(kmax)
  o <- sapply(o, FUN=function(n){
    sapply(seq_len(nrow(pred)), FUN=function(x){
      if(n<1) return(0L)
      as.integer(100*sum(pred[x,seq_len(n)] %in% truth[x,seq_len(n)])/n)
    })
  })
  row.names(o) <- row.names(pred)
  list(precisions=o, AUC=rowSums(o)/kmax)
}
