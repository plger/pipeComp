#' evaluateClustering
#'
#' Evaluates a clustering using 'true' labels
#'
#' @param x The clustering labels
#' @param tl The true labels
#'
#' @return A numeric vector.
#' @importFrom aricode clustComp
#' @export
evaluateClustering <- function(x, tl){
  e <- match_evaluate_multiple(x, tl)
  x <- as.character(x)
  unmatched <- length(x)-sum(e$n_cells_matched)
  c( unlist(e), unmatched.cells=unmatched, 
     unlist(aricode::clustComp(x,tl)) )
}

.compileExcludedCells <- function(before, after){
  tt <- table(before$phenoid)
  tt <- rbind( tt, table(after$phenoid)[names(tt)])
  row.names(tt) <- c("before","after")
  list( table=tt, excluded=setdiff(colnames(before), colnames(after)) )
}

.aggregateExcludedCells <- function(res){
  cls <- sapply(res, FUN=function(x) ncol(x[[1]]$table))
  nout <- do.call(cbind, lapply(names(res), FUN=function(ds){ 
    t(sapply(res[[ds]], FUN=function(x){
      x <- x$table
      y <- as.numeric(x[2,])-as.numeric(x[1,])
      names(y) <- paste0(ds,".",colnames(x))
      y
    }))
  }))
  pc <- do.call(cbind, lapply(names(res), FUN=function(ds){ 
    t(sapply(res[[ds]], FUN=function(x){
      x <- x$table
      y <- (as.numeric(x[2,])-as.numeric(x[1,]))/as.numeric(colSums(x))
      names(y) <- paste0(ds,".",colnames(x))
      round(y*100,2)
    }))
  }))
  colnames(pc) <- paste0("pcOut.", colnames(pc))
  colnames(nout) <- paste0("nOut.", colnames(pc))
  return(cbind(nout, pc))
}

#' @importFrom matrixStats rowMins
#' @importFrom dplyr bind_cols bind_rows
.aggregateClusterEvaluation <- function(res){
  res <- lapply(res, FUN=function(a){
    x <- as.data.frame(t(dplyr::bind_rows(a)))
    colnames(x) <- names(a[[1]])
    for(f in c("pr","re","F1")){
      w <- grep(paste0("^",f,"\\."),colnames(x))
      x[[paste0("min_",f)]] <- matrixStats::rowMins(as.matrix(x[,w,drop=FALSE]))
    }
    y <- as.data.frame(x[,grep("\\.",colnames(x),invert=TRUE)])
    y$true.nbClusts <- length(grep("^pr\\.", colnames(x)))
    y
  })
  for(i in names(res)) colnames(res[[i]]) <- paste(i, colnames(res[[i]]))
  as.data.frame(bind_cols(res), row.names=row.names(res[[1]]))
}

.getVE <- function(x, cl){
  tryCatch({
    summary(lm(x~cl))$r.squared
  }, error=function(e){ print(e); return(NA) })
}



#' evaluateDimRed
#'
#' Gathers evaluation statistics on a reduced space using known cell labels.
#'
#' @param x The matrix of the reduced space, with cells as rows and components 
#' as columns
#' @param clusters The vector indicating each cell's cluster.
#' @param n A numeric vector indiciating the number of top dimensions at which 
#' to gather statistics (default `c(10,20,50)`). Will use all available dimensions
#' if a higher number is given.
#'
#' @return A list with the following components:
#' * silhouettes: a matrix of the silhouette for each cell-cluster pair at each 
#' value of `n`
#' * clust.avg.silwidth: a matrix of the cluster average width at each value of `n`
#' * R2: the proportion of variance in each component (up to `max(n)`) that is 
#' explained by the clusters (i.e. R-squared of a linear model).
#' 
#' @export
evaluateDimRed <- function(x, clusters=NULL, n=c(10,20,50), covars=NULL){
  if(is(x,"Seurat")){
    if(is.null(covars)) covars <- x@meta.data[,c("log10_total_features", "log10_total_counts", "total_features")]
    if(is.null(clusters)) clusters <- x$phenoid
    x <- x@reductions$pca@cell.embeddings
  }else if(is(x, "SingleCellExperiment")){
    if(is.null(covars)) covars <- as.data.frame(colData(x[,c("log10_total_features", "log10_total_counts", "total_features")]))
    if(is.null(clusters)) clusters <- x$phenoid
    x <- reducedDim(x, "PCA")
  }else{
    if(is.null(clusters)) stop("`clusters` must be given!")
  }
  library(cluster)
  clusters <- as.factor(clusters)
  n <- unique(sapply(n, y=ncol(x), FUN=function(x,y) min(x,y)))
  si <- lapply(n, FUN=function(dims){
    silhouette(as.integer(clusters), dist(x[,1:dims]))
  })
  names(si) <- sapply(n, FUN=function(i){
    if(i==ncol(x)) return(paste0("all_", i,"_dims"))
    return(paste0("top_", i,"_dims"))
  })
  # summarize silhouette information
  if(length(n)==1){
    silhouettes <- si[[1]][,1:3]
  }else{
    silhouettes <- cbind(si[[1]][,1:2], do.call(cbind, lapply(si,FUN=function(x) as.numeric(x[,3]))))
  }
  
  # cluster average silhouette width
  csw <- t(sapply(si, FUN=function(x) clus.avg.widths=summary(x)$clus.avg.widths))
  colnames(csw) <- levels(clusters)
  
  # variance in each component explained by clusters
  R2 <- apply(x[,seq_len(n)], 2, cl=clusters, FUN=.getVE)
  
  # correlation of each component with covariates
  covar.cor <- sapply( covars, FUN=function(y) cor(x,y) )
  # correlation of the residuals (after regression on clusters) explained by each covariates
  res <- apply(x[,seq_len(n)], 2, cl=clusters, FUN=function(x, cl){
    lm(x~cl)$residuals
  })
  covar.Rcor <- sapply( covars, FUN=function(y) cor(res,y) )
  
  # each cell's distance to the cluster median
  cs <- split(row.names(x),clusters)
  dists <- lapply(n, FUN=function(i){
    x2 <- lapply(cs, FUN=function(y) x[y,1:i,drop=FALSE])
    sapply(x2, FUN=function(y){
      sqrt(colSums((t(y)-matrixStats::colMedians(y))^2))
    })
  })
  names(dists) <- paste0(n,"dims")
  
  list( silhouettes=silhouettes,
        clust.avg.silwidth=csw,
        cellDistsToMedian=dists,
        covar.cor=covar.cor,
        covar.Rcor=covar.Rcor,
        R2=R2 )
}

.aggregateDR <- function(res, dswise=FALSE){
  res <- lapply(res, FUN=function(x){
    if("dimreduction" %in% names(x)) x <- x$dimreduction
    if("evaluation" %in% names(x[[1]])) x <- lapply(x, FUN=function(x) x$evaluation)
    x
  })
  allsi <- lapply(res, FUN=function(x){
    si <- lapply(x,FUN=function(y) y$silhouettes[,3])
    pp <- parsePipNames(names(x))
    pp <- pp[rep(seq_len(nrow(pp)), sapply(si, length)),]
    pp$silhouette <- unlist(si)
    pp
  })
  allsi <- dplyr::bind_rows(allsi, .id = "dataset")
  
  perDS <- lapply(res, FUN=function(x){
    # check if the dimensions are the same
    ll <- unlist(lapply(x, FUN=function(x){ row.names(x$clust.avg.silwidth) }))
    if(!all(length(unique(table(ll))))){
      x <- lapply(x, FUN=function(x){
        row.names(x$clust.avg.silwidth)[grep("all",row.names(x$clust.avg.silwidth))] <- "all"
      })
      ll <- unlist(lapply(x, FUN=function(x){ row.names(x$clust.avg.silwidth) }))
    }
    sw <- lapply(unique(ll), FUN=function(topX){
      sapply(x, FUN=function(x) unlist(as.data.frame(x$clust.avg.silwidth)[topX,]))
    })
    names(sw) <- unique(ll)
    list( clust.avg.silwidth=sw,
          PC1.covar=sapply(x, FUN=function(x) x$covar.cor[1,]),
          PC1.covarR=sapply(x, FUN=function(x) x$covar.Rcor[1,]),
          PC.R2=sapply(x, FUN=function(x) x$R2) )
  })
  
  if(dswise) return(perDS)
  
  sw <- lapply(names(perDS[[1]]$clust.avg.silwidth), FUN=function(n){
    d <- lapply(perDS, FUN=function(x){ t(x$clust.avg.silwidth[[n]]) })
    for(f in names(perDS)) colnames(d[[f]]) <- paste(f,colnames(d[[f]]))
    do.call(cbind, d)
  })
  names(sw) <- names(perDS[[1]]$clust.avg.silwidth)
  
  PC1.covar <- lapply(row.names(perDS[[1]]$PC1.covar), FUN=function(covar){
    do.call(cbind, lapply(perDS, FUN=function(x) abs(x$PC1.covar[covar,])))
  })
  names(PC1.covar) <- paste0("PC1_covar.",row.names(perDS[[1]]$PC1.covar))
  if(is.null(perDS[[1]]$PC1.covarR)){
    PC1.covarR <- NULL
  }else{
    PC1.covarR <- lapply(row.names(perDS[[1]]$PC1.covarR), FUN=function(covar){
      do.call(cbind, lapply(perDS, FUN=function(x) abs(x$PC1.covarR[covar,])))
    })
    names(PC1.covarR) <- paste0("PC1_covarR.",row.names(perDS[[1]]$PC1.covarR))
  }
  
  PCtop5.R2 <- do.call(cbind, lapply(perDS, FUN=function(x){
    colMeans(x$PC.R2[1:5,,drop=FALSE])
  }))
  
  res <- c( silhouettes=allsi, clust.avg.silwidth=sw,
             PC1.covar, PC1.covarR )
  res$PCtop5.R2 <- PCtop5.R2
  res
}



#' match_evaluate_multiple
#' 
#' Function to match cluster labels with 'true' clusters using the Hungarian algorithm, 
#' and return precision, recall, and F1 score. Written by Lukas Weber in August 2016 as 
#' part of his \href{https://github.com/lmweber/cytometry-clustering-comparison}{cytometry
#' clustering comparison}, with just slight modifications on initial handling of 
#' input arguments.
#'
#' @param clus_algorithm cluster labels from algorithm
#' @param clus_truth true cluster labels. If NULL, will attempt to read them from the names
#' of `clus_algorithm` (expecting the format `clusterName.cellName`)
#'
#' @return A list.
#' @export
match_evaluate_multiple <- function(clus_algorithm, clus_truth=NULL) {
  library(clue)
  
  if(is.null(clus_truth)){
    clus_truth <- .getTrueLabelsFromNames(clus_algorithm)
  }
  if(length(clus_truth)>length(clus_algorithm)){
    if(!is.null(names(clus_truth)) && !is.null(names(clus_algorithm))){
      clus_truth <- clus_truth[names(clus_algorithm)]
    }else{
      stop("Could not match items from the two input vectors!")
    }
  }
  if(!is.numeric(clus_truth)) clus_truth <- as.integer(as.factor(clus_truth))
  if(!is.numeric(clus_algorithm)) clus_algorithm <- as.integer(as.factor(clus_algorithm))
  
  # number of detected clusters
  n_clus <- length(table(clus_algorithm))
  
  # remove unassigned cells (NA's in clus_truth)
  unassigned <- is.na(clus_truth)
  clus_algorithm <- clus_algorithm[!unassigned]
  clus_truth <- clus_truth[!unassigned]
  if (length(clus_algorithm) != length(clus_truth)) warning("vector lengths are not equal")
  
  tbl_algorithm <- table(clus_algorithm)
  tbl_truth <- table(clus_truth)
  
  # detected clusters in rows, true populations in columns
  pr_mat <- re_mat <- F1_mat <- matrix(NA, nrow = length(tbl_algorithm), ncol = length(tbl_truth))
  
  for (i in 1:length(tbl_algorithm)) {
    for (j in 1:length(tbl_truth)) {
      i_int <- as.integer(names(tbl_algorithm))[i]  # cluster number from algorithm
      j_int <- as.integer(names(tbl_truth))[j]  # cluster number from true labels
      
      true_positives <- sum(clus_algorithm == i_int & clus_truth == j_int, na.rm = TRUE)
      detected <- sum(clus_algorithm == i_int, na.rm = TRUE)
      truth <- sum(clus_truth == j_int, na.rm = TRUE)
      
      # calculate precision, recall, and F1 score
      precision_ij <- true_positives / detected
      recall_ij <- true_positives / truth
      F1_ij <- 2 * (precision_ij * recall_ij) / (precision_ij + recall_ij)
      
      if (F1_ij == "NaN") F1_ij <- 0
      
      pr_mat[i, j] <- precision_ij
      re_mat[i, j] <- recall_ij
      F1_mat[i, j] <- F1_ij
    }
  }
  
  # put back cluster labels (note some row names may be missing due to removal of unassigned cells)
  rownames(pr_mat) <- rownames(re_mat) <- rownames(F1_mat) <- names(tbl_algorithm)
  colnames(pr_mat) <- colnames(re_mat) <- colnames(F1_mat) <- names(tbl_truth)
  
  # match labels using Hungarian algorithm applied to matrix of F1 scores (Hungarian
  # algorithm calculates an optimal one-to-one assignment)
  
  # use transpose matrix (Hungarian algorithm assumes n_rows <= n_cols)
  F1_mat_trans <- t(F1_mat)
  
  if (nrow(F1_mat_trans) <= ncol(F1_mat_trans)) {
    # if fewer (or equal no.) true populations than detected clusters, can match all true populations
    labels_matched <- clue::solve_LSAP(F1_mat_trans, maximum = TRUE)
    # use row and column names since some labels may have been removed due to unassigned cells
    labels_matched <- as.numeric(colnames(F1_mat_trans)[as.numeric(labels_matched)])
    names(labels_matched) <- rownames(F1_mat_trans)
    
  } else {
    # if fewer detected clusters than true populations, use transpose matrix and assign
    # NAs for true populations without any matching clusters
    labels_matched_flipped <- clue::solve_LSAP(F1_mat, maximum = TRUE)
    # use row and column names since some labels may have been removed due to unassigned cells
    labels_matched_flipped <- as.numeric(rownames(F1_mat_trans)[as.numeric(labels_matched_flipped)])
    names(labels_matched_flipped) <- rownames(F1_mat)
    
    labels_matched <- rep(NA, ncol(F1_mat))
    names(labels_matched) <- rownames(F1_mat_trans)
    labels_matched[as.character(labels_matched_flipped)] <- as.numeric(names(labels_matched_flipped))
  }
  
  # precision, recall, F1 score, and number of cells for each matched cluster
  pr <- re <- F1 <- n_cells_matched <- rep(NA, ncol(F1_mat))
  names(pr) <- names(re) <- names(F1) <- names(n_cells_matched) <- names(labels_matched)
  
  for (i in 1:ncol(F1_mat)) {
    # set to 0 if no matching cluster (too few detected clusters); use character names 
    # for row and column indices in case subsampling completely removes some clusters
    pr[i] <- ifelse(is.na(labels_matched[i]), 0, pr_mat[as.character(labels_matched[i]), names(labels_matched)[i]])
    re[i] <- ifelse(is.na(labels_matched[i]), 0, re_mat[as.character(labels_matched[i]), names(labels_matched)[i]])
    F1[i] <- ifelse(is.na(labels_matched[i]), 0, F1_mat[as.character(labels_matched[i]), names(labels_matched)[i]])
    
    n_cells_matched[i] <- sum(clus_algorithm == labels_matched[i], na.rm = TRUE)
  }
  
  # means across populations
  mean_pr <- mean(pr)
  mean_re <- mean(re)
  mean_F1 <- mean(F1)
  
  return(list(n_clus = n_clus, 
              pr = pr, 
              re = re, 
              F1 = F1, 
              labels_matched = labels_matched, 
              n_cells_matched = n_cells_matched, 
              mean_pr = mean_pr, 
              mean_re = mean_re, 
              mean_F1 = mean_F1))
}


evaluateNorm <- function(x, clusters=NULL, covars=NULL){
  if(is(x,"Seurat")){
    if(is.null(covars)) covars <- x@meta.data[,c("log10_total_counts", "total_features")]
    if(is.null(clusters)) clusters <- x$phenoid
    x <- x@assays$RNA@data
  }else if(is(x, "SingleCellExperiment")){
    if(is.null(covars)) covars <- as.data.frame(colData(x)[,c("log10_total_counts", "total_features")])
    if(is.null(clusters)) clusters <- x$phenoid
    x <- logcounts(x)
  }
  x <- as.matrix(x)
  res <- data.frame( row.names=row.names(x),
            varExpByCluster=as.numeric(apply(x, 1, cl=clusters, FUN=.getVE)) ) 
  if(!is.null(covars) && length(covars)>0){
    resid <- t(apply(x,1,FUN=function(x){
      lm(x~clusters)$residuals
    }))
    for( f in names(covars) ){
      res[[f]] <- suppressWarnings(as.numeric(cor(t(x),as.numeric(covars[[f]]))))
    }
  }
  return(res)
}

