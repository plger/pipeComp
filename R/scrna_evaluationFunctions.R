#' evaluateClustering
#'
#' Evaluates a clustering using 'true' labels. Entries with missing true labels
#' (i.e. NA) are excluded from calculations. If using `evaluteClustering` in a
#' custom pipeline, you might want to use the corresponding 
#' `pipeComp:::.aggregateClusterEvaluation` aggregation function.
#'
#' @param x The clustering labels
#' @param tl The true labels
#'
#' @return A numeric vector of metrics (see the `pipeComp_scRNA` vignette for 
#' details)
#' @importFrom aricode clustComp
#' @export
#' 
#' @examples
#' # random data
#' dat <- data.frame( 
#'  cluster=rep(LETTERS[1:3], each=10),
#'  x=c(rnorm(20, 0), rnorm(10, 1)),
#'  y=c(rnorm(10, 1), rnorm(20, 0))
#' )
#' # clustering
#' dat$predicted <- kmeans(dist(dat[,-1]),3)$cluster
#' # evaluation
#' evaluateClustering(dat$predicted, dat$cluster)
evaluateClustering <- function(x, tl=NULL){
  if(is.null(tl)) tl <- attr(x,"true.labels")
  if(is.null(tl)) stop("True labels not found!")
  if(!is.null(names(tl)) && !is.null(names(x))){
    if(all(names(x) %in% names(tl))){
      tl <- tl[names(x)]
    }else{
      warning("The names of the clustering vector does not match the names of",
              " the true labels.")
    }
  }
  if(any(is.na(tl))){
    w <- which(!is.na(tl))
    if(length(w)==0) stop("No non-NA true label!")
    tl <- tl[w]
    x <- x[w]
  }
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

#' @importFrom dplyr bind_rows
.aggregateExcludedCells <- function(res){
  pi <- parsePipNames(names(res[[1]]))
  res <- lapply(res, FUN=function(x){
    y <- pi[rep(seq_len(nrow(pi)), each=ncol(x[[1]]$table)),,drop=FALSE]
    row.names(y) <- NULL
    y$subpopulation <- rep(colnames(x[[1]]$table), length(x))
    y$N.before <- unlist(lapply(x, FUN=function(x) as.numeric(x$table[1,])))
    y$N.lost <- unlist(lapply(x, FUN=function(x){
      nl <- as.numeric(x$table[1,]-x$table[2,])
      if(any(is.na(nl))){
        w <- which(is.na(nl))
        nl[w] <- as.numeric(x$table[1,w])
      }
      nl
    }))
    
    if(length(w <- which(is.na(y$N.lost)))>0) y$N.lost[w] <- y$N.before[w]
    y$pc.lost <- round(100*y$N.lost/y$N.before,3)
    y
  })
  res <- dplyr::bind_rows(res, .id="dataset")
  res$dataset <- factor(res$dataset)
  res
}


#' @importFrom matrixStats rowMins
#' @importFrom dplyr bind_cols bind_rows
.aggregateClusterEvaluation <- function(res){
  pi <- row.names(res[[1]])
  if(is.null(pi)) pi <- names(res[[1]])
  pi <- parsePipNames(pi)
  res <- lapply(res, FUN=function(a){
    a <- lapply(a, cn=unique(unlist(lapply(a, names))), FUN=function(x,cn){
      for(f in setdiff(cn,names(x))) x[[f]] <- NA_real_
      x[cn]
    })
    ## accomodate old and new dplyr behaviors:
    x <- dplyr::bind_rows(a)
    if(ncol(x)!=length(a[[1]])) x <- t(x)
    ##
    x <- as.data.frame(x)
    colnames(x) <- names(a[[1]])
    for(f in c("pr","re","F1")){
      w <- grep(paste0("^",f,"\\."),colnames(x))
      x[[paste0("min_",f)]]<- matrixStats::rowMins(as.matrix(x[,w,drop=FALSE]))
    }
    y <- as.data.frame(x[,grep("pr\\.|re\\.|F1\\.",colnames(x),invert=TRUE)])
    y$true.nbClusts <- length(grep("^pr\\.", colnames(x)))
    y
  })
  pi <- pi[rep(seq_len(nrow(pi)), length(res)),,drop=FALSE]
  res <- cbind(pi, dplyr::bind_rows(res, .id="dataset"))
  res$dataset <- factor(res$dataset)
  row.names(res) <- NULL
  res
}

.getVE <- function(x, cl){
  tryCatch({
    summary(lm(x~cl))$r.squared
  }, error=function(e){ print(e); return(NA) })
}



#' evaluateDimRed
#'
#' Gathers evaluation statistics on a reduced space using known cell labels.
#' If using `evaluteDimRed` in a custom pipeline, you will probably want to use
#' `pipeComp:::.aggregateDR` as the corresponding aggregation function.
#'
#' @param x The matrix of the reduced space, with cells as rows and components 
#' as columns
#' @param clusters The vector indicating each cell's cluster.
#' @param n A numeric vector indiciating the number of top dimensions at which 
#' to gather statistics (default `c(10,20,50)`). Will use all available 
#' dimensions if a higher number is given.
#' @param covars A character vectors containing any additional covariates 
#' (column names of `colData`) to track during evalutation. If missing, will
#' attempt to use default covariates. To disable, set `covars=c()`.
#'
#' @return A list with the following components:
#' * silhouettes: a matrix of the silhouette for each cell-cluster pair at each 
#' value of `n`
#' * clust.avg.silwidth: a matrix of the cluster average width at each value of 
#' `n`
#' * R2: the proportion of variance in each component (up to `max(n)`) that is 
#' explained by the clusters (i.e. R-squared of a linear model).
#' 
#' @importFrom stats dist lm
#' @importFrom cluster silhouette
#' @export
#' @examples
#' # random data
#' library(scater)
#' sce <- runPCA(logNormCounts(mockSCE(ngenes = 500)))
#' sce <- addPerCellQC(sce)
#' # random population labels
#' sce$cluster <- sample(LETTERS[1:3], ncol(sce), replace=TRUE)
#' res <- evaluateDimRed(sce, sce$cluster, covars=c("sum","detected"))
#' # average silhouette widths:
#' res$clust.avg.silwidth
#' # adjusted R2 of covariates:
#' res$covar.adjR2
evaluateDimRed <- function(x, clusters=NULL, n=c(10,20,50), covars){
  if(missing(covars)) covars <- c("log10_total_features", "log10_total_counts", 
                                  "total_features")
  if(is(x,"Seurat")){
    if(is.character(covars)) covars <- x[[]][,covars]
    if(is.null(clusters)) clusters <- x$phenoid
    x <- Embeddings(x[["pca"]])
  }else if(is(x, "SingleCellExperiment")){
    if(is.character(covars))
      covars <- as.data.frame(colData(x)[,
                                    intersect(covars,colnames(colData(x)))])
    if(is.null(clusters)) clusters <- x$phenoid
    x <- reducedDim(x, "PCA")
  }else{
    if(is.null(clusters)) stop("`clusters` must be given!")
    if(is.character(covars)) covars <- list()
  }
  if(is.data.frame(covars) && ncol(covars)==0) covars <- list()
  clusters <- as.factor(clusters)
  n <- unique(vapply(n, y=ncol(x), FUN=function(x,y) min(x,y), numeric(1)))
  si <- lapply(n, FUN=function(dims){
    cluster::silhouette(as.integer(clusters), dist(x[,seq_len(dims)]))
  })
  names(si) <- vapply(n, FUN=function(i){
    if(i==ncol(x)) return(paste0("all_", i,"_dims"))
    return(paste0("top_", i,"_dims"))
  }, FUN.VALUE=character(1))
  # summarize silhouette information
  if(length(n)==1){
    silhouettes <- si[[1]][,seq(1,3)]
  }else{
    silhouettes <- lapply(si,FUN=function(x) as.numeric(x[,3]))
    silhouettes <- cbind(si[[1]][,seq(1,2)], do.call(cbind, silhouettes))
  }
  # cluster average silhouette width
  csw <- do.call(rbind, lapply(si, FUN=function(x) summary(x)$clus.avg.widths))
  colnames(csw) <- levels(clusters)
  # variance in each component explained by clusters
  R2 <- apply(x[,seq_len(max(n)),drop=FALSE], 2, cl=clusters, FUN=.getVE)
  # correlation of each component with covariates
  covar.cor <- vapply( covars, FUN=function(y) cor(x,y), numeric(ncol(x)) )
  covar.adjR2 <- vapply(covars, FUN=function(co){
    apply(x[,seq_len(min(5,ncol(x)))], 2, FUN=function(x){ 
      tryCatch({
        summary(lm(x~co+factor(clusters)))$adj.r.squared -
          summary(lm(x~factor(clusters)))$adj.r.squared
      }, error=function(e) NA)
    })
  }, numeric(min(5,ncol(x))))
  # per-subpopuluation correlation of each component with covariates
  ii <- split(seq_len(nrow(x)), clusters)
  covar.cor2 <- lapply( covars, FUN=function(y){
    z <- vapply(ii, FUN=function(i){
      t(cor(x[i,],y[i]))
    }, numeric(ncol(x)))
    row.names(z) <- colnames(x)
    z
  })
  # each cell's distance to the cluster median
  cs <- split(row.names(x),clusters)
  dists <- lapply(n, FUN=function(i){
    x2 <- lapply(cs, FUN=function(y) x[y,seq_len(i),drop=FALSE])
    lapply(x2, FUN=function(y){
      sqrt(colSums((t(y)-matrixStats::colMedians(y))^2))
    })
  })
  names(dists) <- paste0(n,"dims")
  list( silhouettes=silhouettes,
        clust.avg.silwidth=csw,
        cellDistsToMedian=dists,
        covar.cor=covar.cor,
        covar.cor2=covar.cor2,
        covar.adjR2=covar.adjR2,
        R2=R2 )
}

#' @importFrom stats as.formula median
.aggregateDR <- function(res){
  res <- lapply(res, FUN=function(x){
    if("evaluation" %in% names(x)) x <- x$evaluation
    if("dimreduction" %in% names(x)) x <- x$dimreduction
    if("evaluation" %in% names(x)) x <- x$evaluation
    x
  })
  
  pi <- parsePipNames(names(res[[1]]))
  
  agfns <- list(minSilWidth=min, meanSilWidth=mean, medianSilWidth=median, 
                maxSilWidth=max)
  
  allsi <- lapply(res, FUN=function(x){
    subpops <- colnames(x[[1]]$clust.avg.silwidth)
    x <- lapply(x,FUN=function(y) y$silhouettes)
    dims <- table(unlist(lapply(x, FUN=function(x) colnames(x)[c(-1,-2)])))
    if(!any(dims==length(x))){
      warning("Silhouettes computed over incompatible dimensionalities.",
              "Will use the first available one.")
      dims <- c(1,1)
    }
    if(length(unique(dims))==1){
      # single dimensionality
      if(length(dims)>1){
        x <- lapply(x, FUN=function(x){ colnames(x)[3] <- "selected"; x })
        dims <- table(unlist(lapply(x, FUN=function(x) colnames(x)[c(-1,-2)])))
      }
    }
    names(dims) <- dims <- intersect( unique(unlist(lapply(x, FUN=colnames))),
                                      names(dims)[dims==length(x)] )
    lapply(dims, FUN=function(dim){
      sils <- dplyr::bind_rows(lapply(x, FUN=function(sil){ 
        data.frame( subpopulation=subpops,
                    vapply(agfns, FUN.VALUE=numeric(length(unique(sil[,1]))),
                           FUN=function(agf) aggregate(sil[,dim], 
                                                       by=list(cluster=sil[,1]),
                                                       FUN=agf)[,2] ),
                    stringsAsFactors=FALSE )
      }))
      sils<- cbind(pi[rep(seq_len(nrow(pi)),each=length(subpops)),,drop=FALSE],
                    sils)
      row.names(sils) <- NULL
      sils
    })
  })
  
  dims <- unlist(lapply(allsi, names))
  dims <- table(dims)[unique(dims)]
  names(dims) <- dims <- names(dims)[dims==length(allsi)]
  allsi <- lapply(dims, FUN=function(x){
    x <- dplyr::bind_rows(lapply(allsi, FUN=function(y) 
                                              y[[x]]), .id = "dataset")
    x$dataset <- factor(x$dataset)
    x
  })

  R2 <- dplyr::bind_rows(lapply(res, FUN=function(x){
    # find common PCs
    x <- lapply(x, FUN=function(x) x$R2)
    nn <- table(unlist(lapply(x, names)))
    nn <- names(nn)[nn==length(x)]
    x <- do.call(rbind, lapply(x, FUN=function(x) x[nn]))
    cbind(pi, x)
  }), .id="dataset")
  
  if(is.null(tmp <- dim(res[[1]][[1]]$covar.cor)) || any(tmp==0)){
    top5 <- covar.adjR2 <- covar <- NULL
  }else{
    covar <- dplyr::bind_rows(lapply(res, FUN=function(x){
      a <- pi[rep(seq_len(nrow(pi)),each=ncol(x[[1]]$covar.cor)),,drop=FALSE]
      a <- cbind(a, dplyr::bind_rows(lapply(x, FUN=function(x){
        y <- t(x$covar.cor)
        colnames(y) <- paste0("PC",seq_len(ncol(y)))
        data.frame(covariate=colnames(x$covar.cor), y, stringsAsFactors=FALSE)
     })))
      row.names(a) <- NULL
      a
    }), .id="dataset")
    
    top5 <- covar.cor2 <- covar.adjR2 <- NULL
    if(!is.null(res[[1]][[1]]$covar.cor2)){
      covar.cor2 <- dplyr::bind_rows(lapply(res, FUN=function(x){
        x <- lapply(x, FUN=function(x){
          i <- seq_len(min(nrow(x),5))
          reshape2::melt(vapply(x$covar.cor2, FUN.VALUE=numeric(length(i)), 
                                FUN=function(x) rowMeans(x[i,,drop=FALSE])),
                         value.name = "meanCor")
        })
        cbind( pi[rep(seq_len(nrow(pi)),vapply(x,nrow,integer(1))),,drop=FALSE], 
               dplyr::bind_rows(x) )
      }), .id="dataset")
      colnames(covar.cor2)[ncol(pi)+2:3] <- c("component","covariate")
      w <- which( covar.cor2$component %in% 
                    paste0(c("PC","PC_"), rep(seq(1,5),each=2)) )
      top5 <- aggregate( covar.cor2$meanCor[w],
                         by=covar.cor2[,c("dataset","covariate", colnames(pi))],
                         FUN=function(x) mean(abs(x)) )
      ff <- paste( paste( setdiff(colnames(top5), c("x","covariate")),
                          collapse="+"),
                   "~covariate" )
      top5 <- reshape2::dcast( top5, as.formula(ff), value.var="x", 
                               fun.aggregate=mean)
    }
    if(!is.null(res[[1]][[1]]$covar.adjR2)){
      covar.adjR2 <- dplyr::bind_rows(lapply(res, FUN=function(x){
        pi <- parsePipNames(names(x))
        x <- do.call(rbind, 
                     lapply(x, FUN=function(x) x$covar.adjR2[1,,drop=FALSE]))
        row.names(x) <- NULL
        cbind(pi, x)
      }), .id="dataset")
    }
  }
  list( silhouette=allsi, varExpl.subpops=R2, corr.covariate=covar,
        meanAbsCorr.covariate2=top5, PC1.covar.adjR2=covar.adjR2 )
}


#' match_evaluate_multiple
#' 
#' Function to match cluster labels with 'true' clusters using the Hungarian 
#' algorithm, and return precision, recall, and F1 score. Written by Lukas 
#' Weber in August 2016 as part of his 
#' \href{https://github.com/lmweber/cytometry-clustering-comparison}{cytometry
#' clustering comparison}, with just slight modifications on initial handling 
#' of input arguments.
#'
#' @param clus_algorithm cluster labels from algorithm
#' @param clus_truth true cluster labels. If NULL, will attempt to read them 
#' from the names of `clus_algorithm` (expecting the format 
#' `clusterName.cellName`)
#'
#' @return A list.
#' 
#' @importFrom clue solve_LSAP
#' @export
#' 
#' @examples
#' # random data
#' dat <- data.frame( 
#'  cluster=rep(LETTERS[1:3], each=10),
#'  x=c(rnorm(20, 0), rnorm(10, 1)),
#'  y=c(rnorm(10, 1), rnorm(20, 0))
#' )
#' # clustering
#' dat$predicted <- kmeans(dist(dat[,-1]),3)$cluster
#' # evaluation
#' match_evaluate_multiple(dat$predicted, dat$cluster)
match_evaluate_multiple <- function(clus_algorithm, clus_truth=NULL){
  if(is.null(clus_truth)){
    warning("Truth not given, trying to guess it from the names...")
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
  if(!is.numeric(clus_algorithm)) 
    clus_algorithm <- as.integer(as.factor(clus_algorithm))
  
  # number of detected clusters
  n_clus <- length(table(clus_algorithm))
  
  # remove unassigned cells (NA's in clus_truth)
  unassigned <- is.na(clus_truth)
  clus_algorithm <- clus_algorithm[!unassigned]
  clus_truth <- clus_truth[!unassigned]
  if (length(clus_algorithm) != length(clus_truth)) 
    warning("vector lengths are not equal")
  
  tbl_algorithm <- table(clus_algorithm)
  tbl_truth <- table(clus_truth)
  
  # detected clusters in rows, true populations in columns
  pr_mat <- re_mat <- F1_mat <- 
    matrix(NA, nrow = length(tbl_algorithm), ncol = length(tbl_truth))
  
  for (i in seq_len(length(tbl_algorithm))) {
    for (j in seq_len(length(tbl_truth))) {
      i_int <- as.integer(names(tbl_algorithm))[i]  # cluster from algorithm
      j_int <- as.integer(names(tbl_truth))[j]  # cluster from true labels
      
      true_positives <- sum(clus_algorithm == i_int & 
                              clus_truth == j_int, na.rm = TRUE)
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
  
  # put back cluster labels (note some row names may be missing due to removal 
  # of unassigned cells)
  rownames(pr_mat) <- rownames(re_mat) <- rownames(F1_mat) <- 
    names(tbl_algorithm)
  colnames(pr_mat) <- colnames(re_mat) <- colnames(F1_mat) <- 
    names(tbl_truth)
  
  # match labels using Hungarian algorithm applied to matrix of F1 scores 
  # (Hungarian algorithm calculates an optimal one-to-one assignment)
  
  # use transpose matrix (Hungarian algorithm assumes n_rows <= n_cols)
  F1_mat_trans <- t(F1_mat)
  
  if (nrow(F1_mat_trans) <= ncol(F1_mat_trans)) {
    # if fewer (or equal no.) true populations than detected clusters, can 
    # match all true populations
    labels_matched <- solve_LSAP(F1_mat_trans, maximum = TRUE)
    # use row and column names since some labels may have been removed due to 
    # unassigned cells
    labels_matched <- 
      as.numeric(colnames(F1_mat_trans)[as.numeric(labels_matched)])
    names(labels_matched) <- rownames(F1_mat_trans)
    
  } else {
    # if fewer detected clusters than true populations, use transpose matrix
    # and assign NAs for true populations without any matching clusters
    labels_matched_flipped <- solve_LSAP(F1_mat, maximum = TRUE)
    # use row and column names since some labels may have been removed due to 
    # unassigned cells
    labels_matched_flipped <- 
      as.numeric(rownames(F1_mat_trans)[as.numeric(labels_matched_flipped)])
    names(labels_matched_flipped) <- rownames(F1_mat)
    
    labels_matched <- rep(NA, ncol(F1_mat))
    names(labels_matched) <- rownames(F1_mat_trans)
    labels_matched[as.character(labels_matched_flipped)] <- 
      as.numeric(names(labels_matched_flipped))
  }
  
  # precision, recall, F1 score, and number of cells for each matched cluster
  pr <- re <- F1 <- n_cells_matched <- rep(NA, ncol(F1_mat))
  names(pr) <- names(re) <- names(F1) <- names(n_cells_matched) <- 
    names(labels_matched)
  
  for (i in seq_len(ncol(F1_mat))) {
    # set to 0 if no matching cluster (too few detected clusters); use 
    # character names for row and column indices in case subsampling completely 
    #  removes some clusters
    im <- labels_matched[i]
    pr[i] <- ifelse(is.na(im), 0, pr_mat[as.character(im), names(im)])
    re[i] <- ifelse(is.na(im), 0, re_mat[as.character(im), names(im)])
    F1[i] <- ifelse(is.na(im), 0, F1_mat[as.character(im), names(im)])
    
    n_cells_matched[i] <- sum(clus_algorithm == labels_matched[i], na.rm=TRUE)
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

#' evaluateNorm
#' 
#' @param x An object of class 'Seurat' or 'SingleCellExperiment' with 
#' normalized data
#' @param clusters A vector of true cluster identities. If missing, will 
#' attempt to fetch it from the `phenoid` colData.
#' @param covars Covariates to include, either as data.frame or as colData
#' columns of `x`
#'
#' @return a data.frame.
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom matrixStats rowMedians
#' @importFrom Matrix rowMeans
#' @importFrom Seurat GetAssayData
#' @examples
#' # random data
#' library(scater)
#' sce <- logNormCounts(mockSCE(ngenes = 500))
#' sce <- addPerCellQC(sce)
#' # random population labels
#' sce$cluster <- sample(LETTERS[1:3], ncol(sce), replace=TRUE)
#' evaluateNorm(sce, sce$cluster, covars="detected")
evaluateNorm <- function(x, clusters=NULL, covars){
  if(missing(covars)) covars <- c("log10_total_counts", "total_features")
  if(is(x,"Seurat")){
    if(is.character(covars) && length(covars)>0) 
      covars <- x[[]][,covars,drop=FALSE]
    if(is.null(clusters)) clusters <- x$phenoid
    meanCount <- Matrix::rowMeans(GetAssayData(x, assay="RNA", slot="counts"))
    x <- GetAssayData(x, assay="RNA", slot="data")
  }else if(is(x, "SingleCellExperiment")){
    if(is.character(covars) && length(covars)>0) 
      covars <- as.data.frame(colData(x)[,covars,drop=FALSE])
    if(is.null(clusters)) clusters <- x$phenoid
    meanCount <- Matrix::rowMeans(counts(x))
    x <- logcounts(x)
  }
  x <- as.matrix(x)
  ve <- round(as.numeric(apply(x, 1, cl=clusters, FUN=.getVE)),2)
  res <- data.frame( row.names=row.names(x), meanCount=meanCount,
                     varExpByCluster=ve )
  if(!is.null(covars) && length(covars)>0){
    ii <- split(seq_len(ncol(x)), clusters)
    for( f in names(covars) ){
      tmp <- vapply(ii, FUN.VALUE=numeric(nrow(x)), FUN=function(i){
        cor(t(x[,i]), as.numeric(covars[[f]][i]))
      })
      res[[paste0(f,".medianCorr")]] <- round(matrixStats::rowMedians(tmp),2)
      res[[paste0(f,".meanCorr")]] <- round(Matrix::rowMeans(tmp),2)
      res[[paste0(f,".meanAbsCorr")]] <- round(Matrix::rowMeans(abs(tmp)),2)
    } 
  }
  return(res)
}


.aggregateNorm <- function(res){
  res <- lapply(res, FUN=function(x){
    if("evaluation" %in% names(x)) x <- x$evaluation
    if("normalization" %in% names(x)) x <- x$normalization
    x
  })
  res2 <- dplyr::bind_rows(lapply(res, FUN=function(x){
    x <- t(vapply(x, FUN.VALUE=numeric(ncol(x[[1]])), FUN=function(x){
      vapply(x, FUN=function(x) mean(abs(x), na.rm=TRUE), numeric(1))
    }))
    parsePipNames(as.data.frame(x))
  }), .id="dataset")
  cor2 <- dplyr::bind_rows(lapply(res, FUN=function(x){
    parsePipNames(as.data.frame(t(vapply(x, FUN=function(x){
      vapply( x[,-1], y=x[,1], method="spearman", use="pairwise",
              FUN=cor, FUN.VALUE=numeric(1) )
    }, FUN.VALUE=numeric(ncol(x[[1]])-1)))))
  }), .id="dataset")
  list( meanAbsCorr=res2, spearman.with.meanCount=cor2 )
}

