scrna_seurat_defAlternatives <- function(){
  list(
    filt="filt.lenient",
    norm=c("norm.seurat","norm.scran.noscale","norm.scran"),
    sel="sel.vst", selnb=3000,
    dr="seurat.pca", maxdim=30, clustmethod="clust.seurat",
    dims = 10, k = 20, steps = 8, 
    resolution = c(0.01, 0.1, 0.5, 0.8, 1), 
    min.size = 50 )
}

none <- function(x) x

#' filt.mad
#' 
#' Filter cells on the basis of MADs.
#'
#' @param x An object of class `SingleCellExperiment`
#' @param nmads The number of MADs above/below which to filter out (default 3)
#' @param min.cells The minimum number of cells expressing a feature (default 10) to keep
#' the feature.
#' @param min.features The minimum number of features detected in a cell (default 100) to 
#' keep the cell.
#' @param vars A named vector of control variables on which to check for deviations, in 
#' the form `variable=direction`.
#' @param outlier.times The number of times a cell should be an outlier for it to be 
#' excluded.
#'
#' @return A Seurat object.
#' @export
filt.mad <- function( x, nmads=3, min.cells=10, min.features=100,
                      vars=c( "pct_mt"="higher", 
                              "log10_total_features"="both", 
                              "log10_total_counts"="both", 
                              "pct_counts_top_50_features"="higher"
                      ),
                      outlier.times=1
){
  if(length(vars)==0) seWrap(x, min.cells=min.cells, min.features=min.features)
  out <- unlist(lapply(names(vars), o=x, v=vars, nm=nmads, FUN=function(x,o,v,nm){
    tryCatch(
      which(isOutlier(o[[x]], 
                      nmads=nm, 
                      log=F, 
                      type=v[[x]]
      )),
      error=function(e){ 
        warning(e)
        return(c())
      }
    )
  }))
  if(length(out)>0){
    out <- table(out)
    out <- as.numeric(names(out)[which(out>=outlier.times)])
    if(length(out)>0) x <- x[,-out]    
  }
  seWrap(x, min.cells=min.cells, min.features=min.features)
}

filt.lenient <- function(x){
  library(scater)
  if(!("featcount_dist" %in% colnames(colData(x)))) x <- add_meta(x)
  filters <- c( "log10_total_counts:both:5",
                "log10_total_features:both:5",
                "log10_total_features:lower:5",
                "pct_counts_in_top_20_features:both:5",
                "featcount_dist:both:5")
  out <- lapply(strsplit(filters,":"), FUN=function(f){
    which(isOutlier(x[[f[1]]], log=FALSE,
                    nmads=as.numeric(f[3]), 
                    type=f[2] ))
  })
  mtout <- isOutlier(x$pct_counts_Mt, nmads=3, type="lower" ) | 
    (isOutlier(x$pct_counts_Mt, nmads=3, type="higher" ) & x$pct_counts_Mt > 0.08)
  out <- c(out, list(mt=which(mtout)))
  out <- table(unlist(out))
  out <- as.numeric(names(out)[which(out>=2)])
  if(length(out)>0) x <- x[,-out]
  seWrap(x)
}

filt.default <- function(x, times=2){
  library(scater)
  if(!("featcount_dist" %in% colnames(colData(x)))) x <- add_meta(x)
  filters <- c( "log10_total_counts:higher:2.5",
                "log10_total_counts:lower:5",
                "log10_total_features:higher:2.5",
                "log10_total_features:lower:5",
                "pct_counts_in_top_20_features:both:5",
                "featcount_dist:both:5")
  out <- lapply(strsplit(filters,":"), FUN=function(f){
    which(isOutlier(x[[f[1]]], log=FALSE,
                    nmads=as.numeric(f[3]), 
                    type=f[2] ))
  })
  mtout <- isOutlier(x$pct_counts_Mt, nmads=3, type="lower" ) | 
    (isOutlier(x$pct_counts_Mt, nmads=2.5, type="higher" ) & x$pct_counts_Mt > 0.08)
  out <- c(out, list(mt=which(mtout)))
  out <- table(unlist(out))
  out <- as.numeric(names(out)[which(out>=times)])
  if(length(out)>0) x <- x[,-out]
  seWrap(x)
}

filt.stringent <- function(x){
  filt.default(x,1)
}

.translateFilterVars <- function(x){
  vars=c(  "mt"="pct_counts_Mt", 
           "feat"="total_features",
           "counts"="total_counts",
           "lfeat"="log10_total_features", 
           "lcounts"="log10_total_counts", 
           "top50"="pct_counts_in_top_50_features",
           "ratiodist"="featcount_dist"
  )
  x <- strsplit(x,".",fixed=T)
  y <- sapply(x,FUN=function(x){ if(length(x)==1) return("both"); x[[2]] })
  names(y) <- sapply(x,vars=vars,FUN=function(x, vars) vars[x[[1]]])
  y
}

filt.pca <- function(x, vars=NULL){
  library(scater)
  x <- runPCA(x, use_coldata=TRUE, detect_outliers=TRUE, selected_variables=vars)
  seWrap(x[,!x$outlier])
}

filt.pca2 <- function(x){
  filt.pca(x, vars=c("log10_total_counts", "log10_total_features", "pct_counts_Mt", "pct_counts_in_top_50_features")))
}

seWrap <- function(sce, min.cells=10, min.features=0){
  library(Seurat)
  se <- CreateSeuratObject( counts=counts(sce), 
                            min.cells=min.cells, 
                            min.features=min.features, 
                            meta.data=as.data.frame(colData(sce)), 
                            project = "scRNAseq" )
  se@misc$rowData <- as.data.frame(rowData(sce))
  se
}

norm.seurat <- function(x, vars=NULL, noscale=FALSE){
  x <- NormalizeData(x, verbose=FALSE)
  if(noscale){
    x <- SetAssayData(x, slot="scale.data", as.matrix(GetAssayData(x)))
  }else{
    x <- ScaleData(x, verbose=FALSE, vars.to.regress=vars)
  }
  x
}

norm.seurat.noscale <- function(x, ...){
  norm.seurat(x, noscale=TRUE, ...)
}

norm.scran <- function(x, vars=NULL, noscale=FALSE, min.mean=1){
  library(scran)
  a <- x@assays$RNA@counts
  a <- SingleCellExperiment(assays=list(counts=a))
  clusters <- quickCluster(a, min.mean=min.mean, min.size=50)
  a <- computeSumFactors(a, min.mean=min.mean, clusters=clusters)
  a <- normalize(a)
  x <- SetAssayData(x, slot="data", new.data=logcounts(a))
  if(noscale){
    x <- SetAssayData(x, slot="scale.data", as.matrix(GetAssayData(x)))
  }else{
    x <- ScaleData(x, verbose=FALSE, vars.to.regress=vars)
  }
  x
}

norm.scran.noscale <- function(x, ...){
  norm.scran(x, noscale=TRUE, ...)
}

sel.vst <- function(dat, n=2000, excl=c()){
  if(!is.null(dat@misc$vst.var.feat)){
    VariableFeatures(dat) <- dat@misc$vst.var.feat[1:min(n,length(dat@misc$vst.var.feat))]
  }else{
    dat <- FindVariableFeatures(dat, nfeatures=n)
  }
  VariableFeatures(dat) <- subsetFeatureByType(VariableFeatures(dat), excl)
  dat
}

seurat.pca <- function(x, dims=50, weight.by.var=TRUE, seed.use=42){
  RunPCA(x, features=VariableFeatures(x), verbose=FALSE, 
         weight.by.var=weight.by.var, npcs=dims, seed.use = seed.use)
}

clust.seurat <- function(x, rd=NULL, k=20, steps=8, dims=50, seed.use=1234, min.size=0, resolution=0.8){
  dims <- min(dims,ncol(x@reductions$pca@cell.embeddings))
  x <- FindNeighbors(x, k.param=k, dims=1:dims)
  x <- FindClusters(x, resolution=resolution, random.seed=seed.use, verbose=FALSE)
  Idents(x)
}

#' subsetFeatureByType
#'
#' @param g A vector of gene names, either official gene symbols or ensembl stable IDs (or
#'  `ensemblgid.symbol`). Currently only mouse or human supported.
#' @param classes A vector classes to filter.
#'
#' @return A filtered vector of gene names.
#' 
#' @export
subsetFeatureByType <- function(g, classes=c("Mt","conding","ribo")){
  if(length(classes)==0) return(g)
  classes <- match.arg(gsub("ribosomal","ribo",classes), c("Mt","coding","ribo"), several.ok=T)
  #data("ctrlgenes")
  ctrlgenes <- readRDS("/home/Shared_taupo/plger/mm_hs_control_genes_for_scRNAseq.rds")
  go <- g
  if(any(grepl("^ENSG|^ENSMUSG",head(g,n=100)))){
    ## we assume ^ENSEMBL\.whatever rownames
    cg <- lapply(classes, FUN=function(x){ union(ctrlgenes[[1]]$ensembl[[x]], ctrlgenes[[2]]$ensembl[[x]]) })
    g <- sapply(strsplit(g,".",fixed=T),FUN=function(x) x[[1]])
  }else{
    # we assume HGNC/MGI symbols
    cg <- lapply(classes, FUN=function(x){ union(ctrlgenes[[1]]$symbols[[x]], ctrlgenes[[2]]$symbols[[x]]) })
  }
  names(cg) <- classes
  if("coding" %in% classes){
    w <- which(g %in% cg$coding)
    go <- go[w]
    g <- g[w]
    cg <- cg[which(names(cg) != "coding")]
  }
  go[which(!(g %in% unlist(cg)))]
}


getDimensionality <- function(se, method, maxDims=50){
  x <- se@reductions$pca@cell.embeddings
  switch(method,
         essLocal.a=essLocalDimEst(x),
         essLocal.b=essLocalDimEst(x, ver="b"),
         pcaLocal.FO=pcaLocalDimEst(x,ver="FO"),
         pcaLocal.fan=pcaLocalDimEst(x, ver="fan"),
         pcaLocal.maxgap=pcaLocalDimEst(x, ver="maxgap"),
         maxLikGlobal=maxLikGlobalDimEst(x, k=20, unbiased=TRUE),
         pcaOtpmPointwise.max=pcaOtpmPointwiseDimEst(x,N=10),
         elbow=farthestPoint(se@reductions$pca@stdev)-1,
         fisherSeparability=FisherSeparability(x),
         scran.denoisePCA=scran.ndims.wrapper(se),
         jackstraw.elbow=js.wrapper(se,n.dims=ncol(x)-1,ret="ndims")
  )
}

#' applyFilterString
#'
#' @param sce A SingleCellExperiment object.
#' @param filterstring A filtering string.
#'
#' @return A Seurat object.
#' @export
applyFilterString <- function(sce, filterstring){
  x <- strsplit(filterstring,"_",fixed=T)[[1]]
  mads <- as.numeric(x[[2]])
  vars <- .translateFilterVars(strsplit(x,";",fixed=T)[[1]])
  otimes <- ifelse(is.null(x[[3]]),1,as.numeric(x[[3]]))
  filt.mad(sce, nmads=mads, vars=vars, outlier.times=otimes)
}

.tmads <- function(x, nbmads=2.5){
  x2 <- nbmads*median(abs(x-median(x)))
  median(x)+c(-x2,x2)
}

filt.stringent <- function(x){
  filt.default(x,1)
}

.translateFilterVars <- function(x){
  vars=c(  "mt"="pct_counts_Mt", 
           "feat"="total_features",
           "counts"="total_counts",
           "lfeat"="log10_total_features", 
           "lcounts"="log10_total_counts", 
           "top50"="pct_counts_in_top_50_features",
           "ratiodist"="featcount_dist"
  )
  x <- strsplit(x,".",fixed=T)
  y <- sapply(x,FUN=function(x){ if(length(x)==1) return("both"); x[[2]] })
  names(y) <- sapply(x,vars=vars,FUN=function(x, vars) vars[x[[1]]])
  y
}

#' getFilterStrings
#'
#' Returns a combination of filtering strings.
#' 
#' @param mads A vector of number of MADs.
#' @param times A vector of outlier times.
#'
#' @return A vector of filtering strings.
#' @export
getFilterStrings <- function(mads=c(2,2.5,3,5), times=1:3){
  varCombs <- c("","mt","lcounts","mt;lcounts","lfeat;lcounts","mt;lfeat","mt;lfeat;lcounts",
                "mt;lfeat;kcounts;top50","mt;top50","lcounts;top50")
  v2 <- varCombs[grep("lfeat|lcounts", varCombs)]
  v2 <- gsub("lfeat","feat", v2)
  v2 <- gsub("lcounts","counts", v2)
  varCombs <- c(varCombs,v2)
  dirs <- c("lower","higher","both")
  mads<- c(2, 2.5, 3, 5)
  varCombs <- unlist(lapply(strsplit(varCombs,";",fixed=T), dir=dirs, mads=mads, FUN=function(x, dir, mads){
    y <- expand.grid(lapply(x, y=dir, sep=".", FUN=paste))
    apply(y,1,collapse=";",FUN=paste)
  }))
  eg <- expand.grid(varCombs, mads, times)
  nbVars <- sapply(strsplit(as.character(eg[,1]),";"),FUN=length)
  eg <- eg[which(nbVars>=eg[,3]),]
  as.character(apply(eg, 1, collapse="_", FUN=paste))
}

doublet.scDblFinder <- function(x){
  library(scDblFinder)
  x <- scDblFinder(x, verbose=FALSE)
  x[,which(x$scDblFinder.class!="doublet")]
}