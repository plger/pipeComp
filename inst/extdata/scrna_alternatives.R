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
