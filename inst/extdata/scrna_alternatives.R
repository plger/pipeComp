scrna_seurat_defAlternatives <- function(x=list()){
  def <- list(
    filt="filt.lenient",
    norm=c("norm.seurat","norm.scran.noscale","norm.scran"),
    sel="sel.vst", selnb=2000,
    dr="seurat.pca", maxdim=30, clustmethod="clust.seurat",
    dims = 10, k = 20, steps = 8, 
    resolution = c(0.01, 0.1, 0.5, 0.8, 1), 
    min.size = 50 )
  for(f in names(x)) def[[f]] <- x[[f]]
  def
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
  filt.pca(x, vars=c("log10_total_counts", "log10_total_features", "pct_counts_Mt", "pct_counts_in_top_50_features"))
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
    if(is.null(vars)) vars <- c()
    for(f in vars){
      if(!(sd(x[[]][[f]])>0)) vars <- setdiff(vars,f)
    }
    if(length(vars)==0) vars <- NULL
    x <- ScaleData(x, verbose=FALSE, vars.to.regress=vars)
  }
  x
}

norm.scran <- function(x, vars=NULL, noscale=TRUE, min.mean=1){
  library(scran)
  if(is(x,"Seurat")){
    a <- GetAssayData(x, assay = "RNA", slot = "counts")
  }else{
    a <- counts(x)
  }
  a <- SingleCellExperiment(assays=list(counts=a))
  clusters <- quickCluster(a, min.mean=min.mean, min.size=50)
  a <- computeSumFactors(a, min.mean=min.mean, clusters=clusters)
  a <- normalize(a)
  if(is(x,"Seurat")){
    x <- SetAssayData(x, slot="data", new.data=logcounts(a))
    if(noscale){
      x <- SetAssayData(x, slot="scale.data", as.matrix(GetAssayData(x)))
    }else{
      x <- ScaleData(x, verbose=FALSE, vars.to.regress=vars)
    }
  }else{
    if(!noscale) a <- t(scale(t(a)))
    logcounts(x) <- a
  }
  x
}
norm.scran.scaled <- function(x, ...){
  norm.scran(x, noscale=FALSE, ...)
}

norm.none <- function(x, vars=NULL, noscale=TRUE){
  x <- SetAssayData(x, slot="data", log1p(GetAssayData(x, slot="counts")))
  if(noscale){
    x <- SetAssayData(x, slot="scale.data", as.matrix(GetAssayData(x)))
  }else{
    x <- ScaleData(x, verbose=FALSE, vars.to.regress=vars)
  }
  x
}
norm.none.scaled <- function(x){
  norm.none(x, noscale=FALSE)
}

#' norm.seuratvst
#'
#' A wrapper around `sctransform` variance stabilizing transformation.
#' 
#' @param x A Seurat object.
#' @param vars A vector of variables to regress when scaling (default none)
#' @param noscale Ignored.
#' @param variable.features.n Passed to `SCTransform`, default 5000
#'
#' @return A Seurat object with updated data slot.
norm.seuratvst <- function(x, vars=NULL, noscale=FALSE, variable.features.n=5000){
  library(sctransform)
  x <- SCTransform(x, vars.to.regress=vars, verbose=FALSE, 
                   variable.features.n=variable.features.n, 
                   return.only.var.genes=FALSE)
  x@misc$vst.var.feat <- VariableFeatures(x)
  x
}
norm.sctransform <- norm.seuratvst


#' norm.scnorm
#'
#' A wrapper around `SCnorm` normalization.
#' 
#' @param x A SCE or Seurat object.
#' @param vars A vector of variables to regress when scaling (default none). Ignored if `noscale`.
#' @param noscale Logical; whether to disable scaling (default FALSE)
#'
#' @return An object of the same class as `x` with updated slots.
norm.scnorm <- function(x, vars=NULL, noscale=TRUE, nthreads=1){
  library(SCnorm)
  if(is(x,"Seurat")){
    a <- Seurat::GetAssayData(x, slot="counts")
    a <- SCnorm(a, rep("A",ncol(a)), NCores=nthreads)
    a <- log1p(assays(a)$normcounts)
    x <- SetAssayData(x, slot="data", new.data=a)
    if(noscale){
      x <- SetAssayData(x, slot="scale.data", as.matrix(GetAssayData(x)))
    }else{
      x <- ScaleData(x, verbose=FALSE, vars.to.regress=vars)
    }
    return(x)
  }
  a <- counts(x)
  a <- SCnorm(a, rep("A",ncol(a)), NCores=nthreads)
  a <- log1p(assays(a)$normcounts)
  if(!noscale) a <- t(scale(t(a)))
  logcounts(x) <- a
  x
}
norm.scnorm.scaled <- function(x, ...){
  norm.scnorm(x, noscale=FALSE, ...)
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

#' applySelString
#' 
#' Applies a given selection string.
#'
#' @param se A Seurat object.
#' @param selstring A rowData variable to use, or a selection string, e.g. 'vst:2000:coding_rmMt_rmribo'.
#' @param n The number of genes to select (ignored if selstring is a full selection string).
#'
#' @return A filtered Seurat object.
#' @export
applySelString <- function(se, selstring, n=2000){
  #vst:2000:coding_rmMt_rmribo
  x <- strsplit(selstring,":",fixed=T)[[1]]
  excl <- c()
  if(length(x)>2 && x[3]!="") excl <- gsub("rm","",strsplit(x[3], "_")[[1]])
  fn <- paste0("sel.",x[1])
  if(exists(x[2])) n <- as.numeric(x[2])
  if(exists(fn) && is.function(get(fn))) return(get(fn)(se, n=n, excl=excl))
  sel.fromField(se, x[1], n, excl)
}

#' sel.fromField
#' 
#' Selection of features based on a given rowData field (in decreasing order).
#'
#' @param dat A `Seurat` object.
#' @param f The field to use.
#' @param n The number of features to select (default 2000)
#' @param excl Feature types to exclude (default none)
#'
#' @return A `Seurat` object with updated `VariableFeatures`
#' @export
sel.fromField <- function( dat, f, n=2000, excl=c() ){
  if(is.null(dat@misc$rowData[[f]])) return(NULL)
  e <- dat@misc$rowData[row.names(dat),f]
  VariableFeatures(dat) <- row.names(dat)[order(e, decreasing=T)[1:min(n,length(e))]]
  VariableFeatures(dat) <- subsetFeatureByType(VariableFeatures(dat), excl)
  dat
}

sel.deviance <- function(x, n=2000, excl=c()){
  sel.fromField(x, "deviance", n=n, excl=excl)
}
sel.expr <- function(x, n=2000, excl=c()){
  sel.fromField(x, "total_counts", n=n, excl=excl)
}

seurat.pca <- function(x, dims=50, weight.by.var=TRUE, seed.use=42){
  RunPCA(x, features=VariableFeatures(x), verbose=FALSE, 
         weight.by.var=weight.by.var, npcs=dims, seed.use = seed.use)
}
seurat.pca.noweight <- function(x, dims=50, weight.by.var=FALSE, seed.use=42){
  RunPCA(x, features=VariableFeatures(x), verbose=FALSE, 
         weight.by.var=weight.by.var, npcs=dims, seed.use = seed.use)
}

scran.denoisePCA <- function(x, dims=50, pca.method=c("exact","irlba"), ...){
  library(scran)
  BSPARAM <- switch(match.arg(pca.method),
                    exact=ExactParam(),
                    irlba=IrlbaParam() )
  if(is(x, "Seurat")){
    sce <- SingleCellExperiment( list( counts=GetAssayData(x, slot="counts"),
                                       logcounts=GetAssayData(x, slot="data")) )
    return(sceDR2seurat(reducedDim(sce, "PCA"), x, "pca"))
  }
  if(packageVersion("scran") >= "1.13"){
    var.stats <- modelGeneVar(sce)
    sce <- denoisePCA(sce, technical=var.stats, min.rank=2, max.rank=dims, BSPARAM=BSPARAM, ...)
  }else{
    td <- trendVar(sce, use.spikes=FALSE)
    sce <- denoisePCA(sce, technical=td$trend, min.rank=2, max.rank=dims, BSPARAM=BSPARAM, ...)
  }
  if(is(x, "Seurat")) return(sceDR2seurat(reducedDim(sce, "PCA"), x, "pca"))
  return(sce)
}

seGlmPCA <- function(x, weight.by.var=TRUE, dims=20){
  library(glmpca)
  dr <- glmpca(as.matrix(GetAssayData(x, assay = "RNA", slot = "counts")[VariableFeatures(x),]), dims)
  e <- as.matrix(dr$factors)
  colnames(e) <- gsub("dim","dim_",colnames(e))
  if(weight.by.var=="both"){
    x[["glmpca"]] <- CreateDimReducObject(embeddings=e, key="dim_", assay="RNA")
    e <- t(t(e)*dr$d)
    x[["glmpca.wt"]] <- CreateDimReducObject(embeddings=e, key="dim_", assay="RNA")
  }else{
    if(weight.by.var) e <- t(t(e)*dr$d)
    x[["pca"]] <- CreateDimReducObject(embeddings=e, key="dim_", assay="RNA")
  }
  x
}
seGlmPCA.noweight <- function(x, ...){
  seGlmPCA(x, weight.by.var=FALSE, ...)
}

sceDR2seurat <- function(embeddings, object, name){
    if (is.null(rownames(embeddings))){
      rownames(embeddings) <- colnames(object)
    }
    key <- gsub(pattern = "[[:digit:]]", replacement = "_", 
                x = colnames(embeddings)[1])
    if (length(x = key) == 0) key <- paste0(name, "_")
    colnames(embeddings) <- paste0(key, 1:ncol(embeddings))
    object[[name]] <- CreateDimReducObject(embeddings = embeddings, key = key, 
                                           assay = DefaultAssay(object))
    object
}

clust.seurat <- function(x, rd=NULL, k=20, steps=8, dims=50, seed.use=1234, min.size=0, resolution=0.8){
  dims <- min(dims,ncol(x[["pca"]]@cell.embeddings))
  x <- FindNeighbors(x, k.param=k, dims=1:dims, verbose=FALSE)
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
  data("ctrlgenes", package="pipeComp")
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
  library(intrinsicDimension)
  x <- se[["pca"]]@cell.embeddings
  x <- switch(method,
         essLocal.a=essLocalDimEst(x),
         essLocal.b=essLocalDimEst(x, ver="b"),
         pcaLocal.FO=pcaLocalDimEst(x,ver="FO"),
         pcaLocal.fan=pcaLocalDimEst(x, ver="fan"),
         pcaLocal.maxgap=pcaLocalDimEst(x, ver="maxgap"),
         maxLikGlobal=maxLikGlobalDimEst(x, k=20, unbiased=TRUE),
         pcaOtpmPointwise.max=pcaOtpmPointwiseDimEst(x,N=10),
         elbow=farthestPoint(Stdev(se, "pca"))-1,
         fisherSeparability=FisherSeparability(x),
         scran.denoisePCA=scran.ndims.wrapper(se),
         jackstraw.elbow=js.wrapper(se,n.dims=ncol(x)-1,ret="ndims")
  )
  if(is.list(x) && "dim.est" %in% names(x)) x <- max(x$dim.est)
  round(x)
}

js.wrapper <- function(so, n.dims=50, n.rep=500, doplot=TRUE, ret=c("Seurat","pvalues","ndims")){
  ret <- match.arg(ret)
  so <- JackStraw(so, dims = n.dims, num.replicate=n.rep, verbose=FALSE)
  so <- ScoreJackStraw(so, dims = 1:n.dims, verbose=FALSE)
  if(ret=="pvalues") return( Reductions(so,"pca")@jackstraw$overall.p.values[,2] )
  if(ret=="Seurat") return( so )
  y <- so[["pca"]]@jackstraw$overall.p.values[,2]
  nzeros <- which(y>0)[1]-1
  y <- y[-1*(1:nzeros)]
  farthestPoint(-log10(y))+nzeros
}

scran.ndims.wrapper <- function(so){
  library(SingleCellExperiment)
  library(scran)
  sce <- SingleCellExperiment(list(counts=GetAssayData(so,"counts"), logcounts=GetAssayData(so, "data")))[VariableFeatures(so),]
  pcs <- getDenoisedPCs(sce, technical=modelGeneVar(sce))
  ncol(pcs$components)
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
  vars <- .translateFilterVars(strsplit(x,",",fixed=T)[[1]])
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
  x <- strsplit(x,"%",fixed=T)
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
#' @param dirs A vector of directions (higher/lower/both)
#'
#' @return A vector of filtering strings.
getFilterStrings <- function(mads=c(2,2.5,3,5), times=1:2, dirs=c("higher","both")){
  varCombs <- c("","mt","lcounts","mt,lcounts","lfeat,lcounts","mt,lfeat","mt,lfeat,lcounts",
                "mt,lfeat,lcounts,top50","mt,top50","lcounts,top50")
  v2 <- varCombs[grep("lfeat|lcounts", varCombs)]
  v2 <- gsub("lfeat","feat", v2)
  v2 <- gsub("lcounts","counts", v2)
  varCombs <- c(varCombs,v2)
  mads<- c(2, 2.5, 3, 5)
  varCombs <- unlist(lapply(strsplit(varCombs,",",fixed=T), dir=dirs, mads=mads, FUN=function(x, dir, mads){
    y <- expand.grid(lapply(x, y=dir, sep="%", FUN=paste))
    apply(y,1,collapse=",",FUN=paste)
  }))
  eg <- expand.grid(varCombs, mads, times)
  nbVars <- sapply(strsplit(as.character(eg[,1]),","),FUN=length)
  eg <- eg[which(nbVars>=eg[,3]),]
  as.character(apply(eg, 1, collapse="_", FUN=paste))
}

doublet.scDblFinder <- function(x){
  library(scDblFinder)
  x <- scDblFinder(x, verbose=FALSE)
  x[,which(x$scDblFinder.class!="doublet")]
}

doublet.scds <- function(x){
  library(scds)
  x <- cxds_bcds_hybrid(x, list(verb=FALSE), list(verb=FALSE))
  dbn <- ceiling((0.01 * ncol(x)/1000)*ncol(x))
  o <- order(x$hybrid_score, decreasing=TRUE)[seq_len(dbn)]
  x[,-o]
}

farthestPoint <- function(y, x=NULL){
  if(is.null(x)) x <- 1:length(y)
  d <- apply( cbind(x,y), 1, 
              a=c(1,y[1]), b=c(length(y),rev(y)[1]), 
              FUN=function(y, a, b){
                v1 <- a-b
                v2 <- y-a
                abs(det(cbind(v1,v2)))/sqrt(sum(v1*v1))
              })
  order(d,decreasing=T)[1]
}


#' clust.scran
#' 
#' A wrapper to use scran-based clustering.
#'
#' @param ds An object of class `SingleCellExperiment` or `Seurat`.
#' @param rd The name of the dimensionality reduction to use, or a logical 
#' indicating whether to use a reduced space.
#' @param method The method, either `fast_greedy` or `walktrap`.
#' @param k number of NN to consider, default 20.
#' @param steps number of steps for the random walk (walktrap method), default 8.
#' @param dims The (maximum) number of dimensions to use.
#' @param nthreads The number of threads, default 1.
#' @param seed.use not used.
#' @param min.size Minimum size of a cluster (default 50)
#' @param resolution Ignored; for consistency with `clust.seurat`
#' @param neighbor.method Passed to scran
#' @param graph.type "snn.rank", "snn.number", or "knn".
#'
#' @return A factor vector of cluster IDs, with cell names as names.
#' 
#' @export
clust.scran <- function(ds, rd=TRUE, method="walktrap", 
                        graph.type=c("snn.rank","snn.number","knn"),
                        neighbor.method=c("Kmknn","Vptree","Annoy","Hnsw"), 
                        k=20, steps=8, dims=50, nthreads=1, seed.use=NULL, 
                        min.size=5, resolution=NULL){
  graph.type <- match.arg(graph.type)
  weighting <- match.arg(weighting)
  neighbor.method <- switsch(match.arg(neighbor.method),
                             Kmknn=KmknnParam(),
                             Vptree=VptreeParam(),
                             Annoy=AnnoyParam(),
                             Hnsw=HnswParam())
  if(is(ds,"Seurat")){
    ds <- SingleCellExperiment(
      assays=list(
        counts=GetAssayData(ds, assay = "RNA", slot="counts"),
        logcounts=GetAssayData(ds, slot="scale.data")
      ),
      colData=ds[[]],
      reducedDims=lapply(ds@reductions, FUN=function(x) x@cell.embeddings)
    )
  }
  if(is.logical(rd)){
    if(rd){
      if(length(reducedDim(ds))>0){
        rd <- names(reducedDim(ds))[[1]]
        dims <- min(dims, ncol(reducedDim(ds, rd)))
      }else{
        rd <- NULL
      }
    }else{
      dims <- NA
      rd <- NULL
    }
  }
  BPPARAM <- ifelse(nthreads>1,MulticoreParam(nthreads),SerialParam())
  if(graph.type=="knn"){
    g <- scran::buildKNNGraph(ds, BPPARAM=BPPARAM, BNPARAM=neighbor.method, 
                              use.dimred=rd, k=k, d=dims)
  }else{
    weighting <- ifelse(graph.type=="snn.rank", "rank", "number")
    g <- scran::buildSNNGraph(ds, BPPARAM=BPPARAM, BNPARAM=neighbor.method, 
                              type=weighting, use.dimred=rd, k=k, d=dims,)    
  }
  if(method=="walktrap"){
    cl <- igraph::cluster_walktrap(g, steps=steps)$membership
  }else{
    cl <- igraph::cluster_fast_greedy(g)$membership
  }
  if(min.size>0) cl <- scran:::.merge_closest_graph(g, cl, min.size=min.size)
  names(cl) <- colnames(ds)
  as.factor(cl)
}

clust.scran.fg <- function(ds, ...){
  clust.scran(ds, method="fastq_greedy", ...)
}
