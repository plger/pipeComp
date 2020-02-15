
#_______________________________________________________________________________
# MISC WRAPPER -----------------------------------------------------------------


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

seWrap <- function(sce, min.cells=0, min.features=0){
  suppressPackageStartupMessages(library(Seurat))
  se <- CreateSeuratObject( counts=counts(sce), 
                            min.cells=min.cells, 
                            min.features=min.features, 
                            meta.data=as.data.frame(colData(sce)), 
                            project = "scRNAseq" )
  se@misc$rowData <- as.data.frame(rowData(sce))
  if("logcounts" %in%  assayNames(sce)){
    se <- ScaleData(se, verbose = FALSE)
    se@assays$RNA@data <- logcounts(sce)
  } 
  if(!is.null(metadata(sce)$VariableFeats)) VariableFeatures(se) <- metadata(sce)$VariableFeats
  if(length(reducedDimNames(sce)) != 0) se[["pca"]] <- CreateDimReducObject(embeddings=reducedDim(sce), key="PC_", assay="RNA")
  se
}

sceWrap <- function(seu) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Seurat)
  })
  sce <- SingleCellExperiment(list(counts = GetAssayData(seu, assay = "RNA", slot = "counts")), 
                              colData = seu[[]])
  rowData(sce) <- seu@misc$rowData
  sce
}

none <- function(x) x


#_______________________________________________________________________________
# FILTERING  -------------------------------------------------------------------

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
filt.mad <- function(x, nmads=3, min.cells=10, min.features=100,
                     vars=c( "pct_mt"="higher", 
                             "log10_total_features"="both", 
                             "log10_total_counts"="both", 
                             "pct_counts_top_50_features"="higher"
                     ),
                     outlier.times=1, pipeClass = "seurat"
){
  vars <- vars[names(vars) %in% colnames(colData(x))]
  if(length(vars)==0){
    if (pipeClass == "seurat") {
      return(seWrap(x, min.cells=min.cells, min.features=min.features))
    } else {
      return(x[Matrix::rowSums(counts(x) > 0) >= min.cells,
        Matrix::colSums(counts(x) > 0) >= min.features])
    }
  } 
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
  if (pipeClass == "seurat") {
    seWrap(x, min.cells=min.cells, min.features=min.features)
  } else {
    x[Matrix::rowSums(counts(x) > 0) >= min.cells,
      Matrix::colSums(counts(x) > 0) >= min.features]
  }
}

#' applyFilterString
#'
#' @param sce A SingleCellExperiment object.
#' @param filterstring A filtering string.
#'
#' @return A Seurat object.
#' @export
applyFilterString <- function(sce, filterstring, pipeClass){
  x <- strsplit(filterstring,"_",fixed=T)[[1]]
  mads <- as.numeric(x[[2]])
  vars <- .translateFilterVars(strsplit(x,",",fixed=T)[[1]])
  otimes <- ifelse(is.null(x[[3]]),1,as.numeric(x[[3]]))
  filt.mad(sce, pipeClass, nmads=mads, vars=vars, outlier.times=otimes)
}


filt.lenient <- function(x, pipeClass = "seurat"){  
  suppressPackageStartupMessages(library(scater))
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
  if (pipeClass == "seurat")  seWrap(x, min.cells = 10) else x[Matrix::rowSums(counts(x) > 0) >= 10, Matrix::colSums(counts(x) > 0) >= 0]
}

filt.default <- function(x, times=2, pipeClass = "seurat"){ 
  suppressPackageStartupMessages(library(scater))
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
  if (pipeClass == "seurat") seWrap(x, min.cells = 10) else x[Matrix::rowSums(counts(x) > 0) >= 10, Matrix::colSums(counts(x) > 0) >= 0]
}

filt.stringent <- function(x, pipeClass){
  filt.default(x,1, pipeClass)
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

filt.pca <- function(x, vars=NULL, pipeClass = "seurat"){ 
  suppressPackageStartupMessages(library(scater))
  x <- runPCA(x, use_coldata=TRUE, detect_outliers=TRUE, selected_variables=vars)
  if (pipeClass == "seurat") seWrap(x[,!x$outlier], min.cells = 10) else x[Matrix::rowSums(counts(x) > 0) >= 10,!x$outlier]
}

filt.pca2 <- function(x, pipeClass = "seurat"){
  filt.pca(x, vars=c("log10_total_counts", "log10_total_features", "pct_counts_Mt", "pct_counts_in_top_50_features"), pipeClass = pipeClass)
}


#_______________________________________________________________________________
# NORMALIZATION  ---------------------------------------------------------------

norm.seurat <- function(dat, vars=NULL, noscale=FALSE){
  suppressPackageStartupMessages(library(Seurat))
  if (!is(dat, "Seurat")) {
    x <- seWrap(dat)
  } else {
    x <- dat
  }
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
  if (is(dat, "Seurat")){
    return(x)
  } else {
    logcounts(dat) <- GetAssayData(x, assay = "RNA", slot = "data")
    return(dat) 
  }
}

norm.scran <- function(x, vars=NULL, noscale=TRUE, min.mean=1){
  suppressPackageStartupMessages({
    library(scran)
    library(Seurat)
  })
  if(is(x,"Seurat")){
    a <- GetAssayData(x, assay = "RNA", slot = "counts")
  }else{
    a <- counts(x)
  }
  a <- SingleCellExperiment(assays=list(counts=a))
  clusters <- quickCluster(a, min.mean=min.mean, min.size=50)
  a <- computeSumFactors(a, min.mean=min.mean, clusters=clusters)
  a <- logNormCounts(a)
  if(is(x,"Seurat")){ 
    x <- SetAssayData(x, slot="data", new.data=logcounts(a))
    if(noscale){
      x <- SetAssayData(x, slot="scale.data", as.matrix(GetAssayData(x)))
    }else{
      x <- ScaleData(x, verbose=FALSE, vars.to.regress=vars)
    }
  }else{
    if(!noscale) a <- t(scale(t(logcounts(a)))) else a <- logcounts(a)
    logcounts(x) <- a
  }
  x
}
norm.scran.scaled <- function(x, ...){
  norm.scran(x, noscale=FALSE, ...)
}

norm.none <- function(x, vars=NULL, noscale=TRUE){
  if(is(x,"Seurat")){
    a <- GetAssayData(x, slot="counts")
  }else{
    a <- counts(x)
  }
  a <- log1p(a)
  if(is(x,"Seurat")){
    x <- SetAssayData(x, slot="data", a)
    if(noscale){
      x <- SetAssayData(x, slot="scale.data", as.matrix(GetAssayData(x)))
    }else{
      x <- ScaleData(x, verbose=FALSE, vars.to.regress=vars)
    }
  } else {
    if(!noscale) a <- t(scale(t(a)))
    logcounts(x) <- a
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
#' @param x A Seurat/ SCE object.
#' @param vars A vector of variables to regress when scaling (default none)
#' @param noscale Ignored.
#' @param variable.features.n Passed to `SCTransform`, default 5000
#'
#' @return A Seurat/ SCE object with updated data slot.
norm.seuratvst <- function(x, vars=NULL, noscale=FALSE, variable.features.n=5000){
  if(!is(x,"Seurat")){
    a <- seWrap(x)
  } else a <- x
  suppressPackageStartupMessages(library(sctransform))
  a <- SCTransform(a, vars.to.regress=vars, verbose=FALSE, 
                   variable.features.n=variable.features.n, 
                   return.only.var.genes=FALSE)
  if(is(x,"Seurat")){
    a@misc$vst.var.feat <- VariableFeatures(a)
    a
  } else {
    metadata(x)$vst.var.feat <- VariableFeatures(a)
    logcounts(x) <- GetAssayData(a, "data", "SCT")
    x
  }
  
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
  suppressPackageStartupMessages(library(SCnorm))
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

# ______________________________________________________________________________
# FEATURE SELECTION ------------------------------------------------------------

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

sel.vst <- function(dat, n=2000, excl=c()){
  if(!is(dat,"Seurat")){
    a <- seWrap(dat)
  } else a <- dat
  if(!is.null(a@misc$vst.var.feat)){
    VariableFeatures(a) <- a@misc$vst.var.feat[1:min(n,length(a@misc$vst.var.feat))]
  }else{
    a <- FindVariableFeatures(a, nfeatures=n)
  }
  VariableFeatures(a) <- subsetFeatureByType(VariableFeatures(a), excl)
  if(is(dat,"Seurat")){
    a
  } else {
    metadata(dat)$VariableFeats <- VariableFeatures(a)
    dat
  }
}

#' applySelString
#' 
#' Applies a given selection string.
#'
#' @param se A Seurat or SCE object.
#' @param selstring A rowData variable to use, or a selection string, e.g. 'vst:2000:coding_rmMt_rmribo'.
#' @param n The number of genes to select (ignored if selstring is a full selection string).
#'
#' @return A filtered Seurat or SCE object.
#' @export
applySelString <- function(dat, selstring, n=2000){
  #vst:2000:coding_rmMt_rmribo
  x <- strsplit(selstring,":",fixed=T)[[1]]
  excl <- c()
  if(length(x)>2 && x[3]!="") excl <- gsub("rm","",strsplit(x[3], "_")[[1]])
  fn <- paste0("sel.",x[1])
  if(exists(x[2])) n <- as.numeric(x[2])
  if(exists(fn) && is.function(get(fn))) return(get(fn)(dat, n=n, excl=excl))
  sel.fromField(dat, x[1], n, excl)
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
  if(is(dat, "Seurat")) {
    if(is.null(dat@misc$rowData[[f]])) return(NULL)
    a <- dat
  } else {
    if(is.null(rowData(dat)[[f]])) return(NULL)
    a <- seWrap(dat)
    a@misc$rowData <- rowData(dat)
  }
  e <- a@misc$rowData[row.names(a),f]
  VariableFeatures(a) <- row.names(a)[order(e, decreasing=T)[1:min(n,length(e))]]
  VariableFeatures(a) <- subsetFeatureByType(VariableFeatures(a), excl)
  if(is(dat, "Seurat")) {
    a
  } else {
    metadata(dat)$VariableFeats <- VariableFeatures(a)
    dat
  }
}

sel.deviance <- function(x, n=2000, excl=c()){
  sel.fromField(x, "deviance", n=n, excl=excl)
}
sel.expr <- function(x, n=2000, excl=c()){
  sel.fromField(x, "total_counts", n=n, excl=excl)
}


#_______________________________________________________________________________
# DIMENSION REDUCTION ----------------------------------------------------------


seurat.pca <- function(x, dims=50, weight.by.var=TRUE, seed.use=42){
  if(is(x, "Seurat")){
    dat <- x
  } else {
    dat <- seWrap(x)
  }
  dat <- RunPCA(dat, features=VariableFeatures(dat), verbose=FALSE, 
                 weight.by.var=weight.by.var, npcs=dims, seed.use = seed.use)
  if(is(x, "Seurat")){
    dat
  } else {
    reducedDim(x, "PCA") <- Reductions(dat, "pca")@cell.embeddings
    x
  }
}

seurat.pca.noweight <- function(x, dims=50, weight.by.var=FALSE, seed.use=42){
  if(is(x, "Seurat")){
    dat <- x
  } else {
    dat <- seWrap(x)
  }
  dat <- RunPCA(dat, features=VariableFeatures(dat), verbose=FALSE, 
                weight.by.var=weight.by.var, npcs=dims, seed.use = seed.use)
  if(is(x, "Seurat")){
    dat
  } else {
    reducedDim(x, "PCA") <- Reductions(dat, "pca")@cell.embeddings
    x
  }
}

scran.denoisePCA <- function(x, dims=50, pca.method=c("exact","irlba"), ...){
  suppressPackageStartupMessages(library(scran))
  suppressPackageStartupMessages(library(BiocSingular))
  BSPARAM <- switch(match.arg(pca.method),
                    exact=ExactParam(),
                    irlba=IrlbaParam() )
  if(is(x, "Seurat")){
    dat <-  SingleCellExperiment( list( counts=GetAssayData(x, slot="counts"),
                                        logcounts=GetAssayData(x, slot="data")) )
  } else {
    dat <- x
  }
  if(packageVersion("scran") >= "1.13"){
    var.stats <- modelGeneVar(dat)
    dat <- denoisePCA(dat, technical=var.stats, min.rank=2, max.rank=dims, BSPARAM=BSPARAM, ...)
  }else{
    td <- trendVar(dat, use.spikes=FALSE)
    dat <- denoisePCA(dat, technical=td$trend, min.rank=2, max.rank=dims, BSPARAM=BSPARAM, ...)
  }
  if(is(x, "Seurat")) return(sceDR2seurat(reducedDim(dat, "PCA"), x, "pca")) else return(dat)
}

GlmPCA <- function(x, weight.by.var=TRUE, dims=20){
  suppressPackageStartupMessages(library(glmpca))
  if(is(x, "Seurat")){
    dat <- x
  } else {
    dat <- seWrap(x)
  }
  dr <- glmpca(as.matrix(GetAssayData(dat, assay = "RNA", slot = "counts")[VariableFeatures(dat),]), dims)
  e <- as.matrix(dr$factors)
  colnames(e) <- gsub("dim","dim_",colnames(e))
  if(weight.by.var=="both" && length(dr$dev) %in% dim(e)){
    if(is(x, "Seurat")) {
      x[["glmpca"]] <- CreateDimReducObject(embeddings=e, key="dim_", assay="RNA")
      e <- t(t(e)*dr$d)
      x[["glmpca.wt"]] <- CreateDimReducObject(embeddings=e, key="dim_", assay="RNA")
    } else {
      reducedDim(x, "glmpca") <- e
      e <- t(t(e)*dr$d)
      reducedDim(x, "glmpca.wt") <- e
    }
  }else{
    if(weight.by.var && length(dr$dev) %in% dim(e)) e <- t(t(e)*dr$d)
    if(is(x, "Seurat")){
      x[["pca"]] <- CreateDimReducObject(embeddings=e, key="dim_", assay="RNA")
    } else {
      reducedDim(x, "PCA") <- e
    }
  }
  x
}
GlmPCA.noweight <- function(x, ...){
  GlmPCA(x, weight.by.var=FALSE, ...)
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



FisherSeparability <- function(PCAdims, py_script = system.file("extdata", "FisherSeparability.py", package="pipeComp")) {
  # Adapted from https://github.com/auranic/FisherSeparabilityAnalysis
  suppressPackageStartupMessages(library(reticulate))
  suppressPackageStartupMessages(library(Seurat))
  trysource <- try(source_python(py_script))
  if (class(trysource) == "try-error") stop("Cannot source 'FisherSeparability.py'. Make sure:\n1) You are running reticulate and redirecting to a valid Python3 bin/conda (use_python, use_conda)\n2) You have the following modules installed: numpy, math, sklearn.decomposition, seaborn, warnings, scipy.special, matplotlib, scipy.io\n3) You have 'FisherSeparability.py' in your current wd") 
  # Saving to numpy for correct coercion 
  np <- import("numpy")
  tdir <- "tempnpy"
  tfile <- tempfile(fileext=".npy", tmpdir = tdir)
  dir.create(tdir, showWarnings = FALSE)
  np$save(tfile, as.matrix(PCAdims))
  val <- SeparabilityAnalysis(tfile, ProducePlots = TRUE)
  val <- as.numeric(val)
  unlink(tdir, recursive = TRUE)
  return(val)
}


getDimensionality <- function(dat, method, maxDims=50){
  suppressPackageStartupMessages(library(intrinsicDimension))
  if(is(dat, "Seurat")){
    x <- dat[["pca"]]@cell.embeddings
    sdv <- Stdev(dat, "pca")
  } else {
    x <- reducedDim(dat, "PCA")
    sdv <- attr(reducedDim(dat, "PCA"), "percentVar")
  }
  x <- switch(method,
              essLocal.a=essLocalDimEst(x),
              essLocal.b=essLocalDimEst(x, ver="b"),
              pcaLocal.FO=pcaLocalDimEst(x,ver="FO"),
              pcaLocal.fan=pcaLocalDimEst(x, ver="fan"),
              pcaLocal.maxgap=pcaLocalDimEst(x, ver="maxgap"),
              maxLikGlobal=maxLikGlobalDimEst(x, k=20, unbiased=TRUE),
              pcaOtpmPointwise.max=pcaOtpmPointwiseDimEst(x,N=10),
              elbow=farthestPoint(sdv)-1,
              fisherSeparability=FisherSeparability(x),
              scran.denoisePCA=scran.ndims.wrapper(se),
              jackstraw.elbow=js.wrapper(dat,n.dims=ncol(dat)-1,ret="ndims")
  )
  if(is.list(x) && "dim.est" %in% names(x)) x <- max(x$dim.est)
  round(x)
}

js.wrapper <- function(dat, n.dims=50, n.rep=500, doplot=TRUE, ret=c("Seurat", "sce","pvalues","ndims")){
  ret <- match.arg(ret)
  if (!is(dat, "Seurat")) x <- seWrap(dat) else x <- dat
  x <- JackStraw(x, dims = n.dims, num.replicate=n.rep, verbose=FALSE)
  x <- ScoreJackStraw(x, dims = 1:n.dims, verbose=FALSE)
  if(ret=="pvalues") return( Reductions(x,"pca")@jackstraw$overall.p.values[,2] )
  if(ret=="Seurat") return(x)
  if(ret=="sce") {
    metadata(dat)$jackstraw <- Reductions(x,"pca")@jackstraw
    return(dat)
  }
  y <- x[["pca"]]@jackstraw$overall.p.values[,2]
  nzeros <- which(y>0)[1]-1
  y <- y[-1*(1:nzeros)]
  farthestPoint(-log10(y))+nzeros
}

scran.ndims.wrapper <- function(dat){
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scran)
  })
  if(is(dat, "Seurat")) {
    x <- SingleCellExperiment(list(counts=GetAssayData(dat,"counts"), 
                                     logcounts=GetAssayData(dat, "data")))[VariableFeatures(dat),]
  } else {
    x <- dat
  }
  pcs <- getDenoisedPCs(x, technical=modelGeneVar(x))
  ncol(pcs$components)
}


.tmads <- function(x, nbmads=2.5){
  x2 <- nbmads*median(abs(x-median(x)))
  median(x)+c(-x2,x2)
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
  suppressPackageStartupMessages(library(scDblFinder))
  x <- scDblFinder(x, verbose=FALSE)
  x[,which(x$scDblFinder.class!="doublet")]
}

doublet.scds <- function(x){
  suppressPackageStartupMessages(library(scds))
  x <- cxds_bcds_hybrid(x, list(verb=FALSE), list(verb=FALSE))
  dbn <- ceiling((0.01 * ncol(x)/1000)*ncol(x))
  o <- order(x$hybrid_score, decreasing=TRUE)[seq_len(dbn)]
  x[,-o]
}





#_______________________________________________________________________________
# CLUSTERING  ------------------------------------------------------------------

clust.seurat <- function(x, rd=NULL, k=20, steps=8, dims=50, seed.use=1234, min.size=0, resolution=0.8){
  if(is(x, "Seurat")) {
    dat <- x 
  } else {
    dat <- seWrap(x)
  }
  dims <- min(dims,ncol(dat[["pca"]]@cell.embeddings))
  dat <- FindNeighbors(dat, k.param=k, dims=1:dims, verbose=FALSE)
  dat <- FindClusters(dat, resolution=resolution, random.seed=seed.use, verbose=FALSE)
  Idents(dat)
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
  suppressPackageStartupMessages(library(BiocNeighbors))
  graph.type <- match.arg(graph.type)
  neighbor.method <- switch(match.arg(neighbor.method),
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
        rd <- reducedDimNames(ds)[1]
        dims <- min(dims, ncol(reducedDim(ds, rd)))
      }else{
        rd <- NULL
      }
    }else{
      dims <- NA
      rd <- NULL
    }
  }
  BPPARAM <- if(nthreads>1) MulticoreParam(nthreads) else SerialParam()
  if(graph.type=="knn"){
    g <- scran::buildKNNGraph(ds, BPPARAM=BPPARAM, BNPARAM=neighbor.method, 
                              use.dimred=rd, k=k, d=dims)
  }else{
    weighting <- ifelse(graph.type=="snn.rank", "rank", "number")
    
    if (is.null(rd)) {
      g <- scran::buildSNNGraph(ds, BPPARAM=BPPARAM, BNPARAM=neighbor.method, 
                                type=weighting, use.dimred=rd, k=k, d=dims)   
    } else {
      g <- scran::buildSNNGraph(ds, BPPARAM=BPPARAM, BNPARAM=neighbor.method, 
                                type=weighting, use.dimred=rd, k=k) 
    }
     
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
