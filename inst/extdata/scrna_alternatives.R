#_______________________________________________________________________________
# DATASET PREPARATION ----------------------------------------------------------


#' add_meta
#'
#' Adds standard metadata.
#'
#' @param ds An object of class `SingleCellExperiment`
#'
#' @return The `SingleCellExperiment` object with updated metadata.
add_meta <- function(ds){
  library(scater)
  library(Matrix)
  data("ctrlgenes", package = "pipeComp")
  ## detect if row.names contain ensembl ids:
  if(any(grepl("^ENSG|^ENSMUSG",head(row.names(ds),n=100)))){
    ## we assume ^ENSEMBL\.whatever rownames
    cg <- lapply(c("Mt", "coding", "ribo"), FUN=function(x){
      union(ctrlgenes[[1]]$ensembl[[x]], ctrlgenes[[2]]$ensembl[[x]])
    })
    g <- sapply(strsplit(row.names(ds),".",fixed=TRUE),FUN=function(x) x[[1]])
  }else{
    # we assume HGNC/MGI symbols
    g <- row.names(ds)
    cg <- lapply(c("Mt", "coding", "ribo"), FUN=function(x){
      union(ctrlgenes[[1]]$symbols[[x]], ctrlgenes[[2]]$symbols[[x]])
    })
  }
  fc <- lapply(cg,g=g, FUN=function(x,g) g %in% x)
  names(fc) <- c("Mt","coding","ribosomal")
  if(exists("addQCPerCell")){
    ds <- addQCPerCell(ds, subsets=fc, percent_top=c(20,50,100,200))
  }else{
    ds <- addPerCellQC(ds, subsets=fc, percent_top=c(20,50,100,200))
  }
  ds$total_features <- ds$detected
  ds$log10_total_features <- log10(ds$detected)
  ds$total_counts <- ds$sum
  ds$log10_total_counts <- log10(ds$sum+1)
  ds$featcount_ratio <- ds$log10_total_counts/ds$log10_total_features
  ds$featcount_dist <- getFeatCountDist(ds)
  ds$pct_counts_top_50_features <- ds$percent_top_50
  for(f in names(fc)) 
    ds[[paste0("pct_",f)]] <- ds[[paste0("subsets_",f,"_percent")]]
  ds
}

#' getFeatCountDist
#' 
#' Returns the difference to the expected ratio of counts and number of features.
#'
#' @param df a cell metadata data.frame, or an object of class `SingleCellExperiment`
#' @param do.plot Logical; whether to plot the count/feature relationship.
#' @param linear Logical; whether to model the relationship with a linear model (default 
#' TRUE), rather than a loess.
#'
#' @return A vector of differences.
getFeatCountDist <- function(df, do.plot=FALSE, linear=TRUE){
  if(is(df,"SingleCellExperiment")) df <- as.data.frame(colData(df))
  if(linear){
    mod <- lm(df$log10_total_features~df$log10_total_counts)
  }else{
    mod <- loess(df$log10_total_features~df$log10_total_counts)
  }
  pred <- predict(mod, newdata=data.frame(log10_total_counts=df$log10_total_counts))
  df$diff <- df$log10_total_features - pred
  if(do.plot){
    library(ggplot2)
    ggplot(df, aes(x=total_counts, y=total_features, colour=diff)) + 
      geom_point() + geom_smooth(method = "loess", col="black")
  }
  df$diff
}


#' compute_all_gene_info
#' 
#' Populates the rowData of a SCE with various measures of variability, as well as the 
#' proportion of variance explained by clusters.
#'
#' @param sce An object of class `SingleCellExperiment`.
#'
#' @return The updated object.
compute_all_gene_info <- function(sce){
  library(variancePartition)
  library(sctransform)
  library(Seurat)
  cd <- as.data.frame(colData(sce))[,c("phenoid"),drop=FALSE]
  en <- log(t(1000*t(1+counts(sce))/colSums(counts(sce))))
  vp1 <- fitExtractVarPartModel(en, ~phenoid, data=cd)
  vst.out <- vst(counts(sce), colData(sce))
  vste <- vst.out$y+min(vst.out$y)
  vp2 <- fitExtractVarPartModel(vste,~phenoid, data=cd)
  vp <- data.frame(row.names=row.names(vp1), lognorm.varExp=vp1[,1], vst.varExp=vp2[row.names(vp1),1])
  source(system.file("extdata", "willtownes_scrna2019_utils.R", package="pipeComp"))
  gi <- compute_gene_info(counts(sce),gmeta=rowData(sce),mod="poisson")
  vp$deviance <- gi[row.names(vp),"deviance"]
  vp$total_counts <- gi[row.names(vp),"total_counts"]
  se <- .seuratFeatureVariability(sce, vst.out=vst.out)
  colnames(se) <- paste("seurat",colnames(se),sep=".")
  vp <- cbind(vp, se[row.names(vp),])
  RD <- as.data.frame(rowData(sce))
  RD <- RD[,setdiff(colnames(RD), colnames(vp))]
  rowData(sce) <- cbind(RD, vp[row.names(RD),])
  sce
}

.seuratFeatureVariability <- function(sce, vst.out=NULL){
  library(Seurat)
  library(SingleCellExperiment)
  seurat <- CreateSeuratObject( counts(sce), min.cells=0, min.features=0, 
                                project="scRNAseq" )
  seurat <- NormalizeData(seurat, display.progress=FALSE)
  seurat <- ScaleData(seurat, display.progress=FALSE)
  seurat <- FindVariableFeatures(seurat, selection.method="dispersion",
                                 nfeatures=nrow(seurat), display.progress=FALSE)
  seurat <- FindVariableFeatures(seurat, selection.method="vst",
                                 nfeatures=nrow(seurat), display.progress=FALSE)
  mf <- seurat@assays$RNA@meta.features
  if(is.null(vst.out)){
    library(sctransform)
    vst.out <- vst(counts(sce))
  }
  rv <- Seurat:::RowVar(vst.out$y)
  names(rv) <- row.names(vst.out$y)
  mf$res.var <- rv[row.names(mf)]
  mf
}


#' getDevianceExplained
#'
#' @param sce An object of class `SingleCellExperiment`
#' @param form.full The formula for the full model
#' @param form.null The formula for the reduced model (default `~lszie`, i.e. library size)
#' @param tagwise Logical; whether to run tagwise dispersion.
#'
#' @return The proportion of deviance (in the reduced model) explained by the full model.
#' @export
getDevianceExplained <- function(sce, form.full=~lsize+phenoid, form.null=~lsize, tagwise=TRUE){
  library(edgeR)
  sce$lsize <- log(colSums(counts(sce)))
  dds <- DGEList(as.matrix(counts(sce)))
  dds$samples$lib.size <- 1
  CD <- as.data.frame(colData(sce))
  mm <- model.matrix(form.full, data=CD)
  mm0 <- model.matrix(form.null, data=CD)
  dds <- estimateDisp(dds, mm, tagwise=tagwise)
  fit <- glmFit(dds, mm)
  fit0 <- glmFit(dds, mm0)
  de <- (deviance(fit0)-deviance(fit))/deviance(fit0)
  de[which(de<0)] <- 0
  return( de )
}



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
                     outlier.times=1
){
  vars <- vars[names(vars) %in% colnames(colData(x))]
  if(length(vars)==0){
    return( x[Matrix::rowSums(counts(x) > 0) >= min.cells,
              Matrix::colSums(counts(x) > 0) >= min.features] )
  } 
  out <- unlist(lapply(names(vars), o=x, v=vars, nm=nmads, FUN=function(x,o,v,nm){
    tryCatch(
      which(isOutlier(o[[x]], 
                      nmads=nm, 
                      log=FALSE, 
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
  x[Matrix::rowSums(counts(x) > 0) >= min.cells,
    Matrix::colSums(counts(x) > 0) >= min.features]
}

#' applyFilterString
#'
#' @param sce A SingleCellExperiment object.
#' @param filterstring A filtering string.
#'
#' @return A Seurat object.
#' @export
applyFilterString <- function(sce, filterstring){
  x <- strsplit(filterstring,"_",fixed=TRUE)[[1]]
  mads <- as.numeric(x[[2]])
  vars <- .translateFilterVars(strsplit(x,",",fixed=TRUE)[[1]])
  otimes <- ifelse(is.null(x[[3]]),1,as.numeric(x[[3]]))
  filt.mad(sce, nmads=mads, vars=vars, outlier.times=otimes)
}


filt.lenient <- function(x){  
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
  x[Matrix::rowSums(counts(x) > 0) >= 10, Matrix::colSums(counts(x) > 0) >= 0]
}

filt.default <- function(x, times=2){ 
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
  x[Matrix::rowSums(counts(x) > 0) >= 10, Matrix::colSums(counts(x) > 0) >= 0]
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
  x <- strsplit(x,".",fixed=TRUE)
  y <- sapply(x,FUN=function(x){ if(length(x)==1) return("both"); x[[2]] })
  names(y) <- sapply(x,vars=vars,FUN=function(x, vars) vars[x[[1]]])
  y
}

filt.pca <- function(x, vars=NULL){ 
  suppressPackageStartupMessages(library(scater))
  x <- runPCA(x, use_coldata=TRUE, detect_outliers=TRUE, selected_variables=vars)
  x[Matrix::rowSums(counts(x) > 0) >= 10,!x$outlier]
}

filt.pca2 <- function(x){
  filt.pca(x, vars=c("log10_total_counts", "log10_total_features", "pct_counts_Mt", "pct_counts_in_top_50_features"))
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
    dat <- dat[row.names(x),]
    logcounts(dat) <- GetAssayData(x, assay = "RNA", slot = "scale.data")
    return(dat) 
  }
}

norm.scran <- function(x, vars=NULL, noscale=TRUE, min.mean=1){
  suppressPackageStartupMessages({
    library(scran)
    library(Seurat)
    library(SingleCellExperiment)
  })
  if(is(x,"Seurat")){
    a <- GetAssayData(x, assay = "RNA", slot = "counts")
  }else{
    a <- counts(x)
  }
  a <- SingleCellExperiment(assays=list(counts=a))
  clusters <- quickCluster(a, min.mean=min.mean, min.size=50)
  a <- computeSumFactors(a, min.mean=min.mean, clusters=clusters)
  if(is.null(logNormCounts)){
    a <- scater::normalize(a)
  }else{
    a <- logNormCounts(a)
  }
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
  }else{
    a <- x
  }
  a <- SetAssayData(a, slot="counts", new.data=round(GetAssayData(a,slot="counts")))
  suppressPackageStartupMessages(library(sctransform))
  a <- SCTransform(a, vars.to.regress=vars, verbose=FALSE, 
                   variable.features.n=variable.features.n, 
                   return.only.var.genes=FALSE)
  a@misc$vst.var.feat <- VariableFeatures(a)
  if(is(x,"Seurat")) return(a)
  metadata(x)$vst.var.feat <- metadata(x)$VariableFeats <- VariableFeatures(a)
  logcounts(x) <- GetAssayData(a, "data", "SCT")
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

#' norm.scVI
#'
#' A function calling a python wrapper (`scVI.py`) around `scVI` normalization, adapted from the the 'Basic usage' Jupyter notebook (https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/basic_tutorial.ipynb). Note that the function will create a temporary csv file for the intermediate storage of the input count matrix, needed by `scVI`.  
#' 
#' 
#' @param x A SCE or Seurat object.
#' @param py_script Location of the python script
#' @param py_path Optional. If scVI was installed in a specific python bin, pass here the path to it. 
#' @param train_size Size of training set. Default to 0.8 but tutorial however recommends to use 1. 
#' @param n_cores N. cores
#' 
#' @return An object of the same class as `x` with updated slots.

norm.scVI <- function(x, py_script = system.file("extdata", "scVI.py", package="pipeComp"), py_path = NULL, n_cores = 1L, train_size = 1) {
  n_cores <- as.integer(n_cores)
  suppressPackageStartupMessages(library(reticulate))
  if (length(py_path)>0) use_python(py_path ,required=TRUE)
  trysource <- try(source_python(py_script))
  if (class(trysource) == "try-error") stop("Cannot source the python wrapper.") 
  tfile <- tempfile(fileext=".csv", tmpdir = ".")
  if (is(x, "Seurat")) dat <- GetAssayData(x, assay = "RNA", slot = "counts") else dat <- counts(x)
  write.csv(dat, tfile)
  out <- scVI_norm(csv_file = tfile, csv_path = ".", n_cores = n_cores, 
                   train_size = train_size)
  val <- t(out[[1]])
  gnames <- as.character(out[[2]])
  file.remove(tfile)
  dimnames(val) <- list(gnames, colnames(dat))
  x <- x[gnames, ]
  if(is(x,"Seurat")){ 
    x <- SetAssayData(x, slot="data", new.data=as.matrix(val))
    x <- SetAssayData(x, slot="scale.data", as.matrix(val))
  }else{
    logcounts(x) <- as.matrix(val)
  }
  x
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
  classes <- match.arg(gsub("ribosomal","ribo",classes), c("Mt","coding","ribo"), several.ok=TRUE)
  data("ctrlgenes", package="pipeComp")
  go <- g
  if(any(grepl("^ENSG|^ENSMUSG",head(g,n=100)))){
    ## we assume ^ENSEMBL\.whatever rownames
    cg <- lapply(classes, FUN=function(x){ union(ctrlgenes[[1]]$ensembl[[x]], ctrlgenes[[2]]$ensembl[[x]]) })
    g <- sapply(strsplit(g,".",fixed=TRUE),FUN=function(x) x[[1]])
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
  if(is(dat,"Seurat")) return(a)
  metadata(dat)$VariableFeats <- VariableFeatures(a)
  dat
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
  x <- strsplit(selstring,":",fixed=TRUE)[[1]]
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
    a@misc$rowData <- as.data.frame(rowData(dat))
  }
  e <- a@misc$rowData
  if(is(e,"list")) e <- as.data.frame(e)
  e <- e[row.names(a),f]
  VariableFeatures(a) <- row.names(a)[order(e, decreasing=TRUE)[1:min(n,length(e))]]
  VariableFeatures(a) <- subsetFeatureByType(VariableFeatures(a), excl)
  if(is(dat, "Seurat")) return(a)
  metadata(dat)$VariableFeats <- VariableFeatures(a)
  dat
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
    dat <- sceWrap(x)
  } else {
    dat <- x
  }
  lcmin <- min(logcounts(dat))
  if(lcmin<0) logcounts(dat) <- logcounts(dat) - lcmin
  if(packageVersion("scran") >= "1.13"){
    var.stats <- modelGeneVar(dat)
    dat <- denoisePCA(dat, technical=var.stats, min.rank=2, max.rank=dims, 
                      subset.row=metadata(dat)$VariableFeats, BSPARAM=BSPARAM, 
                      ...)
  }else{
    td <- trendVar(dat, use.spikes=FALSE)
    dat <- denoisePCA(dat, technical=td$trend, min.rank=2, max.rank=dims, 
                      subset.row=metadata(dat)$VariableFeats, BSPARAM=BSPARAM, 
                      ...)
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
  colnames(e) <- gsub("dim_", "PC_", colnames(e))
  if(weight.by.var=="both" && length(dr$dev) %in% dim(e)){
    if(is(x, "Seurat")) {
      x[["glmpca"]] <- CreateDimReducObject(embeddings=e, key="PC_", assay="RNA")
      e <- t(t(e)*dr$d)
      x[["glmpca.wt"]] <- CreateDimReducObject(embeddings=e, key="PC_", assay="RNA")
    } else {
      reducedDim(x, "glmpca") <- e
      e <- t(t(e)*dr$d)
      reducedDim(x, "glmpca.wt") <- e
    }
  }else{
    if(weight.by.var && length(dr$dev) %in% dim(e)) e <- t(t(e)*dr$d)
    if(is(x, "Seurat")){
      x[["pca"]] <- CreateDimReducObject(embeddings=e, key="PC_", assay="RNA")
    } else {
      reducedDim(x, "PCA") <- e
    }
  }
  x
}
GlmPCA.noweight <- function(x, ...){
  GlmPCA(x, weight.by.var=FALSE, ...)
}

scran.runPCA <- function(x, dims=50){
  suppressPackageStartupMessages(library(scran))
  suppressPackageStartupMessages(library(BiocSingular))
  if(is(x, "Seurat")){
    dat <-  sceWrap(x)
  } else {
    dat <- x
  }
  dat <- scater::runPCA(dat, ncomponents = dims, 
                        subset_row=metadata(dat)$VariableFeats)
  if(is(x, "Seurat")) return(sceDR2seurat(reducedDim(dat, "PCA"), x, "pca")) else return(dat)
}

#' scVI.latent
#'
#' A function calling a python wrapper (`scVI.py`) around `scVI` low-dimensional latent space, adapted from the official tutorial (https://scvi.readthedocs.io/en/stable/tutorials/basic_tutorial.html). Note that the function will create a temporary csv file for the intermediate storage of the input count matrix, needed by `scVI`.  
#' 
#' 
#' @param x A SCE or Seurat object.
#' @param py_script Location of the python script
#' @param py_path Optional. If scVI was installed in a specific python bin, pass here the path to it. 
#' @param dims Number of dimensions to return. 
#' @param learning_rate Learning rate of the model. If the model is not training properly due to too high learning rate, it will be reduced consecutively a few times before early stop. 
#' @param n_cores N. cores
#' 
#' @return An object of the same class as `x` with updated slots. Note that scVI-LD initially returns unordered components. For convenience with the package, they are ordered by sdev and renamed 'PC'. 

scVI.latent <- function(x, dims = 50L, learning_rate = 1e-3, py_script = system.file("extdata", "scVI.py", package="pipeComp"), py_path = NULL, n_cores = 1L) {
  dims <- as.integer(dims)
  n_cores <- as.integer(n_cores)
  suppressPackageStartupMessages(library(reticulate))
  if (length(py_path)>0) use_python(py_path ,required=TRUE)
  trysource <- try(source_python(py_script))
  if (class(trysource) == "try-error") stop("Cannot source the python wrapper.") 
  tfile <- tempfile(fileext=".csv", tmpdir = ".")
  if (is(x, "Seurat")) {
    dat <- GetAssayData(x, assay = "RNA", slot = "counts")
    dat <- dat[VariableFeatures(x), ]
  } else {
    
    dat <- counts(x)
    dat <- dat[metadata(x)$VariableFeat, ]
  } 
  write.csv(dat, tfile)
  val <- try(scVI_latent(csv_file = tfile, csv_path = ".", n_cores = n_cores, lr = learning_rate))
  # Error with some dataset; "Loss was NaN 10 consecutive times: the model is not training properly. Consider using a lower learning rate" --> reduce lr
  if (class(val) == "try-error") {
    ntry <- 0
    while(ntry < 6 & class(val) == "try-error") {
      ntry <- ntry + 1
      learning_rate <- learning_rate/10
      message("Downgrading learning rate to ", learning_rate)
      val <- try(scVI_latent(csv_file = tfile, csv_path = ".", n_cores = n_cores, lr = learning_rate))
    }
  }
  file.remove(tfile)
  rownames(val) <- colnames(dat)
  # rename
  colnames(val) <- paste0("PC_", 1:ncol(val))
  if(is(x, "Seurat")){
    x[["pca"]] <- CreateDimReducObject(embeddings=as.matrix(val), 
                                       key="PC_", assay="RNA")
  } else {
    reducedDim(x, "PCA") <- as.matrix(val)
  }
  x
}


#' scVI.LD
#'
#' A function calling a python wrapper (`scVI.py`) around `scVI` linear decoded, adapted from the the 'Basic usage' Jupyter notebook (https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/linear_decoder.ipynb). Note that the function will create a temporary csv file for the intermediate storage of the input count matrix, needed by `scVI`.  
#' 
#' 
#' @param x A SCE or Seurat object.
#' @param py_script Location of the python script
#' @param py_path Optional. If scVI was installed in a specific python bin, pass here the path to it. 
#' @param dims Number of dimensions to return. 
#' @param learning_rate Learning rate of the model. If the model is not training properly due to too high learning rate, it will be reduced consecutively a few times before early stop. 
#' @param n_cores N. cores
#' 
#' @return An object of the same class as `x` with updated slots. Note that scVI-LD initially returns unordered components. For convenience with the package, they are ordered by sdev and renamed 'PC'. 

scVI.LD <- function(x, dims = 50L, learning_rate = 1e-3, py_script = system.file("extdata", "scVI.py", package="pipeComp"), py_path = NULL, n_cores = 1L) {
  dims <- as.integer(dims)
  n_cores <- as.integer(n_cores)
  suppressPackageStartupMessages(library(reticulate))
  if (length(py_path)>0) use_python(py_path ,required=TRUE)
  trysource <- try(source_python(py_script))
  if (class(trysource) == "try-error") stop("Cannot source the python wrapper.") 
  tfile <- tempfile(fileext=".csv", tmpdir = ".")
  if (is(x, "Seurat")) {
    dat <- GetAssayData(x, assay = "RNA", slot = "counts")
    dat <- dat[VariableFeatures(x), ]
  } else {
    
    dat <- counts(x)
    dat <- dat[metadata(x)$VariableFeat, ]
  } 
  write.csv(dat, tfile)
  val <- try(scVI_ld(csv_file = tfile, ndims = dims, csv_path = ".", n_cores = n_cores, lr = learning_rate))
  # Error with some dataset; "Loss was NaN 10 consecutive times: the model is not training properly. Consider using a lower learning rate" --> reduce lr
  if (class(val) == "try-error") {
    ntry <- 0
    while(ntry < 6 & class(val) == "try-error") {
      ntry <- ntry + 1
      learning_rate <- learning_rate/10
      message("Downgrading learning rate to ", learning_rate)
      val <- try(scVI_ld(csv_file = tfile, ndims = dims, csv_path = ".", n_cores = n_cores, lr = learning_rate))
    }
  }
  file.remove(tfile)
  rownames(val) <- colnames(dat)
  # rename
  colnames(val) <- paste0("PC_", 1:ncol(val))
  if(is(x, "Seurat")){
    x[["pca"]] <- CreateDimReducObject(embeddings=as.matrix(val), 
                                       key="PC_", assay="RNA")
  } else {
    reducedDim(x, "PCA") <- as.matrix(val)
  }
  x
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


FisherSeparability <- function(PCAdims, py_script = system.file("extdata", "FisherSeparability.py", package="pipeComp")) {
  if(is(PCAdims, "Seurat")){
    PCAdims <- PCAdims[["pca"]]@cell.embeddings
  } else if(is(PCAdims, "SingleCellExperiment")){
    PCAdims <- reducedDim(dat, "PCA")
  }
  # Adapted from https://github.com/auranic/FisherSeparabilityAnalysis
  suppressPackageStartupMessages(library(reticulate))
  suppressPackageStartupMessages(library(Seurat))
  trysource <- try(source_python(py_script))
  if (is(trysource, "try-error")) stop("Cannot source 'FisherSeparability.py'. Make sure:\n1) You are running reticulate and redirecting to a valid Python3 bin/conda (use_python, use_conda)\n2) You have the following modules installed: numpy, math, sklearn.decomposition, seaborn, warnings, scipy.special, matplotlib, scipy.io\n3) You have 'FisherSeparability.py' in your current wd") 
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


js.wrapper <- function(dat, n.dims=NULL, n.rep=500, doplot=FALSE, ret=c("ndims", "Seurat", "sce","pvalues")){
  ret <- match.arg(ret)
  if (!is(dat, "Seurat")) x <- seWrap(dat) else x <- dat
  if(is.null(n.dims)) n.dims <- ncol(Reductions(x, "pca")@cell.embeddings)-1
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
  x <- strsplit(x,"%",fixed=TRUE)
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
  varCombs <- unlist(lapply(strsplit(varCombs,",",fixed=TRUE), dir=dirs, mads=mads, FUN=function(x, dir, mads){
    y <- expand.grid(lapply(x, y=dir, sep="%", FUN=paste))
    apply(y,1,collapse=",",FUN=paste)
  }))
  eg <- expand.grid(varCombs, mads, times)
  nbVars <- sapply(strsplit(as.character(eg[,1]),","),FUN=length)
  eg <- eg[which(nbVars>=eg[,3]),]
  as.character(apply(eg, 1, collapse="_", FUN=paste))
}

doublet.scDblFinder <- function(x, trans = "scran"){
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
    ds <- sceWrap(ds)
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
    if(!is.null(rd)) reducedDim(ds, rd) <- reducedDim(ds, rd)[,seq_len(dims)]
    g <- scran::buildKNNGraph(ds, BPPARAM=BPPARAM, BNPARAM=neighbor.method, 
                              use.dimred=rd, k=k)
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
