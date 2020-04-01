#' getDimensionality
#' 
#' Returns the estimated intrinsic dimensionality of a dataset.
#'
#' @param dat A Seurat or SCE object with a pca embedding.
#' @param method The dimensionality method to use.
#' @param maxDims Deprecated and ignored.
#'
#' @return An integer.
#'
#' @importFrom Seurat Stdev
#' @importFrom SingleCellExperiment reducedDim
#' @import intrinsicDimension
getDimensionality <- function(dat, method, maxDims=NULL){
  if(is(dat, "Seurat")){
    x <- dat[["pca"]]@cell.embeddings
    sdv <- Stdev(dat, "pca")
  } else {
    if(method=="jackstraw.elbow")
      stop("the jackstraw.elbow method is only available for Seurat objects.")
    x <- reducedDim(dat, "PCA")
    sdv <- attr(reducedDim(dat, "PCA"), "percentVar")
  }
  # conversion for backward compatibility:
  conv <- c( "fisherSeparability"="FisherSeparability",
             "scran.denoisePCA"="scran.ndims.wrapper", 
             "jackstraw.elbow"="js.wrapper"
  )
  if(method %in% names(conv)) method <- as.character(conv[method])
  
  x <- switch(method,
              essLocal.a=essLocalDimEst(x),
              essLocal.b=essLocalDimEst(x, ver="b"),
              pcaLocal.FO=pcaLocalDimEst(x,ver="FO"),
              pcaLocal.fan=pcaLocalDimEst(x, ver="fan"),
              pcaLocal.maxgap=pcaLocalDimEst(x, ver="maxgap"),
              maxLikGlobal=maxLikGlobalDimEst(x, k=20, unbiased=TRUE),
              pcaOtpmPointwise.max=pcaOtpmPointwiseDimEst(x,N=10),
              elbow=farthestPoint(sdv)-1,
              ifelse( !is.function(method) && 
                        !( is.character(method) && 
                             is.function(method <- get(method)) ),
                      stop("Unknown dimensionality method!"),
                      method(dat) )
  )
  if(is.list(x) && "dim.est" %in% names(x)) x <- max(x$dim.est)
  as.integer(round(x))
}


# sce2se conversion
# not exported
#' @import SingleCellExperiment Seurat
#' @importFrom SummarizedExperiment assayNames
seWrap <- function(sce, min.cells=10, min.features=0){
  if(is(sce,"Seurat")) return(sce)
  if(!is(sce,"SingleCellExperiment")) stop("not a SingleCellExperiment!")
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
  if(!is.null(metadata(sce)$VariableFeats)) 
    VariableFeatures(se) <- metadata(sce)$VariableFeats
  if(length(reducedDimNames(sce)) != 0) 
    se[["pca"]] <- CreateDimReducObject(embeddings=reducedDim(sce), key="PC_", 
                                        assay="RNA")
  se
}

# se2sce conversion
# not exported
#' @import SingleCellExperiment Seurat
#' @importFrom SummarizedExperiment rowData<-
sceWrap <- function(seu) {
  if(is(seu,"SingleCellExperiment")) return(seu)
  if(!is(seu,"Seurat")) stop("not a Seurat object!")
  sce <- SingleCellExperiment(
    list(counts=GetAssayData(seu, assay="RNA", slot="counts")), 
    colData = seu[[]] )
  if(nrow(norm <- GetAssayData(seu, slot="scale.data"))>0){
    sce <- sce[row.names(norm),]
    logcounts(sce) <- norm
  }
  rowData(sce) <- seu@misc$rowData[row.names(sce),]
  if(length(VariableFeatures(seu)))
    metadata(sce)$VariableFeats <- VariableFeatures(seu)
  if(length(Reductions(seu))>0){
    reducedDims(sce) <- lapply( seu@reductions, 
                                FUN=function(x) x@cell.embeddings )
  }
  sce
}
