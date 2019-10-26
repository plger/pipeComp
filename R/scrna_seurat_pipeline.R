#' scrna_seurat_pipeline
#' 
#' The `PipelineDefinition` for the default scRNAseq clustering pipeline, with 
#' steps for doublet identification, filtering, normalization, feature selection,
#' dimensionality reduction, and clustering.
#' Alternative arguments should be character, numeric or logical vectors of
#' length 1 (e.g. the function name for a method, the number of dimensions, etc).
#' The default pipeline has the following steps and arguments:
#' * doublet: `doubletmethod` (name of the doublet removal function)
#' * filtering: `filt` (name of the filtering function, or filter string)
#' * normalization: `norm` (name of the normalization function)
#' * selection: `sel` (name of the selection function, or variable of rowData on which to select) and 
#'  `selnb` (number of features to select)
#' * dimreduction: `dr` (name of the dimensionality reduction function) and `maxdim` (maximum number 
#' of components to compute)
#' * clustering: `clustmethod` (name of the clustering function), `dims` (number of dimensions to use),
#'  `k` (number of nearest neighbors to use, if applicable), `steps` (number of steps in the random 
#'  walk, if applicable), `resolution` (resolution, if applicable), `min.size` (minimum cluster size,
#'  if applicable)
#' 
#' @param saveDimRed Logical; whether to save the dimensionality reduction for 
#' each analysis (default FALSE)
#' 
#' @importFrom mclust adjustedRandIndex
#' @export
scrna_seurat_pipeline <- function(saveDimRed=FALSE){
  # description for each step
  desc <- list( 
    doublet=
"Takes a SCE object with the `phenoid` colData column, passes it through the 
function `doubletmethod`, and outputs a filtered SCE.",
    filtering=
"Takes a SCE object, passes it through the function `filt`, and outputs a 
filtered Seurat object.",
    normalization=
"Passes the object through function `norm` to return the object with the 
normalized and scale data slots filled.",
    selection=
"Returns a seurat object with the VariableFeatures filled with `selnb` features 
using the function `sel`.",
    dimreduction=
"Returns a seurat object with the PCA reduction with up to `maxdim` components
using the `dr` function.",
    clustering=
"Uses function `clustmethod` to return a named vector of cell clusters." )
  
  # we prepare the functions for each step
  if(saveDimRed){
    DRfun <- function(x, dr, maxdim){ 
      x <- get(dr)(x, dims=maxdim)
      list( x=x, intermediate_return=list( cell.embeddings=x@reductions$pca@cell.embeddings,
                                           evaluation=evaluateDimRed(x) ) )
    }
  }else{
    DRfun <- function(x, dr, maxdim){ 
      get(dr)(x, dims=maxdim)
    }
  }
  f <- list(
        doublet=function(x, doubletmethod){ 
          x2 <- pipeComp:::.runf(doubletmethod, x)
          list( x=x2, intermediate_return=pipeComp:::.compileExcludedCells(x,x2) )
        },
        filtering=function(x, filt){
          x2 <- pipeComp::.runf(filt, x, alt=applyFilterString)
          list( x=x2, intermediate_return=pipeComp:::.compileExcludedCells(x,x2) )
        },
        normalization=function(x, norm){ 
          get(norm)(x)
        },
        selection=function(x, sel, selnb){ 
          x <- .runf(sel, x, n=selnb, alt=applySelString)
          list( x=x, intermediate_return=Seurat::VariableFeatures(x) )
        },
        dimreduction=DRfun,
        # dimensionality=function(x, dims){
        #   if(!is.na(suppressWarnings(as.numeric(dims)))){
        #     dims <- as.integer(dims)
        #   }else{
        #     dims <- getDimensionality(x, dims)
        #   }
        #   x@reductions$pca@cell.embeddings <- x@reductions$pca@cell.embeddings[,seq_len(dims)]
        # },
        clustering=function(x, clustmethod, dims, k, steps, resolution, min.size){
          tl <- x$phenoid
          dims <- .parseDims(dims)
          x <- get(clustmethod)(x, dims=dims, resolution=resolution, k=k, steps=steps, min.size=min.size)
          e <- match_evaluate_multiple(x, tl)
          unmatched <- length(x)-sum(e$n_cells_matched)
          e <- c( unlist(e), unmatched.cells=unmatched, 
                  ARI=mclust::adjustedRandIndex(as.character(x),tl) )
          list( x=x, intermediate_return=e )
        }
  )
  
  eva <- list( doublet=NULL, filtering=NULL, normalization=evaluateNorm,
               selection=NULL, dimreduction=evaluateDimRed , clustering=NULL )
  # functions to aggregate the intermediate_return of the different steps
  agg <- list( doublet=.aggregateExcludedCells,
               filtering=.aggregateExcludedCells,
               normalization=NULL,
               selection=NULL,
               dimreduction=.aggregateDR,
               clustering=.aggregateClusterEvaluation )
  
  PipelineDefinition(functions=f, descriptions=desc, evaluation=eva,
                            aggregation=agg, verbose=FALSE)
}

.parseDims <- function(x, dims){
  if(!is.na(suppressWarnings(as.numeric(dims)))) return(as.integer(dims))
  getDimensionality(x, dims)
}