#' dea_pipeline
#' 
#' The `PipelineDefinition` for bulk RNAseq differential expression analysis 
#' (DEA).
#' 
#' @return A `PipelineDefinition` object to be used with `runPipeline`.
#' 
#' @export
#' @examples
#' pip <- dea_pipeline()
#' pip
dea_pipeline <- function(){
  # description for each step
  desc <- list( 
    filtering=
paste("Takes a SE object, passes it through the function `filt` (with param",
" `minCount`), and outputs a filtered SE object."),
    sva=
paste("Takes a SE object, passes it through the function `sva.method`, and",
      "returns a SE object."),
    dea=
paste("Takes a SE object, passes it through the function `dea.method`, and",
      "returns a DEA data.frame.")
  )
  
  # functions list
  f <- list(
    filtering=function(x, filt, minCount=10) get(filt)(x, minCount=minCount),
    sva=function(x, sva.method, k=1){ 
      get(sva.method)(x, k=k)
    },
    dea=function(x, dea.method){
      mm <- pipeComp:::.getMM(x)
      x2 <- pipeComp:::.homogenizeDEA(get(dea.method)(x,mm))
      metadata(x2)$truth <- metadata(x)$truth
      metadata(x2)$mm <- mm
      x2
    }
  )
  
  agg <- eva <- lapply(f, FUN=function(x) NULL)
  eva$dea <- evaluateDEA
  agg$dea <- aggregateDEAeval

  # default arguments
  def <- list(minCount=10, k=1)

  # initiation function
  initf <- function(x){
    if(is.character(x) && length(x)==1) x <- readRDS(x)
    metadata(x)$truth <- rowData(x)[,c("expected.beta","isDE")]
    cn <- grep("^SV[[:digit:]]+$", colnames(colData(x)), invert=TRUE)
    colData(x) <- colData(x)[,cn,drop=FALSE]
    x
  }

  PipelineDefinition(functions=f, descriptions=desc, evaluation=eva,
                     aggregation=agg, initiation=initf, 
                     defaultArguments=def, verbose=FALSE)
}

#' @importFrom S4Vectors DataFrame
.homogenizeDEA <- function(x, g=NULL){
  colnames(x)[which(colnames(x) %in% c("FDR","padj","adj.P.Val"))] <- "FDR"
  colnames(x)[which(colnames(x) %in% c("P.Value","pvalue","PValue"))] <- "PValue"
  colnames(x)[which(colnames(x) %in% c("log2FoldChange","logFC"))] <- "logFC"
  if(is.null(g)) g <- row.names(x)
  x <- x[g,c("logFC","PValue","FDR")]
  row.names(x) <- g
  if(!is(x,"DataFrame")) x <- DataFrame(x)
  x
}

#' @importFrom stats model.matrix
.getMM <- function(se){
  v <- c( grep("^SV[[:digit:]]+$", colnames(colData(se)), value=TRUE), 
          "condition" )
  f <- as.formula(paste0("~",paste(v,collapse="+")))
  model.matrix(f, data=as.data.frame(colData(se)))
}
