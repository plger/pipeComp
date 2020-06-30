def.filter <- function(x, minCounts=10){
  minCounts <- as.numeric(minCounts)
  x[filterByExpr(assay(x), model.matrix(~x$condition), min.count=minCounts),]
}

dea.edgeR <- function(x, mm, return.fit=FALSE){
  library(edgeR)
  dds <- calcNormFactors(DGEList(assay(x)))
  dds <- estimateDisp(dds, mm)
  fit <- glmFit(dds, mm)
  if(return.fit) return(fit)
  as.data.frame(topTags(glmLRT(fit, "condition"), Inf))
}

dea.edgeR.QLF <- function(x, mm){
  library(edgeR)
  dds <- calcNormFactors(DGEList(assay(x)))
  dds <- estimateDisp(dds, mm)
  res <- topTags(glmQLFTest(glmQLFit(dds, mm), "condition"), Inf)
  as.data.frame(res)
}

dea.voom <- function(x, mm){
  library(limma)
  library(edgeR)
  dds <- calcNormFactors(DGEList(assay(x)))
  v <- voom(dds, mm)
  fit <- eBayes(lmFit(v, mm), robust=TRUE)
  as.data.frame(topTable(fit, "condition", Inf))
}

dea.DESeq2 <- function(x, mm){
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(round(assay(x)), colData(x), design=mm)
  dds <- DESeq(dds)
  as.data.frame(results(dds, name = "condition"))
}

none <- function(x, ...)  return(x)

sva.vstsva <- function(x, k=NULL){
  library(sva)
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix( round(assay(x)), as.data.frame(colData(x)), 
                                 design=~condition )
  dds <- estimateSizeFactors(dds)
  en <- as.matrix(assay(vst(dds, blind = FALSE)))
  sv <- .svawrap(en, model.matrix(~x$condition), k, fn=sva::sva)
  if(!is.null(sv)) colData(x) <- cbind(colData(x), sv)
  x
}

sva.svaseq <- function(x, k=NULL){
  library(sva)
  sv <- .svawrap(.normcounts(assay(x)), model.matrix(~x$condition), k=k, 
                 fn=sva::svaseq)
  if(!is.null(sv)) colData(x) <- cbind(colData(x), sv)
  x
}

.svawrap <- function(x, mm, k=NULL, fn=sva::sva){
  if(is.null(k) || is.na(suppressWarnings(as.integer(k)))){
    if(is.null(k)){
      k <- "be"
    }else{
      k <- match.arg(k, c("be","leek"))
    }
    sv <- fn(x, mm, numSVmethod = k)
  }else{
    sv <- fn(x, mm, n.sv = as.integer(k))
  }
  if(sv$n.sv==0) return(NULL)
  colnames(sv$sv) <- paste0("SV", seq_len(ncol(sv$sv)))
  sv$sv
}

sva.RUVr <- function(x, k=NULL){
  library(RUVSeq)
  if(is.null(k) || is.na(suppressWarnings(as.integer(k)))) k <- 1
  mm <- model.matrix(~x$condition)
  resid <- residuals(dea.edgeR(x, mm, return.fit=TRUE), type="deviance")
  ruv <- RUVSeq::RUVr(round(.normcounts(assay(x))), row.names(x), resid, 
              k=as.integer(k))
  ruv <- ruv$W
  colnames(ruv) <- paste0("SV", seq_len(ncol(ruv)))
  colData(x) <- cbind(colData(x), ruv)
  x
}

sva.RUVs <- function(x, k=NULL){
  library(RUVSeq)
  if(is.null(k) || is.na(suppressWarnings(as.integer(k)))) k <- 1
  ri <- split(seq_len(ncol(x)), x$condition)
  ri <- lapply(ri, dr=max(sapply(ri, length)), FUN=function(x, dr){
    c(x,rep(-1,length(x)-dr))
  })
  ri <- do.call(rbind, ri)
  ruv <- RUVs(round(.normcounts(assay(x))), row.names(x), ri, k=as.integer(k))
  ruv <- ruv$W
  colnames(ruv) <- paste0("SV", seq_len(ncol(ruv)))
  colData(x) <- cbind(colData(x), ruv)
  x
}

.normcounts <- function(x){
  library(edgeR)
  if(is(x,"SummarizedExperiment")) x <- assay(x)
  d <- calcNormFactors(DGEList(x))
  return( edgeR::cpm(d, normalized.lib.sizes=TRUE) * 
            mean(d$samples$lib.size)/1000000 )
}
