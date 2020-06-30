# ==============================================================================
# Author: Anthony Sonrel 
# Email:  anthony.sonrel@uzh.ch
# Created date: Wed Nov 13 2019
# ------------------------------------------------------------------------------
# Project: pipeComp
#
# Description: individual wrappers of each imputation function used in the 
#   pipeComp project. For all of these function, both the input and the output
#   is an SCE object. 
#
# ------------------------------------------------------------------------------
# TODO: 
# ==============================================================================


#' wrapp.saverx
#' 
#' Wrapper for SAVERX imputation. 
#' 
#' @param sce An SCE object to perform the imputation on. 
#' @param organism "Mouse" or "Human". Origin of the dataset.
#' @param python_path Path to the python binary that has an installtion of SAVERX. 
#' @param n_cores Number of cores to use. 
#' 
#' @return An SCE object with the imputed data in the 'counts' assay. 
#' @export

wrapp.saverx <- function(sce, 
                         organism,
                         python_path = '/usr/bin/python3', 
                         n_cores = 1L
                         ){
  
  
  suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(purrr)
    library(SAVERX)
    library(reticulate)
    reticulate::use_python(python = python_path, required = TRUE)
    library(tensorflow)
    library(keras)
    config <- tf$ConfigProto(intra_op_parallelism_threads = n_cores,
                             inter_op_parallelism_threads = n_cores)
    session = tf$Session(config = config)
    k_set_session(session)
  })
  
  saveRDS(counts(sce), file = "temp.rds")
  
  res <- saverx(input.file.name = "temp.rds", data.species = organism,
                clearup.python.session = TRUE, ncores = n_cores)
  file.remove("temp.rds")

  denois <- readRDS(res) %>% .$estimate
  file.remove(res)
  
  .out_check(sce, denois)
  denois <- .prep(sce, denois, out = "output")
  counts(sce) <- denois
  # unlink(gsub("\\/.*", "", res), recursive = TRUE)
  cat("If you want to run SAVERX a second time, please restart the R session (open issue of the package).\n")
  
  return(sce)

}


#' wrapp.scimpute
#' 
#' Wrapper for scImpute imputation. 
#' 
#' @param sce An SCE object to perform the imputation on. 
#' @param n_cores Number of cores to use. 
#' @param scImpute_kcluster: integer, k clusters expected in the data. Can be set.By default, this number is automatically detected if a `sce` if given with a 'phenoid' column in the colData. 
#' @param kcluster_modif Integer, number of clusters to add or substract to the estimated number of clusters.

#' 
#' @return An SCE object with the imputed data in the 'counts' assay.
#' @export
 
wrapp.scimpute <- function(sce, 
                           scImpute_kcluster = NULL, 
                           kcluster_modif = 0, 
                           n_cores = 1L
                           ) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(purrr)
    library(penalized)
    library(scImpute)
    library(doParallel)
  })
  
  rand_suffix <- sample(1:1000000, 1)
  rand_dir <- paste0("temp_scImpute", rand_suffix, "/")
  dir.create(rand_dir)
  
  if(!"phenoid" %in% colnames(colData(sce))) stop("No 'phenoid' in the colData names. Please add it. It should contain the identity of the cells.")
  infile = "rds"
  if (length(scImpute_kcluster) == 0) kcluster <- length(unique(sce$phenoid)) else kcluster <- scImpute_kcluster
  
  count_path <- paste0("temp_scImpute", rand_suffix, ".rds")
  saveRDS(counts(sce), file = paste0(rand_dir, count_path))
  input_dir <- rand_dir
  
  cat(kcluster, "clusters recognized in the input data.")
  
  if (kcluster_modif != 0) {
    
    kcluster <- kcluster + kcluster_modif
    cat("...modified to " , kcluster, "by 'kcluster_modif parameter.\n")
    
  }
  
  scimpute(paste0(input_dir, count_path), ncores = n_cores, Kcluster = kcluster , 
           out_dir =  paste0(getwd(), "/", rand_dir), infile = infile, outfile = "rds")
  denois <- readRDS(paste0(rand_dir, "scimpute_count.rds"))
  
  unlink(rand_dir)
  .out_check(sce, denois)
  denois <- .prep(sce, denois, out = "output")
  counts(sce) <- denois
  
  return(sce)
}


#' wrapp.alra
#' 
#' Wrapper for ALRA imputation. 
#' 
#' @param sce An SCE object to perform the imputation on. 
#' @param alra_norm logical, whether to normalize the data with log counts before imputing with alra (advised in the vignettes but not stated as required).
#' @param alra_path path to the `alra.R` script that contains the required functions for the imputation. Can be dl from https://github.com/KlugerLab/ALRA.
#' 
#' @return An SCE object with the imputed data in the 'counts' assay. 
#' @export

wrapp.alra <- function(sce,                              
                       alra_norm = FALSE,
                       alra_path){
  
  suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(purrr)
    library(rsvd)
    source(alra_path)
  })
  
  data <- as.matrix(counts(sce))
  data <- t(data)
  if (alra_norm) data <- normalize_data(log2(data + 1))
  k_choice <- choose_k(data)

  denois <- alra(data, k = k_choice$k)[[3]]
  rownames(denois) <- rownames(data)
  denois <- t(denois)
  
  .out_check(sce, denois)
  denois <- .prep(sce, denois, out = "output")
  
  counts(sce) <- denois
  
  return(sce)
  
}


#' wrapp.dca
#' 
#' Wrapper for DCA imputation. 
#' 
#' @param sce An SCE object to perform the imputation on. 
#' @param dca_path path to the python installation of DCA. Alternatively, can be set to "dca" if the appropriate virtual environment has been set with `reticulate`` on top of the function. 
#' @param n_cores Number of cores to use. 
#' 
#' @return An SCE object with the imputed data in the 'counts' assay. 
#' @export

wrapp.dca <- function(sce, 
                      dca_path = "/usr/bin/python3/dca",                              
                      n_cores = 1L
                      ){
  
  suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(purrr)
  })
  
  rand_suffix <- sample(1:1000000, 1)
  temp_file <- paste0("temp_DCA", rand_suffix, ".csv")
  temp_fold <- paste0("temp_DCA", rand_suffix)
  
  temp <- counts(sce)
  
  # DCA only accepts integers
  # --> floor the count if decimals
  if (all(colSums(temp) %% 1 == 0) == FALSE) temp <- floor(temp)
  write.csv(as.matrix(temp), file = temp_file, quote = FALSE)

  command <- paste(dca_path, temp_file, temp_fold, "--threads", n_cores)
  
  system(command = command)
  denois <- read.csv(paste0(temp_fold, "/mean.tsv"),sep = "\t", row.names = 1 )
  unlink(temp_fold, recursive = TRUE)
  file.remove(temp_file)
  
  sce <- .prep(sce, denois, out = "input")
  counts(sce) <- denois
  
  return(sce)
}

#' wrapp.empiricalbose
#' 
#' Wrapper for EmpiricalBose imputation. 
#' 
#' @param sce An SCE object to perform the imputation on. 
#' @param python_path Path to a python3 binary that has `numpy` installed.
#' @param n_cores Number of cores to use. 
#' @param restr_to_hvgs Logical; restrict the imputation to HVGs only ? If `TRUE`, the output will only contain these genes.
#' @param eb_out EmpiricalBose method: which output to return. Can be "dirichlet" (default) similar to raw counts or "log_p" similar to log(tmp). 
#' 
#' @return An SCE object with the imputed data in the 'counts' assay. 
#' @export

wrapp.empiricalbose <- function(sce, 
                                python_path = '/usr/bin/python3',
                                restr_to_hvgs = FALSE,
                                eb_out = "dirichlet", 
                                n_cores = 1L) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(purrr)
    library(reticulate)
    library(scran)
    use_python(python_path, required = TRUE)
    library('EmpiricalBose')
    library('future')
    plan("multicore", workers = n_cores)
  })
  
  logcounts(sce) <- log2(counts(sce) + 1)
  
  # restrict to highly variable genes
  if (restr_to_hvgs) sce <- restrict_to_hvgs(sce)
  
  data_denois = cancel_noise(sce = sce)
  
  denois <- assay(data_denois, eb_out)
  if (!restr_to_hvgs) {
    .out_check(sce, denois)
    denois <- .prep(sce, denois, out = "output")
  } else {
    sce <- .prep(sce, denois, out = "input")
  }
  
  
  counts(sce) <- denois
    
  return(sce)
  
}

#' wrapp.enhance
#' 
#' Wrapper for ENHANCE imputation. 
#' 
#' @param sce An SCE object to perform the imputation on. 
#' @param enhance_path Path to the R script containing the ENHANCE functions. Can be dl from https://github.com/yanailab/enhance.
#' 
#' @return An SCE object with the imputed data in the 'counts' assay. 
#' @export

wrapp.enhance <- function(sce, 
                          enhance_path) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(purrr)
    source(enhance_path)
  })
  
  data <- as.matrix(counts(sce))
 
  denois <- try(enhance(data))
  if (is(denois, "try-error"))  denois <- enhance(data, k_nn = 2)
  
  .out_check(sce, denois)
  denois <- .prep(sce, denois, out = "output")
  
  counts(sce) <- denois
  
  return(sce)
  
}


#' wrapp.drimpute
#' 
#' Wrapper for ENHANCE imputation. 
#' 
#' @param sce An SCE object to perform the imputation on. 
#' @param n_cores Number of cores to use. 
#' @param DrImpute_prepross Logical, filter the lowly expressed genes/ cells before running DrImpute (advised in the vignette but not stated as required) ?
#' 
#' @return An SCE object with the imputed data in the 'counts' assay. 
#' @export

wrapp.drimpute <- function(sce, 
                           DrImpute_prepross = TRUE, 
                           n_cores = 1L) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(purrr)
    library(DrImpute)
  })
  
  data <- counts(sce)
  
  if (DrImpute_prepross == TRUE) X <- .drimp_preprocess(data, min.expressed.gene = 0) else X <- data    
  Xlog <- log(X + 1)
  denois <- DrImpute(Xlog, mc.cores = n_cores)
  denois <- exp(denois) - 1
  
  sce <- .prep(sce, denois, out = "input")
  counts(sce) <- denois
  
  return(sce)
  
}



#' imp.scVI
#'
#' A function calling a python wrapper (`scVI.py`) around `scVI` imputation, adapted from the the 'Basic usage' Jupyter notebook (https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/basic_tutorial.ipynb). Note that the function will create a temporary csv file for the intermediate storage of the input count matrix, needed by `scVI`.  
#' 
#' 
#' @param x A SCE object.
#' @param py_script Location of the python script
#' @param py_path Optional. If scVI was installed in a specific python bin, pass here the path to it. 
#' @param train_size Size of training set. Default to 0.8 but tutorial however recommends to use 1. 
#' @param n_cores N. cores
#' 
#' 
#' @return An object of the same class as `x` with updated slots.

imp.scVI <- function(x, py_script = system.file("extdata", "scVI.py", package="pipeComp"), py_path = NULL, n_cores = 1L, train_size = 1) {
  n_cores <- as.integer(n_cores)
  suppressPackageStartupMessages(library(reticulate))
  if (length(py_path)>0) use_python(py_path ,required=TRUE)
  trysource <- try(source_python(py_script))
  if (class(trysource) == "try-error") stop("Cannot source the python wrapper.") 
  tfile <- tempfile(fileext=".csv", tmpdir = ".")
  write.csv(counts(x), tfile)
  out <- scVI_imput(csv_file = tfile, csv_path = ".", n_cores = n_cores,
                    train_size = train_size)
  val <- t(out[[1]])
  gnames <- as.character(out[[2]])
  file.remove(tfile)
  dimnames(val) <- list(gnames, colnames(counts(x)))
  x <- x[gnames, ]
  counts(x) <- as.matrix(val)
}



# MISC. ========================================================================

.out_check <- function(input, output) {
  # Check that input dimensions = output dimensions in some tools
  # else stop
  if (!all(dim(input) == dim(output))) {
    stop("Output dimension differ from input dimension.\nPlease verify. Dimensions from input and output should be the same")
  }
}


.prep <- function(input, output, out) {
  # Reordering of output rows/ columns
  if (out == "input") {
    input <- input[rownames(input) %in% rownames(output), 
                   colnames(input) %in% colnames(output)]
    input <- input[match(rownames(output), rownames(input)), 
                   match(colnames(output), colnames(input))]
    return(input)
  } else {
    output <- output[match(rownames(input), rownames(output)), 
                     match(colnames(input), colnames(output))]
    return(output)
  }
  
}

.drimp_preprocess <- function(x, min.expressed.gene = 0, min.expressed.cell = 2, 
                              max.expressed.ratio = 1, 
                              normalize.by.size.effect = FALSE){
  # preprocessign function from DrImpute 
  if (is(x, 'SummarizedExperiment'))
    X <- assays(x)$count
  else if (is(x, 'matrix'))
    X <- x
  else if (is(x, 'sparseMatrix'))
    X <- x
  else
    stop(sprintf('unknown class(x): %s', class(x)))
  
  M <- ncol(X)
  N <- nrow(X)
  m <- Matrix::colSums(X > 1) >= min.expressed.gene	# cells that have at least min.expressed.gene expreseed genes
  n <- Matrix::rowSums(X > 1) <= max.expressed.ratio * M & Matrix::rowSums(X > 1) >= min.expressed.cell	# genes that are detected in at least min.expressed.cell or at most max.expressed.ratio cells
  if (normalize.by.size.effect){
    sf <- apply((X[n, m] + 1) / exp(Matrix::rowMeans(log(X[n, m] + 1))), 2, median)
    X <- t(t(X[n, m]) / sf)
  }else
    X <- X[n, m]
  
  if (is(x, 'SummarizedExperiment')){
    x <- x[n, m]
    assays(x)$count <- X
  }else if (is(x, 'matrix')){
    x <- as.matrix(X)
  }else if (is(x, 'sparseMatrix')){
    x <- X
  }
  
  x
}



