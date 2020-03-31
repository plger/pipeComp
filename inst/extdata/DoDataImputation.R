# ==============================================================================
# Author: Anthony Sonrel 
# Email:  anthony.sonrel@uzh.ch
# Created date: -
# ------------------------------------------------------------------------------
# Project: 
#   pipeComp
#
# Description: 
#   Wrapper of data imputation functions, optimized for .rds or .csv files.  
#   Arguments of the wrapper below the function start.
#   Designed to be run locally or on a server (use the correct python/ tool paths!)
#   For scImpute and DCA, the wrapper will create temporary folders for 
#   intermediate results. 
#
# ------------------------------------------------------------------------------
# TODO: 
#
# ==============================================================================


DoDataImputation <- function(count, 
                             method = c("SAVERX", "scImpute", "enhance", 
                                        "dca", "EmpiricalBose",  "DrImpute", 
                                        "alra"),
                             organism = c("auto", "Human", "Mouse"), 
                             out_dir = getwd(), 
                             n_cores = 3L, 
                             enhance_path = NULL, 
                             alra_path = NULL, 
                             alra_norm = FALSE,
                             DrImpute_prepross = TRUE, 
                             dca_path = "/home/asonrel/miniconda3/bin/dca", 
                             eb_out = "dirichlet",  
                             scImpute_kcluster = NULL, 
                             kcluster_modif = 0, 
                             restr_to_hvgs = FALSE,
                             python_path = '/usr/bin/python3'){
  # count: path to the .rds (SCE) or .csv file containing the dataset to impute on.
  # method: imputation package to use. 
  # organism: can keep 'auto' for the datasets we use in the pipeComp project. Else, please specify if "Human" or "Mouse".
  # out_dir: deprecated, not used for the moment
  # n_cores: number of cores to use for SAVERX, scImpute, DCA, EmpiricalBose. 
  # enhance_path: path to the enhance R script tool, to dl from https://github.com/yanailab/enhance.
  # alra_path: path to the 'alra.R' script, to dl from https://github.com/KlugerLab/ALRA.
  # alra_norm: logical, whether to normalize the data with log counts before imputing with alra (advised in the vignettes but not stated as required).
  # DrImpute_prepross: logical, filter the lowly expressed genes/ cells before running DrImpute (advised in the vignette but not stated as required) ?
  # dca_path: path to the installation of dca. Alternatively, can be set to "dca" if the appropriate virtual environment has been set with 'reticulate'. 
  # python_path: to specify if you installed 'sctransfer' on a particular python path
  # eb_out: EmpiricalBose method, which output to return: 1) "dirichlet" (default) similar to raw counts, 2) "log_p" similar to log(tmp). 
  # kcluster_modif : add or substract the estimated number of cluster by this parameter when using scimpute.
  # restr_to_hvgs : EmpiricalBose, logical, restrict the imputation to HVGs only ? 
  # python_path: path to python bin, passed to reticulate for SAVERX, EmpiricalBose
  # scImpute_kcluster: integer, number of kcluster expected in the data. Can be set. But by default, is is automatically inferred if a 'sce' if given with a 'phenoid' column in the colData or in the colnames of csv files (warning, suffix numbers are removed as it often use for cell id). 
  
  suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(purrr)
  })
  
  if (all(c("auto", "Human", "Mouse") == organism)) organism <- "auto"
  if (length(method) != 1) stop("Please provide 1 method.")
  
  # SAVERX   ===================================================================
  if (method == "SAVERX") {
    
    if (length(grep("\\.rds$|\\.txt$|\\.csv$|.rds$", count)) != 1) stop("SAVERX needs either a txt, csv or rds file.")
    if (organism == "auto") {
      if (length(grep("Kumar|simMix2", count)) == 1) organism <- "Mouse" else organism <- "Human"
    } 
    
    suppressPackageStartupMessages({
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
    
 
    if (length(grep("\\.rds$", count)) == 1) {
      
      library(SingleCellExperiment)
      data <- readRDS(count) 
      saveRDS(counts(data), file = "temp.rds")
      
      res <- saverx(input.file.name = "temp.rds", data.species = organism,
                    clearup.python.session = TRUE, ncores = n_cores)
      file.remove("temp.rds")
      
    } else {
      
      res <- saverx(input.file.name = count, data.species = organism,
                    clearup.python.session = TRUE, ncores = n_cores)
      
    }
    
    denois <- readRDS(res) %>% .$estimate
    # file.remove(res, recursive = TRUE)
    # unlink(gsub("\\/.*", "", res), recursive = TRUE)
    cat("If you want to run SAVERX a second time, please restart the R session (open issue of the package).\n")
    
  }
  
  # SCIMPUTE ===================================================================
  
  if (method == "scImpute") {
    
    suppressPackageStartupMessages({
      library("penalized")
      library("scImpute")
    })
    
    rand_suffix <- sample(1:1000000, 1)
    rand_dir <- paste0("temp_scImpute", rand_suffix, "/")
    dir.create(rand_dir)
    
    # file format
    if (length(grep("\\.rds", count)) == 1) {
      
      data <- readRDS(count) 
      if(!"phenoid" %in% colnames(colData(data))) stop("No 'phenoid' in the colData names. Please add it. It should contain the identity of the cells.")
      infile = "rds"
      if (length(scImpute_kcluster) == 0) kcluster <- length(unique(data$phenoid)) else kcluster <- scImpute_kcluster
      
      count <- paste0("temp_scImpute", rand_suffix, ".rds")
      saveRDS(counts(data), file = paste0(rand_dir, count))
      input_dir <- rand_dir
      
    } else if (length(grep("\\.csv", count)) == 1) {
      
      infile <- "csv"
      data <- read.csv(count, row.names = 1) 
      if (length(scImpute_kcluster) == 0) kcluster <- gsub("[0-9]*$", "", colnames(data)) %>% unique(.) %>% length(.) else kcluster <- scImpute_kcluster
      input_dir <- "./"
      
    }
    
    cat(kcluster, "clusters recognized in the input data.")
    
    if (kcluster_modif != 0) {
      
      kcluster <- kcluster + kcluster_modif
      cat("...modified to " , kcluster, "by 'kcluster_modif parameter.\n")
      
    }
    
    library(doParallel)
    scimpute(paste0(input_dir, count), ncores = n_cores, Kcluster = kcluster , 
             out_dir =  paste0(getwd(), "/", rand_dir), infile = infile, outfile = "rds")
    denois <- readRDS(paste0(rand_dir, "scimpute_count.rds"))
    
    # unlink(rand_dir, recursive = TRUE)
    
  }
  
  # ALRA     ===================================================================
  
  if (method == "alra") {
    
    library(rsvd)
    source(alra_path)
    
    if (length(grep("\\.rds", count)) == 1) {
      
      data <- readRDS(count) 
      data <- as.matrix(counts(data))
      
    } else if (length(grep("\\.csv", count)) == 1) {
      
      data <- read.csv(count, row.names = 1) 
      
    }
    
    data <- t(data)
    if (alra_norm) data <- normalize_data(log2(data + 1))
    k_choice <- choose_k(data)
    
    
    denois <- alra(data, k = k_choice$k)[[3]]
    rownames(denois) <- rownames(data)
    if (!all(dim(denois) == dim(data))) stop("Output dims differ from input dims. Please verify!")
    denois <- t(denois)
    
  }

  # DCA ========================================================================

  if (method == "dca") {
    
    if (!is(count, "character")) stop("'dca' needs raw count data in TSV/CSV format as input.")
    
    rand_suffix <- sample(1:1000000, 1)
    temp_file <- paste0("temp_DCA", rand_suffix, ".csv")
    temp_fold <- paste0("temp_DCA", rand_suffix)
    
    if (length(grep("\\.rds", count)) == 1) {
      
      temp <- readRDS(count)
      temp <- counts(temp)
      
      # DCA only accepts integers
      # --> floor the count if decimals
      if (all(colSums(temp) %% 1 == 0) == FALSE) temp <- floor(temp)
      write.csv(as.matrix(temp), file = temp_file, quote = FALSE)
      count <- temp_file
      
    } else {
      
      temp <- read.csv(count, row.names = 1)
      if (all(colSums(temp) %% 1 == 0) == FALSE) temp <- floor(temp)
      write.csv(temp, file = temp_file, quote = FALSE)
      
    }
    
    command <- paste(dca_path, temp_file, temp_fold, "--threads", n_cores)
    
    system(command = command)
    denois <- read.csv(paste0(temp_fold, "/mean.tsv"),sep = "\t", row.names = 1 )
    # unlink(temp_fold, recursive = TRUE)
    # file.remove(temp_file)
    
  }
  
  # EMPIRICALBOSE ==============================================================
  
  if (method == "EmpiricalBose") {
    
    suppressPackageStartupMessages({
      library(reticulate)
      library(scran)
      use_python(python_path, required = TRUE)
      library('EmpiricalBose')
      library('future')
      plan("multicore", workers = n_cores)
    })
    
    if (length(grep("\\.rds", count)) == 1) {
        
      data <- readRDS(count)
      logcounts(data) <- log2(counts(data) + 1)
      
    } 
    if (length(grep("\\.csv", count)) == 1) {
      
      data <- read.csv(count, row.names = 1)
      data <- SingleCellExperiment(list(counts = as.matrix(data), 
                                        logcounts = log1p(data)))

    } 
    
    # restrict to highly variable genes
    if (restr_to_hvgs) data <- restrict_to_hvgs(data)

    data_denois = cancel_noise(sce = data)
    
    denois <- assay(data_denois, eb_out)
    
    # reorder 
    if (all(dim(denois) == dim(counts(data)))) {
      denois <- denois[, match(colnames(data), colnames(denois))]
      denois <- denois[match(rownames(data), rownames(denois)), ]
    }

  }
  # ENHANCE ====================================================================
  
  if (method == "enhance") {
    
    cat("Using 'enhance' tool\n")
    source(enhance_path)
    if (length(grep("\\.rds", count)) == 1) {
      
      data <- readRDS(count) 
      data <- as.matrix(counts(data))
      
    } else if (length(grep("\\.csv", count)) == 1) {
      
      data <- read.csv(count, row.names = 1) 
      
    }
    
    denois <- try(enhance(data))
    if(is(denois, "try-error"))  denois <- enhance(data, k_nn = 2)

  }
  
  # DrImpute ===================================================================
  
  if (method == "DrImpute") {
    
    suppressPackageStartupMessages({
      library('DrImpute')
    })
    
    # preprocessign function from DrImpute 
    preprocess <- function(x, min.expressed.gene = 0, min.expressed.cell = 2, max.expressed.ratio = 1, normalize.by.size.effect = FALSE){
      
      if (is(x, 'SummarizedExperiment'))
        X <- assays(x)$count
      else if (is(x,'matrix'))
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
    } # end of preprocess
    
    if (length(grep("\\.rds", count)) == 1) {
      
      data <- readRDS(count)
      data <- counts(data)
      
    } 
    
    if (length(grep("\\.csv", count)) == 1) {
      
      data <- read.csv(count, row.names = 1)
      data <- as.matrix(data)
      
    }
    
    if (DrImpute_prepross == TRUE) X <- preprocess(data, min.expressed.gene = 0) else X <- data    
    Xlog <- log(X + 1)
    denois <- DrImpute(Xlog, mc.cores = n_cores)
    denois <- exp(denois) - 1
    
  }
  
  return(denois)
  
}










