
# ==============================================================================
# Author: Anthony Sonrel 
# Email:  anthony.sonrel@uzh.ch
# Created date: -
# ------------------------------------------------------------------------------
# Project: 
#   pip(e)Comp
#
# Description: 
#   Commands used to run the data imputations from the individual wrappers using
#   SCEs as input. 
#   It is better to run the following tools on a hpc: alra, dca, DrImpute, 
#   EmpiricalBose, enhance.
#   For the all-in-one wrapper, replace the individual wrappers with 
#   'DoDataImputation()' using file paths as input and specifying the method
#   with the 'method' arg. 
#
# ------------------------------------------------------------------------------
# TODO: 
#
# ==============================================================================

# individual wrappers
system.file("extdata", "impWrappers.R", 
            package = "pipeComp")
# all-in-one wrapper
system.file("extdata", "DoDataImputation.R", 
            package = "pipeComp")

#  SCIMPUTE --------------------------------------------------------------------
# Change the selected lines if modified k.  

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/scImpute/"   # <-----
dir.create(out_dir, showWarnings = FALSE)
for (i in input) {
  
  cat("Imputing file ", i, "\n")
  input <- readRDS(i)
  res <- wrapp.scimpute(input,
                        organism = "auto",
                        n_cores = 1L, 
                        kcluster_modif = 0)    # <-----
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}


# SAVERX  ----------------------------------------------------------------------
# WARNING: for some reason, may not work in for loop AND can not run multiple imputations in the same r session (open issue of the package). 
# ==> change manually the selected lines to run on all input files if this happens. 

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/SAVERX/"
dir.create(out_dir, showWarnings = FALSE)
for (i in input) {
  
  if (length(grep("Kumar|simMix2", i)) == 1) organism <- "Mouse" else organism <- "Human"
  cat("Imputing file ", i, "\n")
  input <- readRDS(i)
  res <- wrapp.saverx(input,
                      organism = organism, 
                      n_cores = 1L)   
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}


# DCA --------------------------------------------------------------------------

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/DCA/"
dir.create(out_dir, showWarnings = FALSE)
for (i in input) {
  
  cat("Imputing file ", i, "\n")
  input <- readRDS(i)
  res <- wrapp.dca(input, 
                   dca_path = "/home/asonrel/miniconda3/bin/dca",
                   n_cores = 1L)
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

# Empirical Bose ---------------------------------------------------------------
# HVGs selection on the selected line. 

input = list.files("datasets", full.names = T, pattern = "\\.csv")
out_dir <- "imputations/EmpiricalBose/"
dir.create(out_dir, showWarnings = FALSE)

for (i in input) {
  
  cat("Imputing file ", i, "\n")
  input <- readRDS(i)
  res <- wrapp.empiricalbose(input, 
                             python_path = "/home/asonrel/miniconda3/bin/python",
                             n_cores = 1L,
                             restr_to_hvgs = FALSE)  # <-------
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

# ALRA -------------------------------------------------------------------------

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/alra/"
dir.create(out_dir, showWarnings = FALSE)

for (i in input) {
  
  cat("Imputing file ", i, "\n")
  input <- readRDS(i)
  res <- wrapp.alra(input,
                    alra_norm = FALSE, # <-----
                    alra_path = "/home/asonrel/softwares/ALRA-master/alra.R" # <-----
  )
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

### ENHANCE ------------------------------------------------------------------

enhance_path = "/home/asonrel/softwares/enhance-R-master/enhance.R"

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds") 
out_dir <- "imputations/enhance/"
dir.create(out_dir, showWarnings = FALSE)

for (i in input) {
  
  cat("Imputing file ", i, "\n")
  input <- readRDS(i)
  res <- wrapp.enhance(input,
                       enhance_path = enhance_path) # <------
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

# DCA --------------------------------------------------------------------------

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds") 
out_dir <- "imputations/dca/"
dir.create(out_dir, showWarnings = FALSE)

for (i in input) {
  
  cat("Imputing file ", i, "\n")
  input <- readRDS(i)
  res <- wrapp.dca(input,
                   dca_path = "/home/asonrel/miniconda3/bin/dca",
                   n_cores = 1L)
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

# DRimpute ---------------------------------------------------------------------

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/DrImpute_noprocess/"
dir.create(out_dir, showWarnings = FALSE)

for (i in input) {
  
  cat("Imputing file ", i, "\n")
  input <- readRDS(i)
  res <- wrapp.drimpute(input, 
                        n_cores = 1L, 
                        DrImpute_prepross = FALSE) # <-------
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

