
# ==============================================================================
# Author: Anthony Sonrel 
# Email:  anthony.sonrel@uzh.ch
# Created date: -
# ------------------------------------------------------------------------------
# Project: 
#   pip(e)Comp
#
# Description: 
#   Commands used to run the data imputations from 'DoDataImputation.R' wrapper. 
#   It is better to run the following tools on a hpc: alra, dca, DrImpute, 
#   EmpiricalBose, enhance.
#
# ------------------------------------------------------------------------------
# TODO: 
#
# ==============================================================================

source("DoDataImputation.R")

#  SCIMPUTE ====================================================================
# Change the selected lines if modified k.  

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/scImpute_plus5/"   # <-----
dir.create(out_dir, showWarnings = FALSE)
for(i in input) {
  
  cat("Imputing file ", i, "\n")
  res <- DoDataImputation(count = i, method = "scImpute", organism = "auto",
                          n_cores = 4, kcluster_modif = 5)    # <-----
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

# SAVERX ======================================================================
# WARNING: for some reason, does not work in for loop AND can not run multiple imputations in the same r session (open issue of the package). 
# ==> change manually the selected lines to run on all input files. 

input = list.files("filtCountMatrices", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/SAVERX/"
dir.create(out_dir, showWarnings = FALSE)

res <- DoDataImputation(count = input[1], method = "SAVERX", organism = "auto")   # <-----
out_file <- gsub(".*\\/", "", input[1]) %>% gsub("\\.csv", ".rds", .)    # <------
saveRDS(res, file = paste0(out_dir, out_file))



# DCA ==========================================================================

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/DCA/"
dir.create(out_dir, showWarnings = FALSE)
for(i in input[-c(1:2)]) {
  
  cat("Imputing file ", i, "\n")
  res <- DoDataImputation(count = i, method = "dca", organism = "auto",
                          n_cores = 6L)
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

# BayesBos =====================================================================
# HVGs selection on the selected line. 

input = list.files("filtCountMatrices", full.names = T, pattern = "\\.csv")
out_dir <- "imputations/EmpiricalBose/"
dir.create(out_dir, showWarnings = FALSE)

for(i in input) {
  
  cat("Imputing file ", i, "\n")
  res <- DoDataImputation(count = i, method = "EmpiricalBose",
                          python_path = "/usr/bin/python3",
                          n_cores = 2L,
                          restr_to_hvgs = FALSE)  # <-------
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}

### ALRA -----------------------------------------------------------------------

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/alra/"
dir.create(out_dir, showWarnings = FALSE)

for (i in input) {
  
  cat("Imputing file ", i, "\n")
  res <- DoDataImputation(count = i, method = "alra",
                          organism = "auto",
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
  res <- DoDataImputation(count = i, method = "enhance", 
                          organism = "auto", 
                          enhance_path = enhance_path) # <------
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}


#### DCA ---------------------------------------------------------------------

input = list.files("filtCountMatrices", full.names = T, pattern = "\\.csv|\\.rds") 
out_dir <- "imputations/dca/"
dir.create(out_dir, showWarnings = FALSE)

for (i in input) {
  
  cat("Imputing file ", i, "\n")
  res <- DoDataImputation(count = i, method = "dca",
                          organism = "auto",
                          n_cores = 6L)
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}


### DRimpute ---------------------------------------------------------------------

input = list.files("datasets", full.names = T, pattern = "\\.csv|\\.rds")
out_dir <- "imputations/DrImpute_noprocess/"
dir.create(out_dir, showWarnings = FALSE)

for (i in input) {
  
  cat("Imputing file ", i, "\n")
  res <- DoDataImputation(count = i, 
                          method = "DrImpute",
                          organism = "auto",
                          n_cores = 3L, 
                          DrImpute_prepross = FALSE) # <-------
  
  out_file <- gsub(".*\\/", "", i) %>% gsub("\\.csv$", ".rds", .)
  saveRDS(res, file = paste0(out_dir, out_file))
  
}


