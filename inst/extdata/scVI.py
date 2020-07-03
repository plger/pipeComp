
from scvi.dataset import LoomDataset, CsvDataset, Dataset10X, DownloadableAnnDataset
import urllib.request
import os
from scvi.dataset import BrainLargeDataset, CortexDataset, PbmcDataset, RetinaDataset, HematoDataset, CbmcDataset, BrainSmallDataset, SmfishDataset
import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scvi.models import VAE, AutoZIVAE, SCANVI, VAEC
from scvi.inference import UnsupervisedTrainer
from scvi import set_seed
import torch

def scVI_imput(csv_file, csv_path, vae_model = VAE, train_size = 1.0, n_labels = 0, seed = 1234, n_cores = 1, lr = 1e-3, use_cuda = False): 
  set_seed(seed)
  dat = CsvDataset(csv_file, 
                  save_path=csv_path, 
                  new_n_genes=None) 
  dat.subsample_genes(1000, mode="variance")
  # Based on recommendations in basic_tutorial.ipynb             
  n_epochs = 400 if (len(dat) < 10000) else 200 
  # trainer and model 
  vae = vae_model(dat.nb_genes, n_labels = n_labels)
  trainer = UnsupervisedTrainer(
    vae,
    dat,
    train_size=train_size, # default to 0.8, documentation recommends 1
    use_cuda=use_cuda    
    )
  # limit cpu usage
  torch.set_num_threads(n_cores) 
  trainer.train(n_epochs=n_epochs, lr=lr)
  full = trainer.create_posterior(trainer.model, dat, indices=np.arange(len(dat)))
  # Updating the "minibatch" size after training is useful in low memory configurations
  full = full.update({"batch_size":32})
  imp_val = full.sequential().imputation()
  return[imp_val, dat.gene_names]
  
def scVI_norm(csv_file, csv_path, vae_model = VAE, train_size = 1.0, n_labels = 0, seed = 1234, n_cores = 1, lr = 1e-3, use_cuda = False): 
  set_seed(seed)
  dat = CsvDataset(csv_file, 
                  save_path=csv_path, 
                   new_n_genes=None) 
  dat.subsample_genes(1000, mode="variance")
  # Based on recommendations in basic_tutorial.ipynb    
  n_epochs = 400 if (len(dat) < 10000) else 200 
  # trainer and model 
  vae = vae_model(dat.nb_genes, n_labels = n_labels)
  trainer = UnsupervisedTrainer(
    vae,
    dat,
    train_size=train_size, # default to 0.8, documentation recommends 1
    use_cuda=use_cuda    
    )
  # limit cpu usage
  torch.set_num_threads(n_cores) 
  trainer.train(n_epochs=n_epochs, lr=lr)
  full = trainer.create_posterior(trainer.model, dat, indices=np.arange(len(dat)))
  # Updating the "minibatch" size after training is useful in low memory configurations
  normalized_values = full.sequential().get_sample_scale()
  return[normalized_values, dat.gene_names]


def scVI_latent(csv_file, csv_path, vae_model = VAE, train_size = 1.0, n_labels = 0, seed = 1234, n_cores = 1, lr = 1e-3, use_cuda = False): 
  set_seed(seed)
  dat = CsvDataset(csv_file, 
                  save_path=csv_path, 
                   new_n_genes=None) 
  # Based on recommendations in basic_tutorial.ipynb    
  n_epochs = 400 if (len(dat) < 10000) else 200 
  # trainer and model 
  vae = vae_model(dat.nb_genes, n_labels = n_labels)
  trainer = UnsupervisedTrainer(
    vae,
    dat,
    train_size=train_size, # default to 0.8, documentation recommends 1
    use_cuda=use_cuda    
    )
  # limit cpu usage
  torch.set_num_threads(n_cores) 
  trainer.train(n_epochs=n_epochs, lr=lr)
  full = trainer.create_posterior(trainer.model, dat, indices=np.arange(len(dat)))
  # Updating the "minibatch" size after training is useful in low memory configurations
  Z_hat = full.sequential().get_latent()[0]
  adata = anndata.AnnData(dat.X)
  for i, z in enumerate(Z_hat.T):
      adata.obs[f'Z_{i}'] = z
  # reordering for convenience and correspondance with PCA's ordering
  cellLoads = adata.obs.reindex(adata.obs.std().sort_values().index, axis = 1)
  return(cellLoads)



from scvi.models import LDVAE
import anndata

def scVI_ld(csv_file, csv_path, ndims, vae_model = VAE, n_labels = 0, n_cores=1, seed= 1234, lr = 1e-3, use_cuda = False): 
  set_seed(seed)
  dat = CsvDataset(csv_file, 
                   save_path=csv_path, 
                   new_n_genes=None) 
  # Based on recommendations in linear_decoder.ipynb
  n_epochs = 250
  # trainer and model 
  ldvae = LDVAE(
        dat.nb_genes,
        n_batch = dat.n_batches,
        n_latent = ndims, 
        n_labels = n_labels
        )
  trainerLD = UnsupervisedTrainer(ldvae, dat, use_cuda=use_cuda)
  # limit cpu usage
  torch.set_num_threads(n_cores) 
  trainerLD.train(n_epochs=n_epochs, lr=lr)
  # extract mean value for the ld
  full = trainerLD.create_posterior(trainerLD.model, dat, indices=np.arange(len(dat)))
  Z_hat = full.sequential().get_latent()[0]
  adata = anndata.AnnData(dat.X)
  for i, z in enumerate(Z_hat.T):
      adata.obs[f'Z_{i}'] = z
  # reordering for convenience and correspondance with PCA's ordering
  cellLoads = adata.obs.reindex(adata.obs.std().sort_values().index, axis = 1)
  return(cellLoads)


