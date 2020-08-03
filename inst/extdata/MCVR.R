# Most of this code is from the R implementation of molecular cross-validation
# by Matthew Young ( https://github.com/constantAmateur/MCVR/blob/master/code.R 
#, commit c4caf93), with few modifications

library(Matrix)
library(Seurat)

MCVR <- function(se, verbose=FALSE, tFracs=1, nSplits=2, BPPARAM=NULL, ...){
  if(!is(se, "Seurat") || length(VariableFeatures(se))==0)
    stop("`se` should be a Seurat object with computed variable features.")
  se <- SetAssayData(se, slot="counts", new.data=floor(GetAssayData(se, slot="counts")))
  if(is.null(Misc(se, "norm.function"))){
    if(verbose) warning("No normalization function specified, will use standard Seurat normalization")
    nf <- function(x,...){
      ScaleData(NormalizeData(x,verbose=FALSE),verbose=FALSE)
    }
  }else{
    nf <- Misc(se, "norm.function")
  }
  if(verbose){
    if(is.null(Misc(se, "pca.function")))
      warning("No PCA function specified, will use irlba-approximated PCA")
    if(is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar=TRUE)
  }else{
    if(is.null(BPPARAM)) BPPARAM <- SerialParam()
  }
  mcv <- suppressWarnings({
    molecularCrossValidation(se, nf, Misc(se, "pca.function"), tFracs=tFracs, 
                             BPPARAM=BPPARAM, nSplits=nSplits, verbose=verbose)
  })
  if(verbose) plotMCV(mcv)
  k <- max(mcv$minimas, na.rm=TRUE)
  if(is.na(k) || is.infinite(k)) k <- ncol(Embeddings(se))
  if(k<2) k <- 2
  k
}

#' Create test and training data splits.
#'
#' Splits sample of full molecular pool into two samples from larger molecular pool (for test and training), accounting for the overlap between samples (p).
#'
#' @param src The input source matrix.
#' @param alpha Ratio of training to test fraction.
#' @param p Fraction of total molecules sampled.
#' @param tFrac Titration fraction.
#' @return Two sparse matrices, one for the test, one for the training.
partitionData = function(src,alpha,p,tFrac){
  #Do we need to titrate?  If so do that first.
  src = as(src,'dgTMatrix')
  x = src@x
  i = src@i
  j = src@j
  if(tFrac<1)
    x = rbinom(length(x),x,tFrac)
  #Adjust p by titration fraction
  p = p *tFrac
  #Now work out p' and p''
  p1 = (1+alpha-sqrt((1+alpha)**2-4*alpha*p))/(2*alpha)
  p2 = p1*alpha
  #Construct the training sample
  trx = rbinom(length(x),x,p1/p)
  #Construct the test sample, including an overlap adjustment.
  tsx = x-trx + rbinom(length(trx),trx,p2)
  #Reformat the training data into the matrix
  w = which(trx>0)
  train = sparseMatrix(x=trx[w],i=i[w]+1,j=j[w]+1,dims=dim(src),dimnames=dimnames(src))
  #Reformat the test data i nto the matrix
  w = which(tsx>0)
  tst = sparseMatrix(x=tsx[w],i=i[w]+1,j=j[w]+1,dims=dim(src),dimnames=dimnames(src))
  return(list(tst,train))
}


#' Finds the optimal number of PCs to use
#'
#' Uses molecular cross validation method (https://www.biorxiv.org/content/10.1101/786269v1), which must be applied to raw count data, to determine the optimal number of PCs to use for a given data-set.  This is intended to be run as part of a Seurat v2 workflow, but is written so that it can be used in a general context.  If supplying a seurat object, FindVariableGenes must have been run, otherwise a set of genes on which to perform the principal component analysis must be supplied.
#'
#' Arbitrary normalisation requires equal splits of data.  To check that 50/50 split does not unde-estimate the number of PCs, it is useful to perform a series of "titrations" of tha data, to check the number of PCs is not sensative to the sampling depth.  There is no relationship between the optimum number of PCs for un-normalised data and normalised data, so it is best to live with the limitations of normalising data (50/50 split, need to titrate) than to do find the optimum for a type of normalisation you will not use in practice.
#' 
#' This function was taken from the code of Matthew Young ( https://github.com/constantAmateur/MCVR/blob/master/code.R , commit c4caf93), with some modifications to:
#' i) accept a Seurat object as input, ii) multithread, iii) use custom PCA implementations, iv) fix irlba error.
#'
#' @param se A Seurat object
#' @param normalisation Normalization function to use, with a Seurat object as both input and output
#' @param pcamethod PCA function to use, with a Seurat object as both input and output 
#' @param trainFrac Fraction of data to be used for training.
#' @param p Fraction of total molecules sampled in experiment.
#' @param tFracs Titration fractions.  A vector indicating how to sub-sample the data to test if results are sensative to the depth of sampling.  If NULL, just run one version with all the data.
#' @param nSplits The number of random splits of the data to performe.
#' @param maxPCs Check up to this many PCs.
#' @param errorMetric Type of error to calculate.  Options are mse for mean squared erorr, or poisson, which calculates something proportional to the mean poisson log likelihood.
#' @param poissonMin If the PCs predict a number of counts below this value, use this instead to prevent infinite log-likelihood and penalise predicting zeros/negative values. 
#' @param confInt Used in quantifying the sampling error.  The function returns any PC that is within this confidence interval of the  minimum mean error.
#' @param verbose Logical; whether to print info on progress
#' @param BPPARAM Used to multithread; defaults to serial.
#' @param ... Extra parameters passed to PCA function.
#' @return A list.  Contains the average error across all cells for each titration and split of the data and various summaries thereof. 
molecularCrossValidation = function(se, normalisation, pcamethod=NULL, trainFrac=0.5, p=0.01, tFracs=c(1,0.8,0.5), nSplits=5,
                                    maxPCs=NULL,errorMetric=c('mse','poisson'),poissonMin=1e-6,confInt=0.95,verbose=FALSE,
                                    BPPARAM=SerialParam(progressbar=verbose), ...){
  errorMetric = match.arg(errorMetric)
  if(!is(se, "Seurat") || length(VariableFeatures(se))==0)
    stop("`se` should be a Seurat object with computed variable features.")
  if(is.null(maxPCs)) maxPCs <- min(100, ncol(Embeddings(se)))
  if(is.null(tFracs))
    tFracs=1
  titrate = length(tFracs)>1 || tFracs<1
  #Convert to number std. errors
  nSEs = qnorm(confInt/2 +0.5)
  #Convert to sparse Matrix of a useful type
  dat <- GetAssayData(se, slot="counts")
  dat = as(dat,'dgTMatrix')
  #Check that input data are integer counts
  if(any(dat@x%%1!=0))
    stop("Input matrix must contain only integers.")
  #If normalisation is on, have to do equal splits so we don't need to calculate complicated scale factor
  if(!identical(normalisation,identity)){
    #Arbitrary normalisation, have to have 50/50 split.
    if(trainFrac!=0.5){
      warning("Arbitrary normalisation requires equal splits of data. Setting trainFrac=0.5")
      trainFrac=0.5
    }
    if(!titrate){
      warning("When performing arbitrary normalisation it is useful to perform several titrations (by setting tFracs) to ensure results are not sensative to depth.  Consider re-running with tFracs set.")
    }
    if(errorMetric=='poisson')
      warning("Poisson error is not appropriate for non-count data.  Ensure your normalisation produces counts if you wish to proceed with this error profile.")
  }
  #Work out alpha
  alpha = (1-trainFrac)/trainFrac
  titrates = list()
  for(tFrac in tFracs){
    if(verbose && titrate)
      message(sprintf("Running with %d%% of data",tFrac*100))
    #Calculate correction factor alpha
    mse = matrix(NA,nrow=nSplits,ncol=maxPCs-1)
    colnames(mse) = seq(2,maxPCs)
    rownames(mse) = paste0("Split",seq(nSplits))
    for(i in seq(nSplits)){
      if(verbose) message(sprintf("Performing split %d of %d",i,nSplits)) 
      tst <- partitionData(dat,alpha,p,tFrac)
      #Normalise data
      train <- normalisation(CreateSeuratObject(tst[[2]]))
      tst <- GetAssayData(normalisation(CreateSeuratObject(tst[[1]])),slot="scale.data")
      # in case genes get kicked out by the normalization procedure (e.g. sctransform):
      varGenes <- intersect(VariableFeatures(se), row.names(tst))
      varGenes <- intersect(varGenes, row.names(GetAssayData(train, slot="scale.data")))
      VariableFeatures(train) <- varGenes
      tst <- tst[varGenes,]
      #Run PCA on training data.
      if(verbose) message("Running PCA and cross-validating")
      if(is.null(pcamethod)){
        pca = irlba::prcomp_irlba(GetAssayData(train, slot="scale.data")[varGenes,],center=FALSE,n=maxPCs)
        rot = t(pca$rotation)
        embed = pca$x
      }else{
        pca <- pcamethod(train, dims=maxPCs, weight.by.var=TRUE, ...)
        rot <- t(Loadings(pca))
        embed <- Embeddings(pca)
      }
      mse[i,] <- unlist(bplapply(seq(maxPCs,2),BPPARAM=BPPARAM,FUN=function(k){
        tmp = embed[,seq(k)] %*% rot[seq(k),]
        if(!all(dim(tmp)==dim(tst))) tmp <- t(tmp)
        #Mean squared error
        if(errorMetric=='mse'){
          return(mean((tmp*alpha-tst)**2))
        }else if(errorMetric=='poisson'){
          a = as.vector(tmp*alpha)
          b = as.vector(tst)
          #Fix those that are below minimum
          a[a<poissonMin]=poissonMin
          #Calculate poisson negative log likelihood
          return(-1*mean(dpois(b,a,log=TRUE)))
        }
      }))
    }
    titrates[[length(titrates)+1]] = mse
  }
  #Now work out which PC did the best
  #Work out mean and limits
  lower = do.call(cbind,lapply(titrates,apply,2,min))
  means = do.call(cbind,lapply(titrates,apply,2,mean))
  upper = do.call(cbind,lapply(titrates,apply,2,max))
  #Get the minima for each
  mins = seq(2,maxPCs)[apply(means,2,which.min)]
  #Get the standard error of each estimate
  ses = lapply(titrates,function(e) apply(e,2,sd)/sqrt(nSplits))
  ses = do.call(cbind,ses)
  #Find out the PCs that are within x standard errors of the minimum mean error
  seMins = t(t(means-nSEs*ses) <= apply(means,2,min))
  seMins = apply(seMins,2,which)
  #Ensure formating the same
  if(length(titrates)==1){
    seMins = list(as.numeric(rownames(seMins)))
  }else{
    seMins = lapply(seMins,function(e) as.numeric(names(e)))
  }
  if(length(titrates)==1)
    titrates=titrates[[1]]
  out = list(errors = titrates,
             means = means,
             lower = lower,
             upper = upper,
             minimas = mins,
             mins1sd = seMins,
             titrations = tFracs
             )
  return(out)
}

#' Plots summary of MCV run
#' 
#' Plots normalised loss curves for each titration of data as returned by molecularCrossValidation.
#' 
#' @param mcv Output of molecularCrossValidation.
#' @param cols Colours to apply to each titration line.
#' @return Nothing.  Produces a plot.
plotMCV = function(mcv,cols=seq(ncol(mcv$lower))){
  pcs = as.numeric(rownames(mcv$means))
  titrates = mcv$errors
  if(!is.list(titrates))
    titrates = list(titrates)
  #Work out normalisation factors
  normFacs = apply(mcv$lower,2,min)
  #Define plot area
  yRange = c(1,max(t(mcv$upper)/normFacs))
  layout(matrix(c(1,1,2),nrow=3))
  plot(0,0,type='n',xlab='Number of PCs',ylab='avg error / min avg error',ylim=yRange,xlim=range(pcs),frame.plot=FALSE,main='Error profiles')
  #Make lower,middle and upper plots
  for(i in seq_along(titrates)){
    lines(pcs,mcv$lower[,i]/normFacs[i],col=cols[i],lty=2)
    lines(pcs,mcv$means[,i]/normFacs[i],col=cols[i])
    lines(pcs,mcv$upper[,i]/normFacs[i],col=cols[i],lty=2)
  }
  #abline(v=mcv$minimas,col=seq(ncol(x)))
  legend(mean(mcv$minimas),yRange[2],legend=paste0(round(mcv$titrations*100),'% Data (Min = ',mcv$minimas,')'),col=cols,lty=1)
  tryCatch({
    plot(mcv$titrations*100,mcv$minimas,pch=19,xlab='Data percentage',ylab='Optimal #PCs',main='Convergance',ylim=range(pcs[unlist(mcv$mins1sd)]))
    lines(mcv$titrations*100,mcv$minimas)
    #Add error bars
    arrows(mcv$titrations*100, pcs[sapply(mcv$mins1sd,min)], mcv$titrations*100,pcs[sapply(mcv$mins1sd,max)], length=0.05, angle=90, code=3)
  }, error=function(e) NULL)
}
