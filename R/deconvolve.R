##' MCMC sampling for deconvolution
##'
##' This generates a stan sample from the model
##' 
##' @title deconvole
##' @param X the signature: a genes x cells matrix
##' @param Y the signal: a genes x replicates matrix
##' @param tol the sd in the likelihood
##' @param iter number of iterations in the chains
##' @param cores number of cores used (e.g. = no. cores)
##' @param chains number of MCMC chains
##' @param ... 
##' @return a stan sample object
##' @author Pete Dodd
##' @import rstan
##' @export
deconvolve <- function(X,Y,tol=0.05,iter=1e3,cores=4,chains=4,...){
  ## ouput checks
  cat("No. outputs/genes = ", nrow(Y),"\n")
  cat("No. sources/cells = ", ncol(X), "\n")
  cat("No. replicates = ", ncol(Y), "\n")
  ## make data
  sdata <- list(
    NG = nrow(Y), # number of genes
    NC = ncol(X), # number of cells
    NR = ncol(Y), # number of reps
    Y = Y, # responses
    X = X, # signatures
    tol = tol # tolerance for data likelihood
  )
  ## sampling
  rstan::sampling(stanmodels$bdcF, data = sdata, chains = chains, cores = cores, iter = iter, ...)
}

##' Extract proportions from samples
##'
##' This generates summary statistics from MCMC samples.
##' 
##' @title Post-process
##' @param X stan sample object
##' @return a summary of samples
##' @author Pete Dodd
##' @import rstan
##' @import data.table
##' @export
postprocess.baydec <- function(X) {
  PZ <- rstan::extract(X, "P")
  PS <- data.table::as.data.table(PZ$P)
  names(PS)[names(PS)=="V2"] <- "cells"
  names(PS)[names(PS) == "V1"] <- "expts"
  PSE <- PS[, list(value = mean(value),
                value.lo = quantile(value,0.025),
                value.hi = quantile(value,1-0.025),
                value.sd = sd(value)),
            by = list(cells, expts)]
  PSE
}


