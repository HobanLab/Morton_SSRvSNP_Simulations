rm(list = ls())
library(strataG)
library(adegenet)
library(stringr)
library(parallel)
library(RColorBrewer)
library(scales)
sim.wd <- 'C:/Users/gsalas/Documents/resampling_CIs/Code/'
setwd(sim.wd)
source('Morton_SSRvSNP_Simulations/RScripts/functions_SSRvSNP_Sim.R')
# Parallelism: specify number of cores to use
num_cores <- detectCores() - 4
# Specify number of resampling replicates, to use for all scenarios below
num_reps <- 5
# Pick plot colors (for all plots!)
plotColors <- c('red','red4','darkorange3','coral','purple')
# Flags for processing different datasets (marker=msat, dna; nInd=1200, 4800; different n_to_drop flags; DNA low mutation)
MSAT_Flag <- TRUE
nInd_1200_Flag <- TRUE
if(nInd_1200_Flag==TRUE){
  print('%%% ANALYZING NIND=1200 DATASET %%%')
  # %%% Read in simulations and process results ----
  # Source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_MSAT(paste0(sim.wd,'Datasets/MSAT_N1200/'))
}

gm_resamp_array_function <- function(insert_genind, num_loci, num_reps){
  # Create an empty array named 'resamp_category5loc' to store results
  resamp_category <- array(dim = c(nrow(insert_genind@tab),4,num_reps))
  # loop 25 times for sets (5 loci in one set) of randomly selected loci
  for (i in 1:num_reps) {
    # Randomly sample an amount of loci from the genind object based on user input
    samp_loc <- sample(locNames(insert_genind), size = num_loci, replace = FALSE)
    # Subset the genind object to include only the columns corresponding to the sampled loci
    gm.Wild.genind <- insert_genind[, loc = samp_loc]
    # declare objects
    # access the genind matrix that shows the type of alleles and quantity present among wild individuals
    wildSamp <- gm.Wild.genind@tab
    # calculate the sum of each column in the matrix, ignoring NA values.
    # identify the indices where the sum is not equal to zero, this indicates columns with 
    # variation in allele counts. Subset the original matrix by selecting only the columns identified
    # in the previous step, removing columns with no variation in allele counts
    wildSamp <- wildSamp[, which(colSums(wildSamp, na.rm = TRUE) != 0)]
    # calculating a vector of allele frequencies (sum of allele counts divided by number of haplotypes, or individuals * 2)
    wildComplete <- colSums(wildSamp, na.rm = TRUE) / (nrow(wildSamp) * 2)
    # Subset 'wildComplete' to include only the alleles with non-zero frequency
    wildSubset <- wildComplete[wildComplete > 0]
    # initialize vectors to store results 
    total <- vector(length = nrow(wildSamp))
    common <- vector(length = nrow(wildSamp))
    lowfreq <- vector(length = nrow(wildSamp))
    rare <- vector(length = nrow(wildSamp))
    # Loop for each sample size, ranging from 1 tto the number of loci in the wild population
    for (j in 1:nrow(wildSamp)) {
      # randomly select a subset of rows from the matrix, without replacement
      samp <- sample(nrow(wildSamp), size = j, replace = FALSE)
      # subset the original matrix, to inlcude only the rows randomly selected
      samp <- wildSamp[samp,]
      # Measure the proportion of allelic representation in that sample
      # Check if its the first iteration of the inner loop
      if (j == 1) {
        # Identify the names of alleles in the sample that are also present in 'wildSubset'. 
        # Calculate the proportion of alleles in the sample that are also present in 'wildSubset'.
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        # Identify the names of the alleles in the sample that are also present in the 'wildSubset' for alleles with frequency >0.1
        # Calculate the proportion of common alleles in the sample compared to the 'wildSubset' for alleles with frequency >0.1
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        # Identify the names of alleles in the sample that are also present in the 'wildSubset' for alleles with frequency >0.01 and frequency <0.1
        # Calculate the proportion of low frequency alleles in the sample compared to the 'wildSubset' for allels with frequency >0.01 and frequency <0.1
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        # Identify the names of alleles in the sample that are also present in the 'wildSubset' for allleles with frequency <0.02.
        # Calcualte the proportion of rare alleles in the sample compared to the 'wildSubset' for alleles with frequency <0.01.
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        # Calculate the proportion of alleles in the sample that are also present in the 'wildSubset' for the case where j is not 1
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        # Calculate the proportion of alleles of common alleles in the sample compared toe the 'wildSubset' where j is not 1
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        # Calculate the proportion of alleles of low frequency alleles in the sample compared to the 'wildSubset' for the case where j is not 1
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        # Calculate the proportion of rare alleles in the sample compared ot the 'wildSubset' for the case where j is not 1
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      # Combine the results (proportions) for each sample size into a matrix named 'categorymat_10loc'.
      categorymat <- cbind(total, common, lowfreq, rare)
      # extract the column names (categories) from the matrix
      category <- colnames(categorymat)
      # set row and column names for 'resamp_category10loc'
      dimnames(resamp_category) <- list(paste0("sample ", 1:nrow(categorymat)), category, paste0("replicate ", 1:num_reps))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category[, , i] <- categorymat
    }
  }
  return(resamp_category)
}

list(MSAT_01pop_migHigh.genList[[1]]@tab)
