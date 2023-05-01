# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR SSRvSNP SIMULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used in the Simulations component of the SSRvSNP comparison project.
# Functions are split into three sections:
#   1. Functions used for procesing Arlequin/strataG files. These include functions
#      for converting Arlequin outputs to genind, and for reading in existing strataG
#      params objects and genind objects (so simulations don't need to be run multiple times)
#   2. Commands for measuring ex situ representation in the simulated genind files. These are
#      analogous to the ex situ representation functions used in the Empirical analysis, and 
#      remove missing data from objects (even though simulations generate no missing data).
#   3. Commands for running the resampling analyses, in the simulated genind files. These are
#      analogous to the ex situ representation functions used in the Empirical analysis, and 
#      remove missing data from objects (even though simulations generate no missing data).

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)
library(parallel)

# ---- FUNCTIONS ----
# PROCESSING ARLEQUIN/STRATAG FILES ----
# Function converting Arlequin output to a single genind object (through gtypes format)
strataG_arp2gen <- function(params, repNumber){
  # Extract marker type from params argument
  marker <- params$settings$genetics$fsc.type
  # Read in the Arlequin file, convert it to a gtype object, then to a genind object
  arp <- fscReadArp(params, sim=c(1,repNumber), marker = marker)
  gtype <- df2gtypes(arp, ploidy = 2)
  genind <- gtypes2genind(gtype)
  # In the "other" slot of the genind object, pass the name of the simulation scenario, and return
  genind@other <- list(params$label)
  return(genind)
}

# Function for converting all of the Arlequin files in a directory to genind, generating a list of genind objects
convertAllArp <- function(arp.path, params){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing simulation outputs
  setwd(arp.path)
  # Create an empty list object to receive list of genind.
  # The length of this list is the number of replicates, which is specified as a numeric vector
  gen.List <- vector("list",length=length(dir()[str_detect(dir(), pattern = ".arp")]))
  fscReps <- seq(1, length(gen.List))
  # Move up one directory, in order for the fscReadArp command (within strataG_arp2gen) to work
  setwd("..")
  # Convert all Arlequin files to a list of genind objects
  for(i in 1:length(gen.List)){
    genind.obj <- strataG_arp2gen(params, rep=i)
    gen.List[[i]] <- genind.obj
  }
  # Reset to original working directory, and return a list of genind objects
  setwd(original.wd)
  return(gen.List)
}

# Function for reading in MSAT strataG params files, in specified directory. Prefix specifies how to name the params variables
readParams_MSAT <- function(params.wd, prefix="MSAT"){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing strataG params objects
  setwd(params.wd)
  # Use assign to assign the value (2nd argument) to the string (1st argument) specified. Prefix argument
  # allows variable names to be specifiable, but all variables will end with the same string
  # The pos=1 argument allows the variable to be passed to the global environment (rather than being kept locally)
  # The ^ character in dir allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple params objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.params"), readRDS(
    dir(pattern = "^params.MSAT_01pop_migLow")[length(dir(pattern = "^params.MSAT_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.params"), readRDS(
    dir(pattern = "^params.MSAT_01pop_migHigh")[length(dir(pattern = "^params.MSAT_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.params"), readRDS(
    dir(pattern = "^params.MSAT_04pop_migLow")[length(dir(pattern = "^params.MSAT_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.params"), readRDS(
    dir(pattern = "^params.MSAT_04pop_migHigh")[length(dir(pattern = "^params.MSAT_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.params"), readRDS(
    dir(pattern = "^params.MSAT_16pop_migLow")[length(dir(pattern = "^params.MSAT_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.params"), readRDS(
    dir(pattern = "^params.MSAT_16pop_migHigh")[length(dir(pattern = "^params.MSAT_16pop_migHigh"))]), pos = 1)
  # Reset to original working directory
  setwd(original.wd)
}

# Function for reading in MSAT genind files, in specified directory. Prefix specifies how to name the genind variables
readGeninds_MSAT <- function(geninds.wd, prefix="MSAT"){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing genind objects
  setwd(geninds.wd)
  # Use assign to assign the value (2nd argument) to the string (1st argument) specified. Prefix argument
  # allows variable names to be specifiable, but all variables will end with the same string
  # The pos=1 argument allows the variable to be passed to the global environment (rather than being kept locally)
  # The ^ character in dir allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple genind objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.genList"), readRDS(
    dir(pattern = "^genind.MSAT_01pop_migLow")[length(dir(pattern = "^genind.MSAT_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.genList"), readRDS(
    dir(pattern = "^genind.MSAT_01pop_migHigh")[length(dir(pattern = "^genind.MSAT_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.genList"), readRDS(
    dir(pattern = "^genind.MSAT_04pop_migLow")[length(dir(pattern = "^genind.MSAT_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.genList"), readRDS(
    dir(pattern = "^genind.MSAT_04pop_migHigh")[length(dir(pattern = "^genind.MSAT_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.genList"), readRDS(
    dir(pattern = "^genind.MSAT_16pop_migLow")[length(dir(pattern = "^genind.MSAT_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.genList"), readRDS(
    dir(pattern = "^genind.MSAT_16pop_migHigh")[length(dir(pattern = "^genind.MSAT_16pop_migHigh"))]), pos = 1)
  # Reset to original working directory
  setwd(original.wd)
}

# Function for reading in DNA strataG params files, in specified directory. Prefix specifies how to name the params variables
readParams_DNA <- function(params.wd, prefix="DNA"){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing strataG params objects
  setwd(params.wd)
  # Use assign to assign the value (2nd argument) to the string (1st argument) specified. Prefix argument
  # allows variable names to be specifiable, but all variables will end with the same string
  # The pos=1 argument allows the variable to be passed to the global environment (rather than being kept locally)
  # The ^ character in dir allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple params objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.params"), readRDS(
    dir(pattern = "^params.DNA_01pop_migLow")[length(dir(pattern = "^params.DNA_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.params"), readRDS(
    dir(pattern = "^params.DNA_01pop_migHigh")[length(dir(pattern = "^params.DNA_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.params"), readRDS(
    dir(pattern = "^params.DNA_04pop_migLow")[length(dir(pattern = "^params.DNA_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.params"), readRDS(
    dir(pattern = "^params.DNA_04pop_migHigh")[length(dir(pattern = "^params.DNA_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.params"), readRDS(
    dir(pattern = "^params.DNA_16pop_migLow")[length(dir(pattern = "^params.DNA_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.params"), readRDS(
    dir(pattern = "^params.DNA_16pop_migHigh")[length(dir(pattern = "^params.DNA_16pop_migHigh"))]), pos = 1)
  # Reset to original working directory
  setwd(original.wd)
}

# Function for reading in DNA genind files, in specified directory. Prefix specifies how to name the genind variables
readGeninds_DNA <- function(geninds.wd, prefix="DNA"){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing genind objects
  setwd(geninds.wd)
  # Use assign to assign the value (2nd argument) to the string (1st argument) specified. Prefix argument
  # allows variable names to be specifiable, but all variables will end with the same string
  # The pos=1 argument allows the variable to be passed to the global environment (rather than being kept locally)
  # The ^ character in dir allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple genind objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.genList"), readRDS(
    dir(pattern = "^genind.DNA_01pop_migLow")[length(dir(pattern = "^genind.DNA_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.genList"), readRDS(
    dir(pattern = "^genind.DNA_01pop_migHigh")[length(dir(pattern = "^genind.DNA_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.genList"), readRDS(
    dir(pattern = "^genind.DNA_04pop_migLow")[length(dir(pattern = "^genind.DNA_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.genList"), readRDS(
    dir(pattern = "^genind.DNA_04pop_migHigh")[length(dir(pattern = "^genind.DNA_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.genList"), readRDS(
    dir(pattern = "^genind.DNA_16pop_migLow")[length(dir(pattern = "^genind.DNA_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.genList"), readRDS(
    dir(pattern = "^genind.DNA_16pop_migHigh")[length(dir(pattern = "^genind.DNA_16pop_migHigh"))]), pos = 1)
  # Reset to original working directory
  setwd(original.wd)
}

# EX SITU REPRESENTATION ----
# Function for randomly assigning a proportion of a genind matrix to a "garden" population (the rest get "wild")
assignGardenSamples <- function(genind.obj, proportion=0.2){
  # Create a vector to store population names (start with all samples being "wild")
  popIDs <- rep("wild",length=(nInd(genind.obj)))
  # Get the names of randomly sampled rows of the genind matrix, based on the proportion argument
  gardenSamples <- rownames(genind.obj@tab[sample(nrow(genind.obj@tab), 
                                                  size=nInd(genind.obj)*proportion, replace = FALSE),])
  # Assign randomly sampled rows as "garden"
  popIDs[which(rownames(genind.obj@tab) %in% gardenSamples)] <- "garden"
  # Assign pop values and return genind object
  pop(genind.obj) <- popIDs
  return(genind.obj)
}

# Function for generating a vector of wild allele frequencies from a genind object
getWildFreqs <- function(gen.obj, wholeValues=TRUE){
  # Build a vector of rows corresponding to wild individuals (those that do not have a population of "garden")
  wildRows <- which(pop(gen.obj)!="garden")
  # Build the wild allele frequency vector: colSums of alleles (removing NAs), 
  # divided by number of haplotypes (Ne*2). wholeValues argument determines whether 
  # to return whole percentages or fractions
  if(wholeValues==TRUE){
    wildFreqs <- colSums(gen.obj@tab[wildRows,], na.rm = TRUE)/(length(wildRows)*2)*100
  } else{
    # (Fractions are more useful when trying to generate histograms)
    wildFreqs <- colSums(gen.obj@tab[wildRows,], na.rm = TRUE)/(length(wildRows)*2)
  }
  return(wildFreqs)
}

# Function for generating a vector of total allele frequencies from a genind object
getTotalFreqs <- function(gen.obj){
  # Build allele frequency vector: colSums of alleles (removing NAs), divided by number of haplotypes (Ne*2)
  totalFreqs <- colSums(gen.obj@tab, na.rm = TRUE)/(nInd(gen.obj)*2)*100
  return(totalFreqs)
}

# Exploratory function for reporting the proprtion of alleles of each category, from a (wild) frequency vector
getWildAlleleFreqProportions <- function(gen.obj){
  # Build the wild allele frequency vector, using the getWildFreqs function
  wildFreqs <- getWildFreqs(gen.obj)
  # Very common
  veryCommonAlleles <- wildFreqs[which(wildFreqs > 10)]
  veryCommon_prop <- (length(veryCommonAlleles)/length(wildFreqs))*100
  # Low frequency
  lowFrequencyAlleles <- wildFreqs[which(wildFreqs < 10 & wildFreqs > 1)]
  lowFrequency_prop <- (length(lowFrequencyAlleles)/length(wildFreqs))*100
  # Rare
  rareAlleles <- wildFreqs[which(wildFreqs < 1)]
  rare_prop <- (length(rareAlleles)/length(wildFreqs))*100
  # Build list of proportions, and return
  freqProportions <- c(veryCommon_prop, lowFrequency_prop, rare_prop)
  names(freqProportions) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  return(freqProportions)
}

# Exploratory function for reporting the proprtion of alleles of each category, 
# from a frequency vector (of ALL alleles--garden AND wild)
getTotalAlleleFreqProportions <- function(gen.obj){
  # Build the total allele frequency vector, using the getTotalFreqs function
  totalFreqs <- getTotalFreqs(gen.obj)
  # Very common
  veryCommonAlleles <- totalFreqs[which(totalFreqs > 10)]
  veryCommon_prop <- (length(veryCommonAlleles)/length(totalFreqs))*100
  # Low frequency
  lowFrequencyAlleles <- totalFreqs[which(totalFreqs < 10 & totalFreqs > 1)]
  lowFrequency_prop <- (length(lowFrequencyAlleles)/length(totalFreqs))*100
  # Rare
  rareAlleles <- totalFreqs[which(totalFreqs < 1)]
  rare_prop <- (length(rareAlleles)/length(totalFreqs))*100
  # Build list of proportions, and return
  freqProportions <- c(veryCommon_prop, lowFrequency_prop, rare_prop)
  names(freqProportions) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  return(freqProportions)
}

# Function for summarizing allele frequency proportions across replicates
summarize_alleleFreqProportions <- function(freqProportions){
  # Calculate the mean allele frequency proportions across replicates using apply
  # (Rows are allele frequency categories, columns are replicates. So margin value is 1)
  means <- apply(freqProportions, 1, mean, na.rm=TRUE)
  # Calculate standard deviations
  stdevs <- apply(freqProportions, 1, sd, na.rm=TRUE)
  # Combine statistics into a matrix, and return
  freqPropStats <- cbind(means, stdevs)
  return(freqPropStats)
}

# Function for reporting representation rates, using a vector of allele frequencies and a sample matrix.
# This function assumes that the freqVector represents the absolute allele frequencies
# for the population of interest (typically, the entire wild population). Allele names 
# between the frequency vector and the sample matrix must correspond in order for values to be comparable. 
getAlleleCategories <- function(freqVector, sampleMat){
  # Determine how many total alleles in the sample matrix are found in the frequency vector 
  garden.total_Alleles <- length(which(names(freqVector) %in% colnames(sampleMat)))
  wild.total_Alleles <- length(freqVector)
  total_Percentage <- (garden.total_Alleles/wild.total_Alleles)*100
  # Very common alleles (greater than 10%)
  garden.vCom_Alleles <- length(which(names(which(freqVector > 10)) %in% colnames(sampleMat)))
  wild.vCom_Alleles <- length(which(freqVector > 10))
  vCom_Percentage <- (garden.vCom_Alleles/wild.vCom_Alleles)*100
  # Common alleles (greater than 5%)
  garden.com_Alleles <- length(which(names(which(freqVector > 5)) %in% colnames(sampleMat)))
  wild.com_Alleles <- length(which(freqVector > 5))
  com_Percentage <- (garden.com_Alleles/wild.com_Alleles)*100
  # Low frequency alleles (between 1% and 10%)
  garden.lowFr_Alleles <- length(which(names(which(freqVector < 10 & freqVector > 1)) %in% colnames(sampleMat)))
  wild.lowFr_Alleles <- length(which(freqVector < 10 & freqVector > 1))
  lowFr_Percentage <- (garden.lowFr_Alleles/wild.lowFr_Alleles)*100
  # Rare alleles (less than 1%)
  garden.rare_Alleles <- length(which(names(which(freqVector < 1)) %in% colnames(sampleMat)))
  wild.rare_Alleles <- length(which(freqVector < 1))
  rare_Percentage <- (garden.rare_Alleles/wild.rare_Alleles)*100
  # Concatenate values to vectors
  gardenAlleles <- c(garden.total_Alleles, garden.vCom_Alleles, garden.com_Alleles, garden.lowFr_Alleles, garden.rare_Alleles)
  wildAlleles <- c(wild.total_Alleles, wild.vCom_Alleles, wild.com_Alleles, wild.lowFr_Alleles, wild.rare_Alleles)
  repRates <- c(total_Percentage,vCom_Percentage,com_Percentage,lowFr_Percentage,rare_Percentage) 
  # Bind vectors to a matrix, name dimensions, and return
  exSituValues <- cbind(gardenAlleles, wildAlleles, repRates)
  rownames(exSituValues) <- c("Total","Very common (>10%)","Common (>5%)",
                              "Low frequency (1% -- 10%)","Rare (<1%)")
  colnames(exSituValues) <- c("Garden", "Wild", "Rate (%)")
  return(exSituValues)
}

# This function is a wrapper of getAlleleCategories, and takes as an argument a single genind object
# (containing both garden and wild samples to analyze). It processes that genind object to extract the
# objects it needs to calculate ex situ representation: a vector of wild allele frequencies, and a 
# sample matrix of garden samples. Although simulations don't generate missing data, absent alleles
# (those with frequencies or colSums of 0) are removed from the the frequency vector and sample matrix.
# Depending on the returnMat argument, the function returns either a vector (of just proportions)
# or a matrix (of raw allele counts and proportions)
exSitu_Rep <- function(gen.obj, returnMat = 0){
  # Generate numerical vectors corresponding to garden and wild rows
  garden.Rows <- which(gen.obj@pop == "garden")
  wild.Rows <- which(gen.obj@pop == "wild")
  # Build the wild allele frequency vector: sum the allele counts, and divide by the number of wild samples
  # (which is equal to the number of wild rows in the genind matrix) times 2 (assuming diploid individuals)
  wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(length(wild.Rows)*2)*100
  # Remove any missing alleles (those with frequencies of 0) from the wild allele frequency vector
  wildFreqs <- wildFreqs[which(wildFreqs != 0)]
  # Generate the sample matrix to analyze (i.e. matrix of garden samples).
  gardenMat <- gen.obj@tab[garden.Rows,]
  # Remove any missing alleles (those with colSums of 0) from the matrix of garden samples
  gardenMat <- gardenMat[,which(colSums(gardenMat, na.rm = TRUE) != 0)]
  # Calculate how many alleles (of each category) the garden samples represent
  repValues <- getAlleleCategories(freqVector=wildFreqs, sampleMat = gardenMat)
  if(returnMat==1){
    # If returnMat==1, return a matrix of absolute values (garden alleles, wild alleles, and proportion)
    return(repValues)
  } else {
    # Otherwise, subset matrix returned by getAlleleCategories to just ex situ representation proportions
    repValues <- repValues[,3]
  }
  return(repValues)
}

# Function for summarizing ex situ representation rates across replicates
summarize_exSituRepresentation <- function(repRates){
  # Calculate the mean ex situ representation rate across replicates using apply
  # (Rows are rate categories, columns are replicates. So margin value is 1)
  means <- apply(repRates, 1, mean, na.rm=TRUE)
  # Calculate standard deviations
  stdevs <- apply(repRates, 1, sd, na.rm=TRUE)
  # Combine statistics into a matrix, and return
  repStats <- cbind(means, stdevs)
  return(repStats)
}

# Wrapper function, which generates both allele frequency proportions and ex situ representation rates,
# for a list of genind objects
summarize_simulations <- function(gen.List){
  # Build array to capture allele frequency proportions
  alleleFreqSummaries <- array(dim = c(3, 2, length(gen.List)))
  rownames(alleleFreqSummaries) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Build array to capture ex situ representation rates
  repRateSummaries <- array(dim = c(5, 2, length(gen.List)))
  rownames(repRateSummaries) <- 
    c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  colnames(alleleFreqSummaries) <- colnames(repRateSummaries) <-c("mean", "sd")
  
  # Loop through list of genind objects, calculating metrics for each item
  for (i in 1:length(gen.List)){
    # Calculate and summarize allele frequency scenarios. Each array slot is a different scenario
    alleleFrequencies <- sapply(gen.List[[i]], getWildAlleleFreqProportions)
    alleleFreqSummaries[,,i] <- summarize_alleleFreqProportions(alleleFrequencies)
    # Calculate and summarize ex situ representation rates. Each array slot is a different scenario
    representationRates <- sapply(gen.List[[i]], exSitu_Rep)
    repRateSummaries[,,i] <- summarize_exSituRepresentation(representationRates)
  }
  # Round results to 2 digits
  alleleFreqSummaries <- round(alleleFreqSummaries, 2)  
  repRateSummaries <- round(repRateSummaries, 2)
  # Generate a list of the two arrays, and return
  return(list("alleleFrequencyProportions"=alleleFreqSummaries, "representationRates"=repRateSummaries))
}

# Exploratory function, which creates a histogram of allele frequencies from a genind object
makeAlleleFreqHist <- function(gen.obj){
  # Make a vector of allele frequency values, from the genind object
  wildAlleleFreqs <- getWildFreqs(gen.obj, wholeValues = FALSE)
  # Generate histogram
  hist(wildAlleleFreqs, main=gen.obj@other, freq=FALSE, breaks=seq(0, 1.0, 0.01))
}

# Exploratory function, which creates a histogram of allele frequencies from a list of genind objects
makeAlleleFreqHist_genList <- function(gen.List, outDir="~/Shared/SSRvSNP_Sim/Code/"){
  # Use a loop to process each genind item in the list
  for(i in 1:length(gen.List)){
    # Call the png command, to save histogram outputs. Save to specified directory
    png(filename = paste0(outDir, gen.List[[i]]@other, "_", i, ".png", width = 1262, height = 734))
    # Call makeAlleleFreqHist function, nested within invisible to avoid code printing to standard output
    invisible(makeAlleleFreqHist(gen.List[[i]]))
    dev.off()
  }
}

# RESAMPLING FUNCTIONS ----
# Ex situ sample function, which finds the level of ex situ representation of a sample of individuals.
# It takes a genind object (containing garden and wild samples) and an integer specifying the number of 
# rows (individuals) to draw from the matrix of wild samples. Using the sample and getAlleleCategories 
# functions, it calculates the allelic representation of randomly sampled rows from the  specified matrix. 
# Prior to passing the required frequency vector and sample matrix to getAlleleCategories, 
# missing loci are removed, even though simulations by default do not generate any missing data.
# Just the 3rd column of the matrix generated by getAlleleCategories (i.e. the representation rates) 
# is retained for resampling, since we don't use the absolute allele count values.
exSitu_Sample <- function(gen.obj, numSamples){
  # Build the wild allele frequency vector, using the getWildFreqs function
  freqVector <- getWildFreqs(gen.obj)
  # Remove any missing alleles (those with frequencies of 0) from the frequency vector
  freqVector <- freqVector[which(freqVector != 0)]
  # Create a matrix of wild individuals (those with population "wild") from genind object
  wildMat <- gen.obj@tab[which(pop(gen.obj) == "wild"),]
  # From a matrix of individuals, select a set of random individuals (rows)
  samp <- wildMat[sample(nrow(wildMat), size=numSamples, replace = FALSE),]
  # Remove any missing alleles (those with colSums of 0) from the sample matrix
  samp <- samp[,which(colSums(samp, na.rm = TRUE) != 0)]
  # Calculate how many alleles (of each category) that sample captures, and return
  repRates <- getAlleleCategories(freqVector, samp)
  # Subset matrix returned by getAlleleCategories to just 3rd column (representation rates), and return
  repRates <- repRates[,3]
  return(repRates)
}

# Wrapper for the exSitu_Sample function, iterating that function over all wild samples in a genind object
exSitu_Resample <- function(gen.obj){
  # Apply the exSituSample function number of times equal to number of wild samples,
  # excluding 1 (because we need at least 2 individuals to sample)
  # Resulting matrix needs to be transposed in order to keep columns as different allele categories
  representationMatrix <- t(sapply(2:length(which(pop(gen.obj)=="wild")), 
                                   function(x) exSitu_Sample(gen.obj, x)))
  # Name columns according to categories of allelic representation, and return matrix
  colnames(representationMatrix) <- c("Total","Very common","Common","Low frequency","Rare")
  return(representationMatrix)
}

# Wrapper for exSitu_Resample, which will generate an array of values from a single genind object
Resample_genind <- function(gen.obj, reps=5){
  # Run resampling for all replicates, using sapply and lambda function
  resamplingArray <- sapply(1:reps, function(x) exSitu_Resample(gen.obj = gen.obj), simplify = "array")
  # Rename third array dimension to describe simulation scenario (captured in the genind object), and return
  dimnames(resamplingArray)[[3]] <- rep(unlist(gen.obj@other), dim(resamplingArray)[[3]])
  return(resamplingArray)
}

# Wrapper for Resample_genind, which will generate an array of values from a list of genind objects
Resample_genList <- function(gen.List, reps=5){
  # Run resampling for all replicates, using lapply and lambda function, and return array
  resamplingArray <- lapply(gen.List, function(x) Resample_genind(gen.obj = x, reps=reps))
  return(resamplingArray)
}

# Parallel wrapper for exSitu_Resample, which will generate an array of values from a single genind object
# This function has been made obsolete by Resample_genlist (an alternative), but is still declared here,
# because this parallelized resampling approach can still work for N1200 datasets.
parResample_genind <- function(gen.obj, reps=5, cluster){
  # Run resampling for all replicates, using parSapply and lambda function, and return array
  resamplingArray <- parSapply(cl=cluster, 1:reps, 
                               function(x) exSitu_Resample(gen.obj = gen.obj), simplify = "array")
  # Rename third array dimension to describe simulation scenario (captured in the genind object), and return
  dimnames(resamplingArray)[[3]] <- rep(unlist(gen.obj@other), dim(resamplingArray)[[3]])
  return(resamplingArray)
}

# From resampling arrayy, calculate the mean minimum sample size to represent 95% of the total wild diversity
resample_min95_mean <- function(resamplingArray){
  # resampling array[,1,]: returns the Total column values for each replicate (3rd array dimension)
  # apply(resamplingArray[,1,],1,mean): calculates the average across replicates for each row
  # which(apply(resamplingArray[,1,],1,mean) > 95): returns the rows with averages greater than 95
  # min(which(apply(resamplingArray[,1,],1,mean) > 95)): the lowest row with an average greater than 95
  meanValue <- min(which(apply(resamplingArray[,1,],1, mean, na.rm=TRUE) > 95))
  # Get simulation scenario name from resamplingArray. This is to keep track of values later on
  names(meanValue) <- unique(dimnames(resamplingArray)[[3]])
  return(meanValue)
}

# From resampling array, calculate the standard deviation, at the mean 95% value
resample_min95_sd <- function(resamplingArray){
  # Determine the mean value for representing 95% of allelic diversity
  meanValue <- resample_min95_mean(resamplingArray)
  # Calculate the standard deviation, at that mean value
  sdValue <- apply(resamplingArray[,1,],1,sd)[meanValue]
  # Get simulation scenario name from resamplingArray, and return
  names(sdValue) <- unique(dimnames(resamplingArray)[[3]])
  return(sdValue)
}

# From resampling array, calculate the mean values (across resampling replicates) for each allele frequency category
resample_meanValues <- function(resamplingArray){
  # Declare a matrix to receive average values
  meanValue_mat <- matrix(nrow=nrow(resamplingArray), ncol=ncol(resamplingArray))
  # For each column in the array, average results across resampling replicates (3rd array dimension)
  for(i in 1:ncol(resamplingArray)){
    meanValue_mat[,i] <- apply(resamplingArray[,i,], 1, mean, na.rm=TRUE)
  }
  # Give names to meanValue_mat columns, and return
  colnames(meanValue_mat) <- c("Total","Very common","Common","Low frequency","Rare")
  return(meanValue_mat)
}

# # From resampling array, calculate the mean values (across resampling replicates) for each allele frequency category
# resample_meanValues_OLD <- function(resamplingArray){
#   # Declare a matrix to receive average values
#   meanValue_mat <- matrix(nrow=nrow(resamplingArray), ncol=ncol(resamplingArray))
#   # For each column in the array, average results across resampling replicates (3rd array dimension)
#   meanValue_mat[,1] <- apply(resamplingArray[,1,], 1, mean, na.rm=TRUE)
#   meanValue_mat[,2] <- apply(resamplingArray[,2,], 1, mean, na.rm=TRUE)
#   meanValue_mat[,3] <- apply(resamplingArray[,3,], 1, mean, na.rm=TRUE)
#   meanValue_mat[,4] <- apply(resamplingArray[,4,], 1, mean, na.rm=TRUE)
#   meanValue_mat[,5] <- apply(resamplingArray[,5,], 1, mean, na.rm=TRUE)
#   # Give names to meanValue_mat columns, and return
#   colnames(meanValue_mat) <- c("Total","Very common","Common","Low frequency","Rare")
#   return(meanValue_mat)
# }

# Summary plotting function, from array
resample_Plot <- function(resamplingArray, colors, largePopFlag=FALSE){
  # Create two vectors for colors. This is to show points on the graph and in the legend clearly
  fullColors <- colors
  fadedColors <- c(colors[1], alpha(colors[2:5], 0.2))
  # Generate the average values (across replicates) for each allele frequency category 
  averageValueMat <- resample_meanValues(resamplingArray)
  # Generate the minimum sample size to represent 95% of allelic diversity (across replicates)
  min95_Value <- resample_min95_mean(resamplingArray)
  # Use the matplot function to plot the matrix of average values, with specified settings
  matplot(averageValueMat, ylim=c(0,110), col=fadedColors, pch=16,
          xlab="Number of Individuals", ylab="Percent Diversity Representation",
          main=unique(dimnames(resamplingArray)[[3]]))
  # Mark the 95% threshold line, as well as the 95% minimum sampling size
  abline(h=95, col="black", lty=3); abline(v=min95_Value, col="black")
  # Depending on total population size, plot text and legends in different areas of the graph
  if(largePopFlag==FALSE){
    # nInd = 1200 graphs
    min95_Line <- min95_Value+200
    xLeg <- 950
    yLeg <- 87
  } else{
    # nInd = 4800 graphs
    min95_Line <- min95_Value+700
    xLeg <- 3500
    yLeg <- 60
  }
  # Add text for the minimum sampling size line
  mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value),
        side=1, line=-1.5, at=min95_Line)
  # Add legend
  legend(x=xLeg, y=yLeg, inset = 0.05,
         legend = c("Total","Very common","Common","Low frequency", "Rare"),
         col=fullColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
}

# Summary plotting function, from array. This function saves the plots generate to the disk (in .png format)
resample_Plot_PNG <- function(resamplingArray, colors, largePopFlag=FALSE, data.dir){
  # Extract the scenario name from the resampling array
  scenName <- unique(dimnames(resamplingArray)[[3]])
  # Call png command, to save resampling plot to disk. To determine a unique file name to save plot, 
  # find out how many PNG files already exist in specified folder.
  png(filename = paste0(data.dir, scenName, "_", length(list.files(path=data.dir, pattern = ".png"))+1, ".png"), 
      width = 1262, height = 734)
  # Create two vectors for colors. This is to show points on the graph and in the legend clearly
  fullColors <- colors
  fadedColors <- c(colors[1], alpha(colors[2:5], 0.2))
  # Generate the average values (across replicates) for each allele frequency category 
  averageValueMat <- resample_meanValues(resamplingArray)
  # Generate the minimum sample size to represent 95% of allelic diversity (across replicates)
  min95_Value <- resample_min95_mean(resamplingArray)
  # Use the matplot function to plot the matrix of average values, with specified settings
  matplot(averageValueMat, ylim=c(0,110), col=fadedColors, pch=16,
          xlab="Number of Individuals", ylab="Percent Diversity Representation", main=scenName)
  # Mark the 95% threshold line, as well as the 95% minimum sampling size
  abline(h=95, col="black", lty=3); abline(v=min95_Value, col="black")
  # Depending on total population size, plot text and legends in different areas of the graph
  if(largePopFlag==FALSE){
    # nInd = 1200 graphs
    min95_Line <- min95_Value+200
    xLeg <- 950
    yLeg <- 87
  } else{
    # nInd = 4800 graphs
    min95_Line <- min95_Value+700
    xLeg <- 3500
    yLeg <- 60
  }
  # Add text for the minimum sampling size line
  mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value),
        side=1, line=-1.5, at=min95_Line)
  # Add legend
  legend(x=xLeg, y=yLeg, inset = 0.05,
         legend = c("Total","Very common","Common","Low frequency", "Rare"),
         col=fullColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
  # Turn off plotting device
  dev.off()
}
