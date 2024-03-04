# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SENSE CHECK FSC SIMULATION OUTPUTS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script analyzes the simulation outputs from fastSimcoal2 
# to see if different expectations of simulation results are met.
# It first converts the Arlequin output files from fastSimcoal2 to genind objects that can be processed by adegenet

# Then, it processes these genind objects to examine simulation results in 5 steps:
# 1. Simulations with more populations should have higher numbers of total alleles
# 2. Generate averages of the allele frequency proportions across simulation replicates
# 3. Simulations with higher migration rates should have lower Fst values
# 4. Heterozygosity values for each simulation replicates (these are compared to empirical values)
# 5. Analysis of the allele frequency spectra 

# These checks are made for various sets of simulations: MSAT (nInd=1200 and nInd=4800), and DNA 
# (two different mutation rates: 1e-7, 1e-8; two different population sizes: nInd=1200 and nInd=4800)

# NOTE: DNA datasets take a significiant time to process. Thus, this script is structured to run in the background.

# NOTE: for calculation of allele frequencies (Step 5, histograms), frequency is for  THE ENTIRE species 
# (that is, number of allele counts in the total species divided by total number of individuals). 
# Allele frequencies are NOT calculated according to the frequency strictly WITHIN a population.

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs;
# this is a file located on the RAID1 drive, for space reasons, and linked in the home directory)
sim.wd <- "/home/akoontz/Shared/SSRvSNP_Sim/Code/"
setwd(sim.wd)
# Read in relevant functions
source("RScripts/functions_SSRvSNP_Sim.R")
# Flags for analyzing MSAT or SNP simulations
Flag_MSAT <- TRUE
Flag_DNA <- FALSE

# %%% N1200 %%% ----
print("%%% ANALYZING N1200 SCENARIOS %%%")
# %%% MSAT ----
if(Flag_MSAT == TRUE){
  print("%%% MSAT")
  # Source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N1200_marker/data.MSAT/"))
  
  # 1. MORE ALLELES IN SCENARIOS WITH MORE POPULATIONS ----
  print("%%% TOTAL NUMBER OF ALLELES")
  print("1 population (low and high migration)")
  mean(sapply(MSAT_01pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(MSAT_01pop_migHigh.genList, function(x) ncol(x@tab)))
  print("4 populations (low and high migration)")
  mean(sapply(MSAT_04pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(MSAT_04pop_migHigh.genList, function(x) ncol(x@tab)))
  print("16 populations (low and high migration)")
  mean(sapply(MSAT_16pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(MSAT_16pop_migHigh.genList, function(x) ncol(x@tab)))
  
  # 2. AVERAGE NUMBER OF ALLELES IN EACH FREQUENCY CATEGORY ----
  print("%%% AVERAGE PROPORTION OF ALLELES IN EACH FREQUENCY CATEGORY")
  # Note: this section reports the allele frequency proportions BEFORE garden assignment
  # (random assignment of samples to gardens happens in subsetAndResample.R)
  print("1 population (low and high migration)")
  apply(sapply(MSAT_01pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(MSAT_01pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  print("4 populations (low and high migration)")
  apply(sapply(MSAT_04pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(MSAT_04pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  print("16 populations (low and high migration)")
  apply(sapply(MSAT_16pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(MSAT_16pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  
  # 3. HIGHER FST FOR SCENARIOS WITH LOWER MIGRATION RATES ----
  print("%%% FST")
  print("4 populations (low and high migration)")
  sapply(MSAT_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  sapply(MSAT_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  print("16 populations (low and high migration)")
  sapply(MSAT_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  sapply(MSAT_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  
  # 4. HETEROZYGOSITY ----
  print("%%% HETEROZYGOSITY")
  print("1 population (low and high migration)")
  sapply(MSAT_01pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(MSAT_01pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  print("4 populations (low and high migration)")
  sapply(MSAT_04pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(MSAT_04pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  print("16 populations (low and high migration)")
  sapply(MSAT_16pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(MSAT_16pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  
  # 5. ALLELE FREQUENCY SPECTRA ----
  print("%%% ALLELE FREQUENCIES")
  # Specify the directory to save histograms to
  histDir <- "/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20240304/N1200/MSAT/"
  # Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
  makeAlleleFreqHist_genList(MSAT_01pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_01pop_migHigh.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_04pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_04pop_migHigh.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_16pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_16pop_migHigh.genList, outDir = histDir)
}

# %%% DNA ----
if(Flag_DNA == TRUE){
  print("%%% DNA")
  # Source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N1200_marker/data.DNA/"))
  
  # 1. MORE ALLELES IN SCENARIOS WITH MORE POPULATIONS ----
  print("%%% TOTAL NUMBER OF ALLELES")
  print("1 population (low and high migration)")
  mean(sapply(DNA_01pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(DNA_01pop_migHigh.genList, function(x) ncol(x@tab)))
  print("4 populations (low and high migration)")
  mean(sapply(DNA_04pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(DNA_04pop_migHigh.genList, function(x) ncol(x@tab)))
  print("16 populations (low and high migration)")
  mean(sapply(DNA_16pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(DNA_16pop_migHigh.genList, function(x) ncol(x@tab)))
  
  # 2. AVERAGE NUMBER OF ALLELES IN EACH FREQUENCY CATEGORY ----
  print("%%% AVERAGE PROPORTION OF ALLELES IN EACH FREQUENCY CATEGORY")
  print("1 population (low and high migration)")
  apply(sapply(DNA_01pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(DNA_01pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  print("4 populations (low and high migration)")
  apply(sapply(DNA_04pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(DNA_04pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  print("16 populations (low and high migration)")
  apply(sapply(DNA_16pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(DNA_16pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  
  # 3. HIGHER FST FOR SCENARIOS WITH LOWER MIGRATION RATES ----
  print("%%% FST")
  print("4 populations (low and high migration)")
  sapply(DNA_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  sapply(DNA_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  print("16 populations (low and high migration)")
  sapply(DNA_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE)) 
  sapply(DNA_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  
  # 4. HETEROZYGOSITY ----
  print("%%% HETEROZYGOSITY")
  print("1 population (low and high migration)")
  sapply(DNA_01pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(DNA_01pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  print("4 populations (low and high migration)")
  sapply(DNA_04pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(DNA_04pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  print("16 populations (low and high migration)")
  sapply(DNA_16pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(DNA_16pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  
  # 5. ALLELE FREQUENCY SPECTRA ----
  print("%%% ALLELE FREQUENCIES")
  # Specify the directory to save histograms to
  histDir <- "/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20240304/N1200/DNA/"
  # Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
  makeAlleleFreqHist_genList(DNA_01pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_01pop_migHigh.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_04pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_04pop_migHigh.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_16pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_16pop_migHigh.genList, outDir = histDir)
}

# %%% N4800 %%% ----
print("%%% ANALYZING N4800 SCENARIOS %%%")
# %%% MSAT ----
if(Flag_MSAT == TRUE){
  print("%%% MSAT")
  # Source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N4800_marker/data.MSAT/"))
  
  # 1. MORE ALLELES IN SCENARIOS WITH MORE POPULATIONS ----
  print("%%% TOTAL NUMBER OF ALLELES")
  print("1 population (low and high migration)")
  mean(sapply(MSAT_01pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(MSAT_01pop_migHigh.genList, function(x) ncol(x@tab)))
  print("4 populations (low and high migration)")
  mean(sapply(MSAT_04pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(MSAT_04pop_migHigh.genList, function(x) ncol(x@tab)))
  print("16 populations (low and high migration)")
  mean(sapply(MSAT_16pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(MSAT_16pop_migHigh.genList, function(x) ncol(x@tab)))
  
  # 2. AVERAGE NUMBER OF ALLELES IN EACH FREQUENCY CATEGORY ----
  print("%%% AVERAGE PROPORTION OF ALLELES IN EACH FREQUENCY CATEGORY")
  # Note: this section reports the allele frequency proportions BEFORE garden assignment
  # (random assignment of samples to gardens happens in subsetAndResample.R)
  print("1 population (low and high migration)")
  apply(sapply(MSAT_01pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(MSAT_01pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  print("4 populations (low and high migration)")
  apply(sapply(MSAT_04pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(MSAT_04pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  print("16 populations (low and high migration)")
  apply(sapply(MSAT_16pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(MSAT_16pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  
  # 3. HIGHER FST FOR SCENARIOS WITH LOWER MIGRATION RATES ----
  print("%%% FST")
  print("4 populations (low and high migration)")
  sapply(MSAT_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  sapply(MSAT_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  print("16 populations (low and high migration)")
  sapply(MSAT_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  sapply(MSAT_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  
  # 4. HETEROZYGOSITY ----
  print("%%% HETEROZYGOSITY")
  print("1 population (low and high migration)")
  sapply(MSAT_01pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(MSAT_01pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  print("4 populations (low and high migration)")
  sapply(MSAT_04pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(MSAT_04pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  print("16 populations (low and high migration)")
  sapply(MSAT_16pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(MSAT_16pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  
  # 5. ALLELE FREQUENCY SPECTRA ----
  print("%%% ALLELE FREQUENCIES")
  # Specify the directory to save histograms to
  histDir <- "/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20240304/N4800/MSAT/"
  # Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
  makeAlleleFreqHist_genList(MSAT_01pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_01pop_migHigh.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_04pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_04pop_migHigh.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_16pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(MSAT_16pop_migHigh.genList, outDir = histDir)
}

# %%% DNA ----
if(Flag_DNA == TRUE){
  print("%%% DNA")
  # Source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N4800_marker/data.DNA/"))
  
  # 1. MORE ALLELES IN SCENARIOS WITH MORE POPULATIONS ----
  print("%%% TOTAL NUMBER OF ALLELES")
  print("1 population (low and high migration)")
  mean(sapply(DNA_01pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(DNA_01pop_migHigh.genList, function(x) ncol(x@tab)))
  print("4 populations (low and high migration)")
  mean(sapply(DNA_04pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(DNA_04pop_migHigh.genList, function(x) ncol(x@tab)))
  print("16 populations (low and high migration)")
  mean(sapply(DNA_16pop_migLow.genList, function(x) ncol(x@tab)))
  mean(sapply(DNA_16pop_migHigh.genList, function(x) ncol(x@tab)))
  
  # 2. AVERAGE NUMBER OF ALLELES IN EACH FREQUENCY CATEGORY ----
  print("%%% AVERAGE PROPORTION OF ALLELES IN EACH FREQUENCY CATEGORY")
  print("1 population (low and high migration)")
  apply(sapply(DNA_01pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(DNA_01pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  print("4 populations (low and high migration)")
  apply(sapply(DNA_04pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(DNA_04pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  print("16 populations (low and high migration)")
  apply(sapply(DNA_16pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
  apply(sapply(DNA_16pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
  
  # 3. HIGHER FST FOR SCENARIOS WITH LOWER MIGRATION RATES ----
  print("%%% FST")
  print("4 populations (low and high migration)")
  sapply(DNA_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  sapply(DNA_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  print("16 populations (low and high migration)")
  sapply(DNA_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE)) 
  sapply(DNA_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
  
  # 4. HETEROZYGOSITY ----
  print("%%% HETEROZYGOSITY")
  print("1 population (low and high migration)")
  sapply(DNA_01pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(DNA_01pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  print("4 populations (low and high migration)")
  sapply(DNA_04pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(DNA_04pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  print("16 populations (low and high migration)")
  sapply(DNA_16pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  sapply(DNA_16pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
  
  # 5. ALLELE FREQUENCY SPECTRA ----
  print("%%% ALLELE FREQUENCIES")
  # Specify the directory to save histograms to
  histDir <- "/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20240304/N4800/DNA/"
  # Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
  makeAlleleFreqHist_genList(DNA_01pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_01pop_migHigh.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_04pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_04pop_migHigh.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_16pop_migLow.genList, outDir = histDir)
  makeAlleleFreqHist_genList(DNA_16pop_migHigh.genList, outDir = histDir)
}
