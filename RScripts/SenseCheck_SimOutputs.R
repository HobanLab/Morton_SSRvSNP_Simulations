# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SENSE CHECK FSC SIMULATION OUTPUTS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script analyzes the simulation outputs from fastSimcoal2 
# to see if different expectations of simulation results are met.
# It first converts the Arlequin output files from fastSimcoal2 to genind objects that can be processed by adegenet

# Then, it processes these genind objects to see if the results of simulations meet certain expectations
# 1. Simulations with more populations should have higher numbers of total alleles
# 2. Generate averages of the allele frequency proportions across simulation replicates
# 3. Simulations with higher migration rates should have lower Fst values
# 4. Heterozygosity values for each simulation replicates (these are just compared to empirical values)
# 5. Analysis of the allele frequency spectra 

# These checks are made for various sets of simulations: MSAT (nInd=1200 and nInd=4800), and DNA 
# (two different mutation rates: 1e-7, 1e-8; two different population sizes: nInd=1200 and nInd=4800)
# This script will take a very long time to run (due to large DNA datasets); therefore,
# it's coded to run in the background (hence all of the print commands).

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

# %%% NIND = 1200 %%% ----
print("%%% ANALYZING N1200 SCENARIOS %%%")
# Source the genind objects from previously run simulations, using readGeninds functions
# Microsatellites
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N1200_marker/data.MSAT/"))
# DNA
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N1200_marker/data.DNA/"))

# ---- SENSE CHECK ----
# 1. MORE ALLELES IN SCENARIOS WITH MORE POPULATIONS ----
print("%%% TOTAL NUMBER OF ALLELES")
# MSAT ----
print("---MSAT---")
print("1 population (low and high migration)")
mean(sapply(MSAT_01pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(MSAT_01pop_migHigh.genList, function(x) ncol(x@tab)))
print("4 populations (low and high migration)")
mean(sapply(MSAT_04pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(MSAT_04pop_migHigh.genList, function(x) ncol(x@tab)))
print("16 populations (low and high migration)")
mean(sapply(MSAT_16pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(MSAT_16pop_migHigh.genList, function(x) ncol(x@tab)))

# DNA ----
print("---DNA---")
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
# Note: this section reports the allele frequency proportions BEFORE garden assignment
# (random assignment of samples to gardnes happens in subsetAndResample.R)
# MSAT ----
print("---MSAT---")
print("1 population (low and high migration)")
apply(sapply(MSAT_01pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(MSAT_01pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
print("4 populations (low and high migration)")
apply(sapply(MSAT_04pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(MSAT_04pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
print("16 populations (low and high migration)")
apply(sapply(MSAT_16pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(MSAT_16pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)

# DNA ----
print("---DNA---")
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
# MSAT ----
print("---MSAT---")
print("4 populations (low and high migration)")
sapply(MSAT_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(MSAT_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(MSAT_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(MSAT_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# DNA ----
print("---DNA---")
print("4 populations (low and high migration)")
sapply(DNA_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(DNA_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(DNA_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE)) 
sapply(DNA_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# 4. HETEROZYGOSITY ----
print("%%% HETEROZYGOSITY")
# MSAT ----
print("---MSAT---")
print("1 population (low and high migration)")
sapply(MSAT_01pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(MSAT_01pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
print("4 populations (low and high migration)")
sapply(MSAT_04pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(MSAT_04pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(MSAT_16pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(MSAT_16pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))

# DNA ----
print("---DNA---")
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
# QUESTION: when we calculate allele frequencies, do we divide by the number of individuals in the
# population? Or, in the entire species? Currently, doing the entire species...
# MSAT ----
# Specify the directory to save histograms to
histDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20230524/N1200/MSAT/"
# Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
makeAlleleFreqHist_genList(MSAT_01pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_01pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_04pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_04pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_16pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_16pop_migHigh.genList, outDir = histDir)

# DNA ----
# Specify the directory to save histograms to
histDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20230524/N1200/DNA/"
# Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
makeAlleleFreqHist_genList(DNA_01pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_01pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_04pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_04pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_16pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_16pop_migHigh.genList, outDir = histDir)

# %%% NIND = 4800 %%% ----
print("%%% ANALYZING N4800 SCENARIOS %%%")
# Source the genind objects from previously run simulations, using readGeninds functions
# Microsatellites
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N4800_marker/data.MSAT/"))
# DNA
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N4800_marker/data.DNA/"))

# ---- SENSE CHECK ----
# 1. MORE ALLELES IN SCENARIOS WITH MORE POPULATIONS ----
print("%%% TOTAL NUMBER OF ALLELES")
# MSAT ----
print("---MSAT---")
print("1 population (low and high migration)")
mean(sapply(MSAT_01pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(MSAT_01pop_migHigh.genList, function(x) ncol(x@tab)))
print("4 populations (low and high migration)")
mean(sapply(MSAT_04pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(MSAT_04pop_migHigh.genList, function(x) ncol(x@tab)))
print("16 populations (low and high migration)")
mean(sapply(MSAT_16pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(MSAT_16pop_migHigh.genList, function(x) ncol(x@tab)))

# DNA ----
print("---DNA---")
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
# Note: this section reports the allele frequency proportions BEFORE garden assignment
# (random assignment of samples to gardnes happens in subsetAndResample.R)
print("---MSAT---")
print("1 population (low and high migration)")
apply(sapply(MSAT_01pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(MSAT_01pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
print("4 populations (low and high migration)")
apply(sapply(MSAT_04pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(MSAT_04pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
print("16 populations (low and high migration)")
apply(sapply(MSAT_16pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(MSAT_16pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)

# DNA ----
print("---DNA---")
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
# MSAT ----
print("---MSAT---")
print("4 populations (low and high migration)")
sapply(MSAT_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(MSAT_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(MSAT_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(MSAT_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# DNA ----
print("---DNA---")
print("4 populations (low and high migration)")
sapply(DNA_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(DNA_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(DNA_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE)) 
sapply(DNA_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# 4. HETEROZYGOSITY ----
print("%%% HETEROZYGOSITY")
# MSAT ----
print("---MSAT---")
print("1 population (low and high migration)")
sapply(MSAT_01pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(MSAT_01pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
print("4 populations (low and high migration)")
sapply(MSAT_04pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(MSAT_04pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(MSAT_16pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(MSAT_16pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))

# DNA ----
print("---DNA---")
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
# QUESTION: when we calculate allele frequencies, do we divide by the number of individuals in the
# population? Or, in the entire species? Currently, doing the entire species...
# MSAT ----
# Specify the directory to save histograms to
histDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20230524/N4800/MSAT/"
# Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
makeAlleleFreqHist_genList(MSAT_01pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_01pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_04pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_04pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_16pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(MSAT_16pop_migHigh.genList, outDir = histDir)

# DNA ----
# Specify the directory to save histograms to
histDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20230524/N4800/DNA/"
# Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
makeAlleleFreqHist_genList(DNA_01pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_01pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_04pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_04pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_16pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_16pop_migHigh.genList, outDir = histDir)

# %%% DNA: LOWER MUTATION DATASETS %%% ----
print("%%% ANALYZING DNA LOW MUTATION (1e-8) SCENARIOS %%%")
# Run simulations
# N1200
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N1200_lowMut/data.DNA/"), prefix="DNA_N1200_lowMut")
# N4800
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N4800_lowMut/data.DNA/"), prefix="DNA_N4800_lowMut")

# ---- SENSE CHECK ----
# 1. MORE ALLELES IN SCENARIOS WITH MORE POPULATIONS ----
print("%%% TOTAL NUMBER OF ALLELES")
# N1200 ----
print("---N1200---")
print("1 population (low and high migration)")
mean(sapply(DNA_N1200_lowMut_01pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(DNA_N1200_lowMut_01pop_migHigh.genList, function(x) ncol(x@tab)))
print("4 populations (low and high migration)")
mean(sapply(DNA_N1200_lowMut_04pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(DNA_N1200_lowMut_04pop_migHigh.genList, function(x) ncol(x@tab)))
print("16 populations (low and high migration)")
mean(sapply(DNA_N1200_lowMut_16pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(DNA_N1200_lowMut_16pop_migHigh.genList, function(x) ncol(x@tab)))
# N4800 ----
print("---N4800---")
print("1 population (low and high migration)")
mean(sapply(DNA_N4800_lowMut_01pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(DNA_N4800_lowMut_01pop_migHigh.genList, function(x) ncol(x@tab)))
print("4 populations (low and high migration)")
mean(sapply(DNA_N4800_lowMut_04pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(DNA_N4800_lowMut_04pop_migHigh.genList, function(x) ncol(x@tab)))
print("16 populations (low and high migration)")
mean(sapply(DNA_N4800_lowMut_16pop_migLow.genList, function(x) ncol(x@tab)))
mean(sapply(DNA_N4800_lowMut_16pop_migHigh.genList, function(x) ncol(x@tab)))

# 2. AVERAGE NUMBER OF ALLELES IN EACH FREQUENCY CATEGORY ----
print("%%% AVERAGE PROPORTION OF ALLELES IN EACH FREQUENCY CATEGORY")
# Note: this section reports the allele frequency proportions BEFORE garden assignment
# (random assignment of samples to gardnes happens in subsetAndResample.R)
# N1200 ----
print("---N1200---")
print("1 population (low and high migration)")
apply(sapply(DNA_N1200_lowMut_01pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(DNA_N1200_lowMut_01pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
print("4 populations (low and high migration)")
apply(sapply(DNA_N1200_lowMut_04pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(DNA_N1200_lowMut_04pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
print("16 populations (low and high migration)")
apply(sapply(DNA_N1200_lowMut_16pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(DNA_N1200_lowMut_16pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
# N4800 ----
print("---N4800---")
print("1 population (low and high migration)")
apply(sapply(DNA_N4800_lowMut_01pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(DNA_N4800_lowMut_01pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
print("4 populations (low and high migration)")
apply(sapply(DNA_N4800_lowMut_04pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(DNA_N4800_lowMut_04pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)
print("16 populations (low and high migration)")
apply(sapply(DNA_N4800_lowMut_16pop_migLow.genList, getTotalAlleleFreqProportions), 1, mean)
apply(sapply(DNA_N4800_lowMut_16pop_migHigh.genList, getTotalAlleleFreqProportions), 1, mean)

# 3. HIGHER FST FOR SCENARIOS WITH LOWER MIGRATION RATES ----
print("%%% FST")
# N1200 ----
print("---N1200---")
print("4 populations (low and high migration)")
sapply(DNA_N1200_lowMut_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(DNA_N1200_lowMut_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(DNA_N1200_lowMut_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE)) 
sapply(DNA_N1200_lowMut_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
# N4800 ----
print("---N4800---")
print("4 populations (low and high migration)")
sapply(DNA_N4800_lowMut_04pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(DNA_N4800_lowMut_04pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(DNA_N4800_lowMut_16pop_migLow.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE)) 
sapply(DNA_N4800_lowMut_16pop_migHigh.genList, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# 4. HETEROZYGOSITY ----
print("%%% HETEROZYGOSITY")
# N1200 ----
print("---N1200---")
print("1 population (low and high migration)")
sapply(DNA_N1200_lowMut_01pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(DNA_N1200_lowMut_01pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
print("4 populations (low and high migration)")
sapply(DNA_N1200_lowMut_04pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(DNA_N1200_lowMut_04pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(DNA_N1200_lowMut_16pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(DNA_N1200_lowMut_16pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
# N4800 ----
print("---N4800---")
print("1 population (low and high migration)")
sapply(DNA_N4800_lowMut_01pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(DNA_N4800_lowMut_01pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
print("4 populations (low and high migration)")
sapply(DNA_N4800_lowMut_04pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(DNA_N4800_lowMut_04pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
print("16 populations (low and high migration)")
sapply(DNA_N4800_lowMut_16pop_migLow.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))
sapply(DNA_N4800_lowMut_16pop_migHigh.genList, function(x) mean(c(Hs(x)), na.rm=TRUE))

# 5. ALLELE FREQUENCY SPECTRA ----
print("%%% ALLELE FREQUENCIES")
# QUESTION: when we calculate allele frequencies, do we divide by the number of individuals in the
# population? Or, in the entire species? Currently, doing the entire species...
# N1200 ----
print("---N1200---")
# Specify the directory to save histograms to
histDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20230524/lowMut/N1200/"
# Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
makeAlleleFreqHist_genList(DNA_N1200_lowMut_01pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N1200_lowMut_01pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N1200_lowMut_04pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N1200_lowMut_04pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N1200_lowMut_16pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N1200_lowMut_16pop_migHigh.genList, outDir = histDir)
# N4800 ----
print("---N4800---")
# Specify the directory to save histograms to
histDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/SimulationSummaries_20230524/lowMut/N4800/"
# Generate histograms of each simulation replicate in each genind list, and save to a png file in the directory
makeAlleleFreqHist_genList(DNA_N4800_lowMut_01pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N4800_lowMut_01pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N4800_lowMut_04pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N4800_lowMut_04pop_migHigh.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N4800_lowMut_16pop_migLow.genList, outDir = histDir)
makeAlleleFreqHist_genList(DNA_N4800_lowMut_16pop_migHigh.genList, outDir = histDir)
