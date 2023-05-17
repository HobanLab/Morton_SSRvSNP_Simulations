# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GENERATE FSC PARAMETERS USING STRATAG %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script contains the parameter values used to run MSAT and DNA (SNP) fastSimcoal2 simulations
# using the R package strataG. It declares the simulation parameters, runs the simulations, and then
# converts the Arelquin outputs to genind files.

# This script contains the parameter values that were explored during the simulation paramter
# optimization phase of the project; the simulation outputs generated in this script are contained in 
# the SimulationOutputs/ArchivedScenarios/Archived Scenarios folder. 

# For the simulations used in the study, refer to the GenerateFSCparams.R script.

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

# ---- VARIABLES ----
num_reps <- 5
fscVersion <- "fsc2709"
# DEMES
# Specify number of total individuals, for all simulations
# Since there are 4 deme and 16 deme scenarios, this value must be divisible by 4 and 16
nInd <- 1200
# 1 Population
demeA <- fscDeme(deme.size = nInd, sample.size = nInd)
demes1 <- fscSettingsDemes(demeA)
# 4 Populations
demeB <- fscDeme(deme.size = nInd/4, sample.size = nInd/4)
demes4 <- fscSettingsDemes(demeB, demeB, demeB, demeB)
# 16 Populations
demeC <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demes16 <- fscSettingsDemes(demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,
                            demeC,demeC,demeC,demeC)
# MIGRATION
low_mig <- 0.001
high_mig <- 0.01
# 4 Populations
mig.mat.4.Low <- matrix(low_mig, nrow=4, ncol = 4); diag(mig.mat.4.Low) <- 0
mig.mat.4.High <- matrix(high_mig, nrow=4, ncol = 4); diag(mig.mat.4.High) <- 0
mig.mat.4.Final <- matrix(0, nrow=4, ncol = 4)
mig4Low <- fscSettingsMigration(mig.mat.4.Low, mig.mat.4.Final)
mig4High <- fscSettingsMigration(mig.mat.4.High, mig.mat.4.Final)
# 16 Populations
mig.mat.16.Low <- matrix(low_mig, nrow=16, ncol = 16); diag(mig.mat.16.Low) <- 0
mig.mat.16.High <- matrix(high_mig, nrow=16, ncol = 16); diag(mig.mat.16.High) <- 0
mig.mat.16.Final <- matrix(0, nrow=16, ncol = 16)
mig16Low <- fscSettingsMigration(mig.mat.16.Low, mig.mat.16.Final)
mig16High <- fscSettingsMigration(mig.mat.16.High, mig.mat.16.Final)

# HISTORICAL EVENTS
# 4 Populations
hist.event0 <- fscEvent(event.time = 50000, source = 0, sink = 0, prop.migrants = 0, migr.mat = 1)
hist.event1 <- fscEvent(event.time = 50000, source = 1, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event2 <- fscEvent(event.time = 50000, source = 2, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event3 <- fscEvent(event.time = 50000, source = 3, sink = 0, prop.migrants = 1, migr.mat = 1)
histEvent4 <- fscSettingsEvents(hist.event0, hist.event1, hist.event2, hist.event3)
# 16 Populations
hist.event4 <- fscEvent(event.time = 50000, source = 4, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event5 <- fscEvent(event.time = 50000, source = 5, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event6 <- fscEvent(event.time = 50000, source = 6, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event7 <- fscEvent(event.time = 50000, source = 7, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event8 <- fscEvent(event.time = 50000, source = 8, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event9 <- fscEvent(event.time = 50000, source = 9, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event10 <- fscEvent(event.time = 50000, source = 10, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event11 <- fscEvent(event.time = 50000, source = 11, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event12 <- fscEvent(event.time = 50000, source = 12, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event13 <- fscEvent(event.time = 50000, source = 13, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event14 <- fscEvent(event.time = 50000, source = 14, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event15 <- fscEvent(event.time = 50000, source = 15, sink = 0, prop.migrants = 1, migr.mat = 1)
histEvent16 <- fscSettingsEvents(hist.event0,hist.event1,hist.event2,hist.event3,hist.event4,hist.event5,hist.event6,
                                 hist.event7,hist.event8,hist.event9,hist.event10,hist.event11,hist.event12,
                                 hist.event13,hist.event14,hist.event15)
# GENETIC PARAMETERS
# MSAT
# Higher mutation rate (1e-3) and range constraint (15) to match Empirical distributions
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 1e-3, range.constraint = 15)
MSATgenetics <- fscSettingsGenetics(msats, num.chrom = 20)
# DNA
# Mutation rate set to match Empirical R80 distributions
# Large sequence length/number of blocks/number of chromosomes set to generate 1,000 SNP loci
dna_mutRate <- 1e-7
dna <- fscBlock_dna(sequence.length = 500, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, dna, dna, dna, dna, dna, dna, 
                                   num.chrom = 100)
# ---- MSAT ----
# %%% LOW RANGE (Range constraint = 10) ----
# Outputs are stored within a folder in the parent directory named "MSAT_range10_marker"
msat.wd <- paste0(sim.wd,"SimulationOutputs/ArchivedScenarios/MSAT_range10_marker/")
setwd(msat.wd)
# MSAT Genetic parameters
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 5e-4, range.constraint = 10)
MSATgenetics <- fscSettingsGenetics(msats, num.chrom = 20)

# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
MSAT_01pop_migLow.params <- fscWrite(demes = demes1, genetics = MSATgenetics,
                                     label = "MSAT_01pop_migLow", use.wd=TRUE)
MSAT_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = MSATgenetics, label = "MSAT_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 1 population, low migration")
MSAT_01pop_migLow.params <- fscRun(MSAT_01pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 1 population, high migration")
MSAT_01pop_migHigh.params <- fscRun(MSAT_01pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
MSAT_01pop_migLow.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_01pop_migLow/"),
                                          params = MSAT_01pop_migLow.params)
MSAT_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_01pop_migHigh/"),
                                           params = MSAT_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(MSAT_01pop_migLow.genind, file = paste0("data.MSAT/genind.MSAT_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_01pop_migLow.params, file = paste0("data.MSAT/params.MSAT_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_01pop_migHigh.genind, file = paste0("data.MSAT/genind.MSAT_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(MSAT_01pop_migHigh.params, file = paste0("data.MSAT/params.MSAT_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
MSAT_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                     genetics = MSATgenetics, label = "MSAT_04pop_migLow", use.wd=TRUE)
MSAT_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                      genetics = MSATgenetics, label = "MSAT_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 4 populations, low migration")
MSAT_04pop_migLow.params <- fscRun(MSAT_04pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 4 populations, high migration")
MSAT_04pop_migHigh.params <- fscRun(MSAT_04pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
MSAT_04pop_migLow.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_04pop_migLow/"),
                                          params = MSAT_04pop_migLow.params)
MSAT_04pop_migHigh.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_04pop_migHigh/"),
                                           params = MSAT_04pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(MSAT_04pop_migLow.genind, file = paste0("data.MSAT/genind.MSAT_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_04pop_migLow.params, file = paste0("data.MSAT/params.MSAT_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_04pop_migHigh.genind, file = paste0("data.MSAT/genind.MSAT_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(MSAT_04pop_migHigh.params, file = paste0("data.MSAT/params.MSAT_04pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
MSAT_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16,
                                     genetics = MSATgenetics, label = "MSAT_16pop_migLow", use.wd=TRUE)
MSAT_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16,
                                      genetics = MSATgenetics, label = "MSAT_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 16 populations, low migration")
MSAT_16pop_migLow.params <- fscRun(MSAT_16pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 16 populations, high migration")
MSAT_16pop_migHigh.params <- fscRun(MSAT_16pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
MSAT_16pop_migLow.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_16pop_migLow/"),
                                          params = MSAT_16pop_migLow.params)
MSAT_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_16pop_migHigh/"),
                                           params = MSAT_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(MSAT_16pop_migLow.genind, file = paste0("data.MSAT/genind.MSAT_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_16pop_migLow.params, file = paste0("data.MSAT/params.MSAT_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_16pop_migHigh.genind, file = paste0("data.MSAT/genind.MSAT_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(MSAT_16pop_migHigh.params, file = paste0("data.MSAT/params.MSAT_16pop_migHigh.",Sys.Date(),".Rdata"))

# %%% HIGH RANGE (Range constraint = 15) ----
msat.wd <- paste0(sim.wd,"SimulationOutputs/ArchivedScenarios/MSAT_range15_marker/")
setwd(msat.wd)
# MSAT Genetic parameters
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 5e-4, range.constraint = 15)
MSATgenetics <- fscSettingsGenetics(msats, num.chrom = 20)

# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
MSAT_01pop_migLow.params <- fscWrite(demes = demes1, genetics = MSATgenetics,
                                     label = "MSAT_01pop_migLow", use.wd=TRUE)
MSAT_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = MSATgenetics, label = "MSAT_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 1 population, low migration")
MSAT_01pop_migLow.params <- fscRun(MSAT_01pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 1 population, high migration")
MSAT_01pop_migHigh.params <- fscRun(MSAT_01pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
MSAT_01pop_migLow.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_01pop_migLow/"),
                                          params = MSAT_01pop_migLow.params)
MSAT_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_01pop_migHigh/"),
                                           params = MSAT_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(MSAT_01pop_migLow.genind, file = paste0("data.MSAT/genind.MSAT_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_01pop_migLow.params, file = paste0("data.MSAT/params.MSAT_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_01pop_migHigh.genind, file = paste0("data.MSAT/genind.MSAT_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(MSAT_01pop_migHigh.params, file = paste0("data.MSAT/params.MSAT_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
MSAT_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                     genetics = MSATgenetics, label = "MSAT_04pop_migLow", use.wd=TRUE)
MSAT_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                      genetics = MSATgenetics, label = "MSAT_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 4 populations, low migration")
MSAT_04pop_migLow.params <- fscRun(MSAT_04pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 4 populations, high migration")
MSAT_04pop_migHigh.params <- fscRun(MSAT_04pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
MSAT_04pop_migLow.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_04pop_migLow/"),
                                          params = MSAT_04pop_migLow.params)
MSAT_04pop_migHigh.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_04pop_migHigh/"),
                                           params = MSAT_04pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(MSAT_04pop_migLow.genind, file = paste0("data.MSAT/genind.MSAT_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_04pop_migLow.params, file = paste0("data.MSAT/params.MSAT_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_04pop_migHigh.genind, file = paste0("data.MSAT/genind.MSAT_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(MSAT_04pop_migHigh.params, file = paste0("data.MSAT/params.MSAT_04pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
MSAT_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16,
                                     genetics = MSATgenetics, label = "MSAT_16pop_migLow", use.wd=TRUE)
MSAT_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16,
                                      genetics = MSATgenetics, label = "MSAT_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 16 populations, low migration")
MSAT_16pop_migLow.params <- fscRun(MSAT_16pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 16 populations, high migration")
MSAT_16pop_migHigh.params <- fscRun(MSAT_16pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
MSAT_16pop_migLow.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_16pop_migLow/"),
                                          params = MSAT_16pop_migLow.params)
MSAT_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(msat.wd, "MSAT_16pop_migHigh/"),
                                           params = MSAT_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(MSAT_16pop_migLow.genind, file = paste0("data.MSAT/genind.MSAT_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_16pop_migLow.params, file = paste0("data.MSAT/params.MSAT_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_16pop_migHigh.genind, file = paste0("data.MSAT/genind.MSAT_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(MSAT_16pop_migHigh.params, file = paste0("data.MSAT/params.MSAT_16pop_migHigh.",Sys.Date(),".Rdata"))

# ---- DNA ----
# ---- 1,000 LOCI ----
# MEDIUM HIGH MUTATION ----
# Outputs are stored within a folder on the RAID1 directory, to allow for larger file sizes
dna_TEST.wd <- paste0(sim.wd,"SimulationOutputs/ArchivedScenarios/DNA_marker/")
setwd(dna_TEST.wd)
# GENETIC PARAMETERS
# Medium high mutation rate, increased sequence length/blocks/chromsomes
dna_mutRate <- 5e-6
dna <- fscBlock_dna(sequence.length = 500, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, dna, dna, dna, dna, dna, dna,
                                  num.chrom = 100)
# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migLow", use.wd=TRUE)
DNA_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 1 population, low migration")
DNA_01pop_migLow.params <- fscRun(DNA_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_01pop_migHigh.params <- fscRun(DNA_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_01pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_01pop_migLow/"),
                                        params = DNA_01pop_migLow.params)
DNA_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_01pop_migHigh/"),
                                         params = DNA_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_01pop_migLow.genind, file = paste0("data.DNA/genind.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
# 4 POPULATIONS
# Write parameter files
DNA_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                   genetics = DNAgenetics, label = "DNA_04pop_migLow", use.wd=TRUE)
DNA_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                    genetics = DNAgenetics, label = "DNA_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 4 populations, low migration")
DNA_04pop_migLow.params <- fscRun(DNA_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_04pop_migHigh.params <- fscRun(DNA_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_04pop_migLow.TEST.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_04pop_migLow/"),
                                        params = DNA_04pop_migLow.params)
DNA_04pop_migHigh.TEST.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_04pop_migHigh/"),
                                         params = DNA_04pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_04pop_migLow.TEST.genind, file = paste0("data.DNA/genind.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.TEST.genind, file = paste0("data.DNA/genind.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
# 16 POPULATIONS
# Write parameter files
DNA_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16,
                                   genetics = DNAgenetics, label = "DNA_16pop_migLow", use.wd=TRUE)
DNA_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16,
                                    genetics = DNAgenetics, label = "DNA_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 16 populations, low migration")
DNA_16pop_migLow.params <- fscRun(DNA_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_16pop_migHigh.params <- fscRun(DNA_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_16pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_16pop_migLow/"),
                                        params = DNA_16pop_migLow.params)
DNA_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_16pop_migHigh/"),
                                         params = DNA_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_16pop_migLow.genind, file = paste0("data.DNA/genind.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migLow.params, file = paste0("data.DNA/params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.params, file = paste0("data.DNA/params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))

# MEDIUM MUTATION ----
# Outputs are stored within a folder on the RAID1 directory, to allow for larger file sizes
dna_TEST.wd <- paste0(sim.wd,"SimulationOutputs/ArchivedScenarios/DNA_1000Loci_medMut_marker/")
setwd(dna_TEST.wd)
# GENETIC PARAMETERS
# Medium mutation rate, increased sequence length/blocks/chromosomes
dna_mutRate <- 5e-7
dna <- fscBlock_dna(sequence.length = 500, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, dna, dna, dna, dna, dna, dna,
                                   num.chrom = 100)
# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migLow", use.wd=TRUE)
DNA_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 1 population, low migration")
DNA_01pop_migLow.params <- fscRun(DNA_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_01pop_migHigh.params <- fscRun(DNA_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_01pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_01pop_migLow/"),
                                         params = DNA_01pop_migLow.params)
DNA_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_01pop_migHigh/"),
                                          params = DNA_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_01pop_migLow.genind, file = paste0("data.DNA/genind.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
DNA_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                    genetics = DNAgenetics, label = "DNA_04pop_migLow", use.wd=TRUE)
DNA_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                     genetics = DNAgenetics, label = "DNA_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 4 populations, low migration")
DNA_04pop_migLow.params <- fscRun(DNA_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_04pop_migHigh.params <- fscRun(DNA_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_04pop_migLow.TEST.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_04pop_migLow/"),
                                              params = DNA_04pop_migLow.params)
DNA_04pop_migHigh.TEST.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_04pop_migHigh/"),
                                               params = DNA_04pop_migHigh.params)

# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_04pop_migLow.TEST.genind, file = paste0("data.DNA/genind.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.TEST.genind, file = paste0("data.DNA/genind.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
DNA_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16,
                                    genetics = DNAgenetics, label = "DNA_16pop_migLow", use.wd=TRUE)
DNA_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16,
                                     genetics = DNAgenetics, label = "DNA_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 16 populations, low migration")
DNA_16pop_migLow.params <- fscRun(DNA_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_16pop_migHigh.params <- fscRun(DNA_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_16pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_16pop_migLow/"),
                                         params = DNA_16pop_migLow.params)
DNA_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_16pop_migHigh/"),
                                          params = DNA_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_16pop_migLow.genind, file = paste0("data.DNA/genind.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migLow.params, file = paste0("data.DNA/params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.params, file = paste0("data.DNA/params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))

# LOW MUTATION ----
# Outputs are stored within a folder on the RAID1 directory, to allow for larger file sizes
dna_TEST.wd <- paste0(sim.wd,"SimulationOutputs/ArchivedScenarios/DNA_1000Loci_lowMut_marker/")
setwd(dna_TEST.wd)
# GENETIC PARAMETERS
# Low mutation rate, increased sequence length/blocks/chromosomes
dna_mutRate <- 5e-8
dna <- fscBlock_dna(sequence.length = 500, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, dna, dna, dna, dna, dna, dna,
                                   num.chrom = 100)
# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migLow", use.wd=TRUE)
DNA_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 1 population, low migration")
DNA_01pop_migLow.params <- fscRun(DNA_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_01pop_migHigh.params <- fscRun(DNA_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_01pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_01pop_migLow/"),
                                         params = DNA_01pop_migLow.params)
DNA_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_01pop_migHigh/"),
                                          params = DNA_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_01pop_migLow.genind, file = paste0("data.DNA/genind.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
DNA_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                    genetics = DNAgenetics, label = "DNA_04pop_migLow", use.wd=TRUE)
DNA_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                     genetics = DNAgenetics, label = "DNA_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 4 populations, low migration")
DNA_04pop_migLow.params <- fscRun(DNA_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_04pop_migHigh.params <- fscRun(DNA_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_04pop_migLow.TEST.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_04pop_migLow/"),
                                              params = DNA_04pop_migLow.params)
DNA_04pop_migHigh.TEST.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_04pop_migHigh/"),
                                               params = DNA_04pop_migHigh.params)

# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_04pop_migLow.TEST.genind, file = paste0("data.DNA/genind.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.TEST.genind, file = paste0("data.DNA/genind.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
DNA_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16,
                                    genetics = DNAgenetics, label = "DNA_16pop_migLow", use.wd=TRUE)
DNA_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16,
                                     genetics = DNAgenetics, label = "DNA_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 16 populations, low migration")
DNA_16pop_migLow.params <- fscRun(DNA_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_16pop_migHigh.params <- fscRun(DNA_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_16pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_16pop_migLow/"),
                                         params = DNA_16pop_migLow.params)
DNA_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_TEST.wd, "DNA_16pop_migHigh/"),
                                          params = DNA_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_16pop_migLow.genind, file = paste0("data.DNA/genind.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migLow.params, file = paste0("data.DNA/params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.params, file = paste0("data.DNA/params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))

# ---- LOW NUMBERS OF LOCI (~20) ----
# %%% LOW MUTATION RATE (1E-8) ----
# Outputs are stored within a folder in the parent directory named "DNA_marker"
dna.wd <- paste0(sim.wd,"SimulationOutputs/ArchivedScenarios/DNA_marker/")
setwd(dna.wd)
# DNA Genetic parameters
# Variable for DNA mutation rate
dna_mutRate <- 1e-8
dna <- fscBlock_dna(sequence.length = 25, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, num.chrom = 5)

# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migLow", use.wd=TRUE)
DNA_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 1 population, low migration")
DNA_01pop_migLow.params <- fscRun(DNA_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_01pop_migHigh.params <- fscRun(DNA_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_01pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_01pop_migLow/"),
                                         params = DNA_01pop_migLow.params)
DNA_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_01pop_migHigh/"),
                                          params = DNA_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_01pop_migLow.genind, file = paste0("data.DNA/genind.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
DNA_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                    genetics = DNAgenetics, label = "DNA_04pop_migLow", use.wd=TRUE)
DNA_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                     genetics = DNAgenetics, label = "DNA_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 4 populations, low migration")
DNA_04pop_migLow.params <- fscRun(DNA_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_04pop_migHigh.params <- fscRun(DNA_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_04pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_04pop_migLow/"),
                                         params = DNA_04pop_migLow.params)
DNA_04pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_04pop_migHigh/"),
                                          params = DNA_04pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_04pop_migLow.genind, file = paste0("data.DNA/genind.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migLow.params, file = paste0("data.DNA/params.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.params, file = paste0("data.DNA/params.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
DNA_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16,
                                    genetics = DNAgenetics, label = "DNA_16pop_migLow", use.wd=TRUE)
DNA_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16,
                                     genetics = DNAgenetics, label = "DNA_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 16 populations, low migration")
DNA_16pop_migLow.params <- fscRun(DNA_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_16pop_migHigh.params <- fscRun(DNA_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_16pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_16pop_migLow/"),
                                         params = DNA_16pop_migLow.params)
DNA_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_16pop_migHigh/"),
                                          params = DNA_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_16pop_migLow.genind, file = paste0("data.DNA/genind.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migLow.params, file = paste0("data.DNA/params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.params, file = paste0("data.DNA/params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))

# %%% MEDIUM MUTATION RATE (1E-6) ----
# Outputs are stored within a folder in the parent directory named "DNA_highMut_marker"
dna_medMut.wd <- paste0(sim.wd,"SimulationOutputs/ArchivedScenarios/DNA_medMut_marker/")
setwd(dna_medMut.wd)
# DNA Genetic parameters
# Variable for DNA mutation rate
dna_mutRate <- 1e-6
dna <- fscBlock_dna(sequence.length = 25, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, num.chrom = 5)

# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_medMut_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_medMut_01pop_migLow", use.wd=TRUE)
DNA_medMut_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_medMut_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 1 population, low migration")
DNA_medMut_01pop_migLow.params <- fscRun(DNA_medMut_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_medMut_01pop_migHigh.params <- fscRun(DNA_medMut_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_medMut_01pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_medMut.wd, "DNA_medMut_01pop_migLow/"),
                                                params = DNA_medMut_01pop_migLow.params)
DNA_medMut_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_medMut.wd, "DNA_medMut_01pop_migHigh/"),
                                                 params = DNA_medMut_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_medMut_01pop_migLow.genind, file = paste0("data.DNA/genind.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_01pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_01pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_01pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
DNA_medMut_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                           genetics = DNAgenetics, label = "DNA_medMut_04pop_migLow", use.wd=TRUE)
DNA_medMut_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                            genetics = DNAgenetics, label = "DNA_medMut_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 4 populations, low migration")
DNA_medMut_04pop_migLow.params <- fscRun(DNA_medMut_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_medMut_04pop_migHigh.params <- fscRun(DNA_medMut_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_medMut_04pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_medMut.wd, "DNA_medMut_04pop_migLow/"),
                                                params = DNA_medMut_04pop_migLow.params)
DNA_medMut_04pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_medMut.wd, "DNA_medMut_04pop_migHigh/"),
                                                 params = DNA_medMut_04pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_medMut_04pop_migLow.genind, file = paste0("data.DNA/genind.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_04pop_migLow.params, file = paste0("data.DNA/params.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_04pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_04pop_migHigh.params, file = paste0("data.DNA/params.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
DNA_medMut_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16,
                                           genetics = DNAgenetics, label = "DNA_medMut_16pop_migLow", use.wd=TRUE)
DNA_medMut_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16,
                                            genetics = DNAgenetics, label = "DNA_medMut_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 16 populations, low migration")
DNA_medMut_16pop_migLow.params <- fscRun(DNA_medMut_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_medMut_16pop_migHigh.params <- fscRun(DNA_medMut_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_medMut_16pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_medMut.wd, "DNA_medMut_16pop_migLow/"),
                                                params = DNA_medMut_16pop_migLow.params)
DNA_medMut_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_medMut.wd, "DNA_medMut_16pop_migHigh/"),
                                                 params = DNA_medMut_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_medMut_16pop_migLow.genind, file = paste0("data.DNA/genind.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_16pop_migLow.params, file = paste0("data.DNA/params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_16pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_medMut_16pop_migHigh.params, file = paste0("data.DNA/params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))


# %%% HIGH MUTATION RATE (1E-5) ----
# Outputs are stored within a folder in the parent directory named "DNA_highMut_marker"
dna_highMut.wd <- paste0(sim.wd,"SimulationOutputs/ArchivedScenarios/DNA_highMut_marker/")
setwd(dna_highMut.wd)
# DNA Genetic parameters
# Variable for DNA mutation rate
dna_mutRate <- 1e-5
dna <- fscBlock_dna(sequence.length = 25, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, num.chrom = 5)

# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_highMut_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_highMut_01pop_migLow", use.wd=TRUE)
DNA_highMut_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_highMut_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 1 population, low migration")
DNA_highMut_01pop_migLow.params <- fscRun(DNA_highMut_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_highMut_01pop_migHigh.params <- fscRun(DNA_highMut_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_highMut_01pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_highMut.wd, "DNA_highMut_01pop_migLow/"),
                                                 params = DNA_highMut_01pop_migLow.params)
DNA_highMut_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_highMut.wd, "DNA_highMut_01pop_migHigh/"),
                                                  params = DNA_highMut_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_highMut_01pop_migLow.genind, file = paste0("data.DNA/genind.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_01pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_01pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_01pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
DNA_highMut_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                            genetics = DNAgenetics, label = "DNA_highMut_04pop_migLow", use.wd=TRUE)
DNA_highMut_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                             genetics = DNAgenetics, label = "DNA_highMut_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 4 populations, low migration")
DNA_highMut_04pop_migLow.params <- fscRun(DNA_highMut_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_highMut_04pop_migHigh.params <- fscRun(DNA_highMut_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_highMut_04pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_highMut.wd, "DNA_highMut_04pop_migLow/"),
                                                 params = DNA_highMut_04pop_migLow.params)
DNA_highMut_04pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_highMut.wd, "DNA_highMut_04pop_migHigh/"),
                                                  params = DNA_highMut_04pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_highMut_04pop_migLow.genind, file = paste0("data.DNA/genind.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_04pop_migLow.params, file = paste0("data.DNA/params.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_04pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_04pop_migHigh.params, file = paste0("data.DNA/params.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
DNA_highMut_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16,
                                            genetics = DNAgenetics, label = "DNA_highMut_16pop_migLow", use.wd=TRUE)
DNA_highMut_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16,
                                             genetics = DNAgenetics, label = "DNA_highMut_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 16 populations, low migration")
DNA_highMut_16pop_migLow.params <- fscRun(DNA_highMut_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_highMut_16pop_migHigh.params <- fscRun(DNA_highMut_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_highMut_16pop_migLow.genind <- convertAllArp(arp.path = paste0(dna_highMut.wd, "DNA_highMut_16pop_migLow/"),
                                                 params = DNA_highMut_16pop_migLow.params)
DNA_highMut_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna_highMut.wd, "DNA_highMut_16pop_migHigh/"),
                                                  params = DNA_highMut_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage
saveRDS(DNA_highMut_16pop_migLow.genind, file = paste0("data.DNA/genind.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_16pop_migLow.params, file = paste0("data.DNA/params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_16pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_highMut_16pop_migHigh.params, file = paste0("data.DNA/params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))

