# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GENERATE FSC PARAMETERS USING STRATAG (DNA: LOW MUTATION) %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses strataG to create the fastsimcoal2 (fsc) parameter files used for the simulation component
# of the SSR v. SNP comparison project. Then, those parameter files are used to run fsc simulations. 
# Arlequin outputs are converted to genind, and both the strataG params objects and the 
# genind objects are saved to .Rdata files for long-term storage.

# This script is similar to the N1200/N4800 versions, but instead only simulates DNA markers, and uses
# a lower mutation rate (1e-8). This is to document the impact that mutation rate has on our findings.
# The code is divided into 2 sections: one for N1200 individuals, and one for N4800 individuals.

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
# This section specifies fsc simulation parameters as variables, for running simulations using strataG (below)

# Specify the number of simulation replicates, per fsc simulation scenario
# For each simulation replicate, a corresponding Arlequin (and, ultimately, genind) object will be created
num_reps <- 5
# Generate a string, indicating the version of fsc being utilized
fscVersion <- "fsc2709"

# DEMES
# Demes are specified for each section below (N1200 or N4800)

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
# DNA
# Lower mutation rate (compared to "original" simulations)
dna_mutRate <- 1e-8
dna <- fscBlock_dna(sequence.length = 500, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, dna, dna, dna, dna, dna, dna, 
                                   num.chrom = 100)
# ---- N1200 ----
# SPECIFIC VARIABLES
# Specify number of total individuals, for all simulations
# Since there are 4 deme and 16 deme scenarios, this value must be divisible by 4 and 16
nInd <- 1200
# DEMES
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
# Specify working directory
dna.wd <- paste0(sim.wd,"SimulationOutputs/DNA_N1200_lowMut/")
setwd(dna.wd)

# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_N1200_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, 
                                    label = "DNA_N1200_01pop_migLow", use.wd=TRUE)
DNA_N1200_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, 
                                     label = "DNA_N1200_01pop_migHigh", use.wd=TRUE)
# Run parameter files. 
# For DNA simulations, set all.sites=TRUE (in order to print monomorphic as well as polymorphic sequence sites)
print("DNA: 1 population, low migration")
DNA_N1200_01pop_migLow.params <- fscRun(DNA_N1200_01pop_migLow.params, num.sims = num_reps, 
                                  all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_N1200_01pop_migHigh.params <- fscRun(DNA_N1200_01pop_migHigh.params, num.sims = num_reps, 
                                   all.sites = TRUE, exec = fscVersion)
# Move .log and .par files into respective simulation folders
file.copy(from=paste0(dna.wd,"DNA_N1200_01pop_migLow.log"), to=paste0(dna.wd,"DNA_N1200_01pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_01pop_migLow.par"), to=paste0(dna.wd,"DNA_N1200_01pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_01pop_migHigh.log"), to=paste0(dna.wd,"DNA_N1200_01pop_migHigh/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_01pop_migHigh.par"), to=paste0(dna.wd,"DNA_N1200_01pop_migHigh/"))
file.remove(list.files(pattern = ".log"))
file.remove(list.files(pattern = ".par"))
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_N1200_01pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N1200_01pop_migLow/"), 
                                         params = DNA_N1200_01pop_migLow.params)
DNA_N1200_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N1200_01pop_migHigh/"), 
                                          params = DNA_N1200_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage. Note that genind objects 
# need to have a file name of genind.scenarioName in order to be processed by the readGendinds_* functions
saveRDS(DNA_N1200_01pop_migLow.genind, file = paste0("data.DNA/genind.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_01pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_01pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_01pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
DNA_N1200_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                    genetics = DNAgenetics, label = "DNA_N1200_04pop_migLow", use.wd=TRUE)
DNA_N1200_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                     genetics = DNAgenetics, label = "DNA_N1200_04pop_migHigh", use.wd=TRUE)
# Run parameter files. 
# For DNA simulations, set all.sites=TRUE (in order to print monomorphic as well as polymorphic sequence sites)
print("DNA: 4 populations, low migration")
DNA_N1200_04pop_migLow.params <- fscRun(DNA_N1200_04pop_migLow.params, num.sims = num_reps, 
                                  all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_N1200_04pop_migHigh.params <- fscRun(DNA_N1200_04pop_migHigh.params, num.sims = num_reps, 
                                   all.sites = TRUE, exec = fscVersion)
# Move .log and .par files into respective simulation folders
file.copy(from=paste0(dna.wd,"DNA_N1200_04pop_migLow.log"), to=paste0(dna.wd,"DNA_N1200_04pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_04pop_migLow.par"), to=paste0(dna.wd,"DNA_N1200_04pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_04pop_migHigh.log"), to=paste0(dna.wd,"DNA_N1200_04pop_migHigh/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_04pop_migHigh.par"), to=paste0(dna.wd,"DNA_N1200_04pop_migHigh/"))
file.remove(list.files(pattern = ".log"))
file.remove(list.files(pattern = ".par"))
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_N1200_04pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N1200_04pop_migLow/"),
                                              params = DNA_N1200_04pop_migLow.params)
DNA_N1200_04pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N1200_04pop_migHigh/"),
                                               params = DNA_N1200_04pop_migHigh.params)

# Save genind and params objects to Rdata files, for long term storage. Note that genind objects 
# need to have a file name of genind.scenarioName in order to be processed by the readGendinds_* functions
saveRDS(DNA_N1200_04pop_migLow.genind, file = paste0("data.DNA/genind.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_04pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_04pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_04pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
DNA_N1200_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16, 
                                    genetics = DNAgenetics, label = "DNA_N1200_16pop_migLow", use.wd=TRUE)
DNA_N1200_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16, 
                                     genetics = DNAgenetics, label = "DNA_N1200_16pop_migHigh", use.wd=TRUE)
# Run parameter files. 
# For DNA simulations, set all.sites=TRUE (in order to print monomorphic as well as polymorphic sequence sites)
print("DNA: 16 populations, low migration")
DNA_N1200_16pop_migLow.params <- fscRun(DNA_N1200_16pop_migLow.params, num.sims = num_reps, 
                                  all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_N1200_16pop_migHigh.params <- fscRun(DNA_N1200_16pop_migHigh.params, num.sims = num_reps, 
                                   all.sites = TRUE, exec = fscVersion)
# Move .log and .par files into respective simulation folders
file.copy(from=paste0(dna.wd,"DNA_N1200_16pop_migLow.log"), to=paste0(dna.wd,"DNA_N1200_16pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_16pop_migLow.par"), to=paste0(dna.wd,"DNA_N1200_16pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_16pop_migHigh.log"), to=paste0(dna.wd,"DNA_N1200_16pop_migHigh/"))
file.copy(from=paste0(dna.wd,"DNA_N1200_16pop_migHigh.par"), to=paste0(dna.wd,"DNA_N1200_16pop_migHigh/"))
file.remove(list.files(pattern = ".log"))
file.remove(list.files(pattern = ".par"))
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_N1200_16pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N1200_16pop_migLow/"), 
                                         params = DNA_N1200_16pop_migLow.params)
DNA_N1200_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N1200_16pop_migHigh/"), 
                                          params = DNA_N1200_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage. Note that genind objects 
# need to have a file name of genind.scenarioName in order to be processed by the readGendinds_* functions
saveRDS(DNA_N1200_16pop_migLow.genind, file = paste0("data.DNA/genind.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_16pop_migLow.params, file = paste0("data.DNA/params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_16pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_N1200_16pop_migHigh.params, file = paste0("data.DNA/params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))

# ---- N4800 ----
# SPECIFIC VARIABLES
# Specify number of total individuals, for all simulations
# Since there are 4 deme and 16 deme scenarios, this value must be divisible by 4 and 16
nInd <- 4800
# DEMES
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
# Specify working directory
dna.wd <- paste0(sim.wd,"SimulationOutputs/DNA_N4800_lowMut/")
setwd(dna.wd)

# 1 POPULATION
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_N4800_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, 
                                          label = "DNA_N4800_01pop_migLow", use.wd=TRUE)
DNA_N4800_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, 
                                           label = "DNA_N4800_01pop_migHigh", use.wd=TRUE)
# Run parameter files. 
# For DNA simulations, set all.sites=TRUE (in order to print monomorphic as well as polymorphic sequence sites)
print("DNA: 1 population, low migration")
DNA_N4800_01pop_migLow.params <- fscRun(DNA_N4800_01pop_migLow.params, num.sims = num_reps, 
                                        all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_N4800_01pop_migHigh.params <- fscRun(DNA_N4800_01pop_migHigh.params, num.sims = num_reps, 
                                         all.sites = TRUE, exec = fscVersion)
# Move .log and .par files into respective simulation folders
file.copy(from=paste0(dna.wd,"DNA_N4800_01pop_migLow.log"), to=paste0(dna.wd,"DNA_N4800_01pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_01pop_migLow.par"), to=paste0(dna.wd,"DNA_N4800_01pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_01pop_migHigh.log"), to=paste0(dna.wd,"DNA_N4800_01pop_migHigh/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_01pop_migHigh.par"), to=paste0(dna.wd,"DNA_N4800_01pop_migHigh/"))
file.remove(list.files(pattern = ".log"))
file.remove(list.files(pattern = ".par"))
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_N4800_01pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N4800_01pop_migLow/"), 
                                               params = DNA_N4800_01pop_migLow.params)
DNA_N4800_01pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N4800_01pop_migHigh/"), 
                                                params = DNA_N4800_01pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage. Note that genind objects 
# need to have a file name of genind.scenarioName in order to be processed by the readGendinds_* functions
saveRDS(DNA_N4800_01pop_migLow.genind, file = paste0("data.DNA/genind.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_01pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_01pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_01pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS
# Write parameter files
DNA_N4800_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4,
                                          genetics = DNAgenetics, label = "DNA_N4800_04pop_migLow", use.wd=TRUE)
DNA_N4800_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4,
                                           genetics = DNAgenetics, label = "DNA_N4800_04pop_migHigh", use.wd=TRUE)
# Run parameter files. 
# For DNA simulations, set all.sites=TRUE (in order to print monomorphic as well as polymorphic sequence sites)
print("DNA: 4 populations, low migration")
DNA_N4800_04pop_migLow.params <- fscRun(DNA_N4800_04pop_migLow.params, num.sims = num_reps, 
                                        all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_N4800_04pop_migHigh.params <- fscRun(DNA_N4800_04pop_migHigh.params, num.sims = num_reps, 
                                         all.sites = TRUE, exec = fscVersion)
# Move .log and .par files into respective simulation folders
file.copy(from=paste0(dna.wd,"DNA_N4800_04pop_migLow.log"), to=paste0(dna.wd,"DNA_N4800_04pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_04pop_migLow.par"), to=paste0(dna.wd,"DNA_N4800_04pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_04pop_migHigh.log"), to=paste0(dna.wd,"DNA_N4800_04pop_migHigh/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_04pop_migHigh.par"), to=paste0(dna.wd,"DNA_N4800_04pop_migHigh/"))
file.remove(list.files(pattern = ".log"))
file.remove(list.files(pattern = ".par"))
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_N4800_04pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N4800_04pop_migLow/"),
                                               params = DNA_N4800_04pop_migLow.params)
DNA_N4800_04pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N4800_04pop_migHigh/"),
                                                params = DNA_N4800_04pop_migHigh.params)

# Save genind and params objects to Rdata files, for long term storage. Note that genind objects 
# need to have a file name of genind.scenarioName in order to be processed by the readGendinds_* functions
saveRDS(DNA_N4800_04pop_migLow.genind, file = paste0("data.DNA/genind.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_04pop_migLow.params, file = paste0("data.DNA/params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_04pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_04pop_migHigh.params, file = paste0("data.DNA/params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS
# Write parameter files
DNA_N4800_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16, 
                                          genetics = DNAgenetics, label = "DNA_N4800_16pop_migLow", use.wd=TRUE)
DNA_N4800_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16, 
                                           genetics = DNAgenetics, label = "DNA_N4800_16pop_migHigh", use.wd=TRUE)
# Run parameter files. 
# For DNA simulations, set all.sites=TRUE (in order to print monomorphic as well as polymorphic sequence sites)
print("DNA: 16 populations, low migration")
DNA_N4800_16pop_migLow.params <- fscRun(DNA_N4800_16pop_migLow.params, num.sims = num_reps, 
                                        all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_N4800_16pop_migHigh.params <- fscRun(DNA_N4800_16pop_migHigh.params, num.sims = num_reps, 
                                         all.sites = TRUE, exec = fscVersion)
# Move .log and .par files into respective simulation folders
file.copy(from=paste0(dna.wd,"DNA_N4800_16pop_migLow.log"), to=paste0(dna.wd,"DNA_N4800_16pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_16pop_migLow.par"), to=paste0(dna.wd,"DNA_N4800_16pop_migLow/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_16pop_migHigh.log"), to=paste0(dna.wd,"DNA_N4800_16pop_migHigh/"))
file.copy(from=paste0(dna.wd,"DNA_N4800_16pop_migHigh.par"), to=paste0(dna.wd,"DNA_N4800_16pop_migHigh/"))
file.remove(list.files(pattern = ".log"))
file.remove(list.files(pattern = ".par"))
# Convert Arlequin outputs to genind
print("%%% Convert Arlequin outputs to genind")
DNA_N4800_16pop_migLow.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N4800_16pop_migLow/"), 
                                               params = DNA_N4800_16pop_migLow.params)
DNA_N4800_16pop_migHigh.genind <- convertAllArp(arp.path = paste0(dna.wd, "DNA_N4800_16pop_migHigh/"), 
                                                params = DNA_N4800_16pop_migHigh.params)
# Save genind and params objects to Rdata files, for long term storage. Note that genind objects 
# need to have a file name of genind.scenarioName in order to be processed by the readGendinds_* functions
saveRDS(DNA_N4800_16pop_migLow.genind, file = paste0("data.DNA/genind.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_16pop_migLow.params, file = paste0("data.DNA/params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_16pop_migHigh.genind, file = paste0("data.DNA/genind.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
saveRDS(DNA_N4800_16pop_migHigh.params, file = paste0("data.DNA/params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
