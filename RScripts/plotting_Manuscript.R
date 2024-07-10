# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% MSAT VERSUS SNP SIMULATION STUDY: PLOTTING %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script generates all of the plots included in the manuscript describing the results of 
# simulation study comparing minimum sample size estimates (MSSE) and their variance using
# microsatellite and SNP markers. 

# It's divided into 2 sections: one building the figures used in the main text of the manuscript, 
# and one building the figures used in the supplement. Figures are wrapped in the PNG command to 
# automatically save them to the specified directory (plotDir)

library(strataG)
library(adegenet)
library(stringr)
library(parallel)
library(RColorBrewer)
library(scales)

# %%% FUNCTIONS AND VARIABLES %%% ----
# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
sim.wd <- "/home/akoontz/Shared/SSRvSNP_Sim/Code/"
setwd(sim.wd)
# Read in relevant functions
source("RScripts/functions_SSRvSNP_Sim.R")
# Specify the filepath to the folder to save all of the below graphs to
plotDir <- "/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/MS_1stDraft_20240710/"

# %%%% MAIN TEXT %%%% -----
# %%% ORIGINAL SIMULATIONS: NIND 1200 %%% ----
# Pick plot colors (for all plots!)
plotColors <- c("red","red4","darkorange3","coral","purple")

# N_TO_DROP = 0 %%% ----
# SPECIFY DIRECTORY TO SAVE RESAMPLING PLOTS
N1200_MSAT_N0_plotDir <- 
  "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N1200/MSAT/N0/"
N1200_DNA_N0_plotDir <- 
  "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N1200/DNA/N0/"
# READ IN RESAMPLING ARRAYS
# MSAT
N1200_MSAT_N0_resampArr_filepath <- 
  paste0(sim.wd, "SimulationOutputs/MSAT_N1200_marker/data.MSAT/N0/MSAT_N1200_resampArr.Rdata")
N1200_MSAT_N0_resamplingArrays <- readRDS(file=N1200_MSAT_N0_resampArr_filepath)
# DNA
N1200_DNA_N0_resampArr_filepath <- 
  paste0(sim.wd, "SimulationOutputs/DNA_N1200_marker/data.DNA/N0/DNA_N1200_resampArr.Rdata")
N1200_DNA_N0_resamplingArrays <- readRDS(file=N1200_DNA_N0_resampArr_filepath)
# PLOT RESAMPLING RESULTS
# MSAT
invisible(rapply(N1200_MSAT_N0_resamplingArrays, resample_Plot_PNG, 
                 colors=plotColors, data.dir=N1200_MSAT_N0_plotDir))
# DNA
invisible(rapply(N1200_DNA_N0_resamplingArrays, resample_Plot_PNG, 
                 colors=plotColors, data.dir=N1200_DNA_N0_plotDir))

# %%% ORIGINAL SIMULATIONS: NIND 4800 %%% ----
# N_TO_DROP = 0 %%% ----
# SPECIFY DIRECTORY TO SAVE RESAMPLING PLOTS
N4800_MSAT_N0_plotDir <- 
  "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N4800/MSAT/N0/"
N4800_DNA_N0_plotDir <- 
  "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N4800/DNA/N0/"
# READ IN RESAMPLING ARRAYS
# MSAT
N4800_MSAT_N0_resampArr_filepath <- 
  paste0(sim.wd, "SimulationOutputs/MSAT_N4800_marker/data.MSAT/N0/MSAT_N4800_resampArr.Rdata")
N4800_MSAT_N0_resamplingArrays <- readRDS(file=N4800_MSAT_N0_resampArr_filepath)
# DNA
N4800_DNA_N0_resampArr_filepath <- 
  paste0(sim.wd, "SimulationOutputs/DNA_N4800_marker/data.DNA/N0/DNA_N4800_resampArr.Rdata")
N4800_DNA_N0_resamplingArrays <- readRDS(file=N4800_DNA_N0_resampArr_filepath)
# PLOT RESAMPLING RESULTS
# MSAT
invisible(rapply(N4800_MSAT_N0_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE, 
                 colors=plotColors, data.dir=N4800_MSAT_N0_plotDir))
# DNA
invisible(rapply(N4800_DNA_N0_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE, 
                 colors=plotColors, data.dir=N4800_DNA_N0_plotDir))

# %%%% SUPPLEMENT %%%% -----
# %%% FUNCTIONS %%% ----
# Read in relevant functions from repository
source("RScripts/functions_SSRvSNP_Sim.R")

# Declare custom function for plotting a single histogram of allele frequencies from a genind object
makeAlleleFreqHist <- function(gen.obj, title=gen.obj@other, yMax=25){
  # Make a vector of allele frequency values, from the genind object
  wildAlleleFreqs <- getWildFreqs(gen.obj, wholeValues = FALSE)
  # Generate histogram
  hist(wildAlleleFreqs, main=title, freq=FALSE, breaks=seq(0, 1.0, 0.01), 
       ylim=c(0, yMax), xlab="", ylab="")
}

# %%% DNA LOW MUTATION RATES %%% ----
# SPECIFY DIRECTORY TO SAVE RESAMPLING PLOTS
N1200_DNA_lowMut_plotDir <- 
  "~/Documents/SSRvSNP/Simulations/Documentation/Images/dnaMutationRateTests_052023/N1200/"
N4800_DNA_lowMut_plotDir <- 
  "~/Documents/SSRvSNP/Simulations/Documentation/Images/dnaMutationRateTests_052023/N4800/"
# READ IN RESAMPLING ARRAYS
# N1200
N1200_DNA_lowMut_resampArr_filepath <- 
  paste0(sim.wd, "SimulationOutputs/DNA_N1200_lowMut/data.DNA/DNA_N1200_resampArr.Rdata")
N1200_DNA_lowMut_resamplingArrays <- readRDS(file=N1200_DNA_lowMut_resampArr_filepath)
# N4800
N4800_DNA_lowMut_resampArr_filepath <- 
  paste0(sim.wd, "SimulationOutputs/DNA_N4800_lowMut/data.DNA/DNA_N4800_resampArr.Rdata")
N4800_DNA_lowMut_resamplingArrays <- readRDS(file=N4800_DNA_lowMut_resampArr_filepath)
# PLOT RESAMPLING RESULTS
# N1200
invisible(rapply(N1200_DNA_lowMut_resamplingArrays, resample_Plot_PNG, 
                 colors=plotColors, data.dir=N1200_DNA_lowMut_plotDir))
# N4800
invisible(rapply(N4800_DNA_lowMut_resamplingArrays, resample_Plot_PNG, 
                 colors=plotColors, largePopFlag=TRUE, data.dir=N4800_DNA_lowMut_plotDir))


# %%% READ IN DATA %%% ----
# %%% SIMULATED
# Source the genind objects from previously run simulations, using readGeninds functions
# MSAT
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N1200_marker/data.MSAT/"), prefix="MSAT_N1200")
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N4800_marker/data.MSAT/"), prefix="MSAT_N4800")
# DNA
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N1200_marker/data.DNA/"), prefix="DNA_N1200")
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N4800_marker/data.DNA/"), prefix="DNA_N4800")
# DNA low mutation (1e-8)
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N1200_lowMut/data.DNA/"), prefix="DNA_N1200_lowMut")
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N4800_lowMut/data.DNA/"), prefix="DNA_N4800_lowMut")

# %%% EMPIRICAL (QUAC DATASETS)
# MSAT 
# Declare filepatth to gen object
QUAC.MSAT.filePath <- 
  "/home/akoontz/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/"
# QUAC.MSAT.genind <- 
#   read.genepop(paste0(QUAC.MSAT.filePath, "Adegenet_Files/QUAC_woK_allpop_clean.gen"), ncode = 3)
QUAC.MSAT.genind <- 
  read.genepop(paste0(QUAC.MSAT.filePath, "Adegenet_Files/QUAC_wK_allpop_clean.gen"), ncode = 3)
# Correct popNames: samples with popName pattern QAc-G- are garden 
levels(QUAC.MSAT.genind@pop)[grep(pattern = "QAc-G-", levels(QUAC.MSAT.genind@pop))] <- 
  rep("garden", length(grep(pattern = "QAc-G-", levels(QUAC.MSAT.genind@pop))))
# Correct popNames: samples with popName pattern QAc-W- are wild
levels(QUAC.MSAT.genind@pop)[grep(pattern = "QAc-W-", levels(QUAC.MSAT.genind@pop))] <- 
  rep("wild", length(grep(pattern = "QAc-W-", levels(QUAC.MSAT.genind@pop))))
# SNP
# Declare filepatth to gen object
# genpop.filePath <- 
#   "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops_NoK/"
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
QUAC.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUAC_popmap_GardenWild"), header=FALSE)[,2])

# %%% PLOTTING %%% ----
# %%% MSAT
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# SIMULATED (N1200)
makeAlleleFreqHist(MSAT_N1200_01pop_migHigh.genList[[3]], title = "")
mtext("Simulated (N1200; 4 populations, high migration)", side=1, line=-7, cex=1.2)
# Plot title
title("Microsatellites: Allele frequency distributions", line=1.5, cex.main=2)
# Update margins to reduce space above middle graph
par(mar=c(2,4.5,2,1.5)+0.1)
# SIMULATED (N4800)
makeAlleleFreqHist(MSAT_N4800_01pop_migHigh.genList[[3]], title = "", yMax = 20)
mtext("Simulated (N4800; 4 populations, high migration)", side=1, line=-7, cex=1.2)
# y-axis label
mtext(text="Number of alleles", side=2, line=3, cex=1, srt=90)
# Update margins to allow for more space below bottom graph
par(mar=c(4,4.5,2,1.5)+0.1)
# Empirical
makeAlleleFreqHist(QUAC.MSAT.genind, title = "", yMax = 45)
mtext("Empirical comparison (Quercus acerifolia)", side=1, line=-7, cex=1.2)
# x-axis label
mtext(text="Allele frequency category", side=1, line=2.5, cex=1)

# %%% SNP
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# SIMULATED (N1200)
makeAlleleFreqHist(DNA_N1200_01pop_migHigh.genList[[3]], title = "", yMax = 35)
mtext("Simulated (N1200; 4 populations, high migration)", side=1, line=-7, cex=1.2)
# Plot title
title("SNPs: Allele frequency distributions", line=1.5, cex.main=2)
# Update margins to reduce space above middle graph
par(mar=c(2,4.5,2,1.5)+0.1)
# SIMULATED (N4800)
makeAlleleFreqHist(DNA_N4800_01pop_migHigh.genList[[3]], title = "", yMax = 60)
mtext("Simulated (N4800; 4 populations, high migration)", side=1, line=-7, cex=1.2)
# y-axis label
mtext(text="Number of alleles", side=2, line=3, cex=1, srt=90)
# Update margins to allow for more space below bottom graph
par(mar=c(4,4.5,2,1.5)+0.1)
# Empirical
makeAlleleFreqHist(QUAC.SNP.DN.R80.genind, title = "", yMax = 30)
mtext("Empirical comparison (Quercus acerifolia)", side=1, line=-7, cex=1.2)
# x-axis label
mtext(text="Allele frequency category", side=1, line=2.5, cex=1)

# %%% SNP (Low mutation)
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# SIMULATED (N1200)
makeAlleleFreqHist(DNA_N1200_lowMut_01pop_migHigh.genList[[3]], title = "", yMax = 90)
mtext("Simulated (N1200, low mutation; 4 populations, high migration)", side=1, line=-7, cex=1.2)
# Plot title
title("SNPs (low mutation): Allele frequency distributions", line=1.5, cex.main=2)
# Update margins to reduce space above middle graph
par(mar=c(2,4.5,2,1.5)+0.1)
# SIMULATED (N4800)
makeAlleleFreqHist(DNA_N4800_lowMut_01pop_migHigh.genList[[3]], title = "", yMax = 40)
mtext("Simulated (N4800, low mutation; 4 populations, high migration)", side=1, line=-7, cex=1.2)
# y-axis label
mtext(text="Number of alleles", side=2, line=3, cex=1, srt=90)
# Update margins to allow for more space below bottom graph
par(mar=c(4,4.5,2,1.5)+0.1)
# Empirical
makeAlleleFreqHist(QUAC.SNP.DN.R80.genind, title = "", yMax=30)
mtext("Empirical comparison (Quercus acerifolia)", side=1, line=-7, cex=1.2)
# x-axis label
mtext(text="Allele frequency category", side=1, line=2.5, cex=1)
