# %%%%%%%%%%%%%%%%%%%%%
# %%% LINEAR MODELS %%%
# %%%%%%%%%%%%%%%%%%%%%

# This script reads in resampling arrays generated for MSAT and SNP datasets (using nInd values
# of 1200 and 4800), and derives the 95% allelic diversity minimum sample size estimates (MSSEs)
# for those datasets. It also builds data.frames of the variable parameters used during the simulations
# (number of populations, migration rate, and marker type, MSAT or DNA).

# Using these dataframes, linear models are constructed to measure the relative effects of the simulation
# parameters (our explanatory variables) on MSSEs (our response variables). 

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
sim.wd <- "/home/akoontz/Shared/SSRvSNP_Sim/Code/"
setwd(sim.wd)
# Read in relevant functions
source("RScripts/functions_SSRvSNP_Sim.R")

# %%% NIND = 1200 %%% ----
# %%% Read in resampling arrays and calculate 95% minimum sampling size estimates (MSSEs)
# MSAT
MSAT_N1200_resampArr_filepath <- paste0(sim.wd, "SimulationOutputs/MSAT_marker/data.MSAT/MSAT_N1200_resampArr.Rdata")
MSAT_N1200_resamplingArrays <- readRDS(MSAT_N1200_resampArr_filepath)
MSAT_N1200_min95Values <- rapply(MSAT_N1200_resamplingArrays, resample_min95_mean, how = "list")

# DNA
DNA_N1200_resampArr_filepath <- paste0(sim.wd, "SimulationOutputs/DNA_marker/data.DNA/DNA_N1200_resampArr.Rdata")
DNA_N1200_resamplingArrays <-readRDS(DNA_N1200_resampArr_filepath)
DNA_N1200_min95Values <- rapply(DNA_N1200_resamplingArrays, resample_min95_mean, how = "list")

# %%% Build parameters data.frame, from which linear models will be built
params <- data.frame(expand.grid(n.pop=c(1,4,16), mig.Rate=c(0.001,0.01), marker=c("MSAT", "DNA")))
# Need to encode variables as text, versus integers, because we need R
# to recognize our explanatory variables as class factor (not numeric)
params <- data.frame(expand.grid(n.pop=c("one","four","sixteen"), mig.Rate=c("low","high"), marker=c("MSAT", "DNA")))
# Replicate each parameter combination 5 times, for each simulation replicate
results <- params[rep(1:nrow(params), each=5),]
# Append MSSE values to results data.frame
results$MSSE <- c(unlist(MSAT_N1200_min95Values), unlist(DNA_N1200_min95Values))

# %%% Checking Q Q plots (check links that Bela sent)

# %%% Run linear model
n1200_model <- lm(MSSE ~ (n.pop+mig.Rate+marker), data=results)
summary(n1200_model)
# Call ANOVA, to assess whether MSSEs are different for different explanatory variables
anova(n1200_model)
# Note that only the p values for marker type are less than 5%

# %%% Plotting this data
plot(MSSE ~ (n.pop+mig.Rate+marker), data=results)

# %%% NIND = 4800 %%% ----
# %%% Read in resampling arrays and calculate 95% minimum sampling size estimates (MSSEs)
# MSAT
MSAT_N4800_resampArr_filepath <- paste0(sim.wd, "SimulationOutputs/MSAT_N4800_marker/data.MSAT/MSAT_N4800_resampArr.Rdata")
MSAT_N4800_resamplingArrays <- readRDS(MSAT_N4800_resampArr_filepath)
MSAT_N4800_min95Values <- rapply(MSAT_N4800_resamplingArrays, resample_min95_mean, how = "list")

# DNA
DNA_N4800_resampArr_filepath <- paste0(sim.wd, "SimulationOutputs/DNA_N4800_marker/data.DNA/DNA_N4800_resampArr.Rdata")
DNA_N4800_resamplingArrays <-readRDS(DNA_N4800_resampArr_filepath)
DNA_N4800_min95Values <- rapply(DNA_N4800_resamplingArrays, resample_min95_mean, how = "list")

# %%% Build parameters data.frame, from which linear models will be built
params <- data.frame(expand.grid(n.pop=c(1,4,16), mig.Rate=c(0.001,0.01), marker=c("MSAT", "DNA")))
# Replicate each parameter combination 5 times, for each simulation replicate
results <- params[rep(1:nrow(params), each=5),]
# Append MSSE values to results data.frame
results$MSSE <- c(unlist(MSAT_N4800_min95Values), unlist(DNA_N4800_min95Values))

# %%% Run linear model
n4800_model <- lm(MSSE ~ (n.pop+mig.Rate+marker), data=results)
summary(n4800_model)
# Call ANOVA, to assess whether MSSEs are different for different explanatory variables
anova(n4800_model)
# In N4800 model, p values are lower for both marker type and migration rate

# %%% ACROSS TOTAL POPULATION SIZES %%% ----
# %%% Build parameters data.frame, from which linear models will be built
params <- data.frame(expand.grid(n.pop=c(1,4,16), mig.Rate=c(0.001,0.01), 
                                 t.pop.size=c(1200,4800), marker=c("MSAT", "DNA")))
# Replicate each parameter combination 5 times, for each simulation replicate
results <- params[rep(1:nrow(params), each=5),]
# Append MSSE values (from both N1200 and N4800 scenarios) to results data.frame
results$MSSE <- c(unlist(MSAT_N1200_min95Values), unlist(MSAT_N4800_min95Values),
                  unlist(DNA_N1200_min95Values), unlist(DNA_N4800_min95Values))

# %%% Run linear model
bothPopSizes_model <- lm(MSSE ~ (n.pop+mig.Rate+t.pop.size+marker), data=results)
summary(bothPopSizes_model)
# Call ANOVA, to assess whether MSSEs are different for different explanatory variables
anova(bothPopSizes_model)
# In model including both population sizes, p values are lower for both marker type and population size
