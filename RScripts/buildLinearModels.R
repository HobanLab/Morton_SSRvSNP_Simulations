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
# Specify numeric explanatory variables as categorical, using as.factor 
params <- data.frame(expand.grid(n.pop=as.factor(c(1,4,16)), mig.Rate=as.factor(c(0.001,0.01)), marker=c("MSAT", "DNA")))
# Replicate each parameter combination 5 times, for each simulation replicate
results.N1200 <- params[rep(1:nrow(params), each=5),]
# Append MSSE values to results data.frame
results.N1200$MSSE <- c(unlist(MSAT_N1200_min95Values), unlist(DNA_N1200_min95Values))

# %%% Run linear model
n1200_model <- lm(MSSE ~ (n.pop+mig.Rate+marker), data=results.N1200)
summary(n1200_model)
confint(n1200_model)
# Call ANOVA, to assess whether MSSEs are different for different explanatory variables
anova(n1200_model)
# Note that only the p values for marker type are less than 5%

# %%% Plotting model results
plot(MSSE ~ (n.pop+mig.Rate+marker), data=results.N1200, ylab="95% minimum sample size estimates", 
     main="MSSE values across explanatory variables (nInd = 1200)")
# Diagnostic model plots: used to assess whether the model meets our assumptions 
# (particularly, that model residuals are Normally distributed). Change mfrow (to allow multiple plots per window)
par(mfrow=c(2,2))
plot(n1200_model)
# Set mfrow back to default value (1 plot per window)
par(mfrow=c(1,1))

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
# Specify numeric explanatory variables as categorical, using as.factor 
params <- data.frame(expand.grid(n.pop=as.factor(c(1,4,16)), mig.Rate=as.factor(c(0.001,0.01)), marker=c("MSAT", "DNA")))
# Replicate each parameter combination 5 times, for each simulation replicate
results.N4800 <- params[rep(1:nrow(params), each=5),]
# Append MSSE values to results data.frame
results.N4800$MSSE <- c(unlist(MSAT_N4800_min95Values), unlist(DNA_N4800_min95Values))

# %%% Run linear model
n4800_model <- lm(MSSE ~ (n.pop+mig.Rate+marker), data=results.N4800)
summary(n4800_model)
# Call ANOVA, to assess whether MSSEs are different for different explanatory variables
anova(n4800_model)
confint(n4800_model)
# In N4800 model, p values are lower for both marker type and migration rate

# %%% Plotting model results
plot(MSSE ~ (n.pop+mig.Rate+marker), data=results.N4800, ylab="95% minimum sample size estimates", 
     main="MSSE values across explanatory variables (nInd = 4800)")
# Diagnostic model plots: used to assess whether the model meets our assumptions 
# (particularly, that model residuals are Normally distributed). Change mfrow (to allow multiple plots per window)
par(mfrow=c(2,2))
plot(n4800_model)
# Set mfrow back to default value (1 plot per window)
par(mfrow=c(1,1))

# %%% ACROSS TOTAL POPULATION SIZES %%% ----
# %%% Build parameters data.frame, from which linear models will be built
# Specify numeric explanatory variables as categorical, using as.factor 
params <- data.frame(expand.grid(n.pop=as.factor(c(1,4,16)), mig.Rate=as.factor(c(0.001,0.01)), 
                                 t.pop.size=as.factor(c(1200,4800)), marker=c("MSAT", "DNA")))
# Replicate each parameter combination 5 times, for each simulation replicate
results.total <- params[rep(1:nrow(params), each=5),]
# Append MSSE values (from both N1200 and N4800 scenarios) to results data.frame
results.total$MSSE <- c(unlist(MSAT_N1200_min95Values), unlist(MSAT_N4800_min95Values),
                  unlist(DNA_N1200_min95Values), unlist(DNA_N4800_min95Values))

# %%% Run linear model
bothPopSizes_model <- lm(MSSE ~ (n.pop+mig.Rate+t.pop.size+marker), data=results.total)
summary(bothPopSizes_model)
# Call ANOVA, to assess whether MSSEs are different for different explanatory variables
anova(bothPopSizes_model)
confint(bothPopSizes_model)
# In model including both population sizes, p values are lower for both marker type and population size

# %%% Plotting model results
plot(MSSE ~ (n.pop+mig.Rate+marker), data=results.total, ylab="95% minimum sample size estimates", 
     main="MSSE values across explanatory variables (both population sizes)")
# Diagnostic model plots: used to assess whether the model meets our assumptions 
# (particularly, that model residuals are Normally distributed). Change mfrow (to allow multiple plots per window)
par(mfrow=c(2,2))
plot(bothPopSizes_model)
# Set mfrow back to default value (1 plot per window)
par(mfrow=c(1,1))
