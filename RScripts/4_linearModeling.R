# %%%%%%%%%%%%%%%%%%%%%
# %%% LINEAR MODELS %%%
# %%%%%%%%%%%%%%%%%%%%%

# This script reads in resampling arrays generated for MSAT and SNP datasets (using nInd values
# of 1200 and 4800; see 3_resampling.R), and builds data.frames of the parameter values used for 
# simulations (number of populations, total population size, migration rate, and marker type, MSAT or DNA).

# Using these data.frames, 2 types of linear models are constructed: 
# 1. One which makes the number of samples the response variable. These models are used, along with
#    the predict function, to estimate the number of samples required, on average, to represent 95%
#    of the total allelic diversity found within a population (the 95% minimum sample size estimate, or MSSE)
# 2. One which sets the Total allelic representation as the response variable. These models are used 
#    to examine how different simulation parameters impact estimates of allelic representation.

# The script is split into 3 different sections based on total population size and mutation rate: 
# N1200 and N4800 (MSAT and DNA; for DNA, mutation rate of 1e-7), and then low mutation rate 
# DNA datasets (1e-8, N1200 and N4800)

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
sim.wd <- '/home/akoontz/Shared/SSRvSNP_Sim/Code/'
setwd(sim.wd)
# Read in relevant functions
source('RScripts/0_functions.R')
# data.table::rbindlist: used below to collapse lists of data.frames into data.frames
library(data.table)
# Libraries for plotting
library(RColorBrewer) ; library(scales)

# %%% NIND = 1200 %%% ----
# %%% Read in resampling arrays and build results data.frames ----
# %%% MSAT
N1200_MSAT_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/MSAT_N1200_marker/data.MSAT/N0/MSAT_N1200_resampArr.Rdata')
N1200_MSAT_resamplingArrays <- readRDS(N1200_MSAT_resampArr_filepath)
# Transform resampling arrays into data.frames (for linear modeling)
N1200_MSAT_DF <- rapply(N1200_MSAT_resamplingArrays, resample_array2dataframe, how = 'list')
# Collapse multiple data.frames into single large data.frame (with all replicates of all scenarios;
# this step takes multiple calls because we have a list of lists, and rbindlist cannot be rapply-ed)
N1200_MSAT_DF <- lapply(N1200_MSAT_DF, rbindlist)
N1200_MSAT_DF <- rbindlist(N1200_MSAT_DF)
# %%% DNA
N1200_DNA_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/DNA_N1200_marker/data.DNA/N0/DNA_N1200_resampArr.Rdata')
N1200_DNA_resamplingArrays <-readRDS(N1200_DNA_resampArr_filepath)
# Transform resampling arrays into data.frames (for linear modeling)
N1200_DNA_DF <- rapply(N1200_DNA_resamplingArrays, resample_array2dataframe, how = 'list')
# Collapse multiple data.frames into single large data.frame (with all replicates of all scenarios;
# this step takes multiple calls because we have a list of lists, and rbindlist cannot be rapply-ed)
N1200_DNA_DF <- lapply(N1200_DNA_DF, rbindlist)
N1200_DNA_DF <- rbindlist(N1200_DNA_DF)
# %%% Combine DNA and MSAT results into single data.frame
N1200_DF <- rbind(N1200_MSAT_DF, N1200_DNA_DF)

# %%% Predict 95% MSSEs ----
# Declare a variable that's used as the prediction point (i.e. 95% "Total" allelic representation)
predictPoint <- data.frame(Total=95)
# MSAT: build model and predict
N1200_MSAT_MSSEmodel <- lm(nSamples ~ Total, data=N1200_MSAT_DF)
predict(N1200_MSAT_MSSEmodel, predictPoint, interval = "prediction")
# DNA: build model and predict
N1200_DNA_MSSEmodel <- lm(nSamples ~ Total, data=N1200_DNA_DF)
predict(N1200_DNA_MSSEmodel, predictPoint, interval = "prediction")

# %%% Build linear models to capture the impact of simulation parameters ----
# With markers
N1200_model <- lm(Total ~ nSamples + nPops + migRate + marker, data=N1200_DF)
summary(N1200_model)
summary.aov(N1200_model)
# Without markers
N1200_noMarker_model <- lm(Total ~ nSamples + nPops + migRate, data=N1200_DF)
summary(N1200_noMarker_model)
summary.aov(N1200_noMarker_model)
# ANOVA demonstrating the value of introducing markers (note p value and large F statistic value)
anova(N1200_model, N1200_noMarker_model)

# %%% Plotting model results ----
# Use layout to generate multiple plots in a single window
layout(matrix(c(1,1,2,2,3,3), nrow=3, byrow=TRUE))
# Standard resampling curve, MSAT data
plot(Total ~ nSamples, data=N1200_DF[which(N1200_DF$marker == "MSAT"),], 
     ylab='', xlab='Number of samples', pch=16, col=alpha("purple", 0.01))
# Standard resampling curve, DNA data
points(Total ~ nSamples, data=N1200_DF[which(N1200_DF$marker == "DNA"),],
       pch=16, col=alpha("darkgreen", 0.01))
# Legend and title
legend(x=625, y=156, inset = 0.05, xpd=TRUE,
       legend = c('MSAT', 'DNA'), col=c('purple','darkgreen'), 
       pch = c(20,20), cex=1.2, pt.cex = 2, bty='n', y.intersp = 0.15)
mtext("N1200: Allelic representation across scenario parameters", line=2, cex=1.5)

# Boxplots
# Migration rate
boxplot(Total ~ marker+migRate, data=N1200_DF, ylab='Allelic representation (%)', 
        xlab = '', col=c('darkgreen','purple','darkgreen','purple'),
        names=c('DNA','MSAT','DNA','MSAT'))
# Line below doesn't work...seems like it's tricky to add lines when you're using a custom plotting window
lines(2.5, 60, col = "red")
# Number of populations
boxplot(Total ~ marker+nPops, data=N1200_DF, ylab='Allelic representation (%)', 
        xlab = '', col=c('darkgreen','purple','darkgreen','purple','darkgreen','purple'),
        names=c('DNA','MSAT','DNA','MSAT','DNA','MSAT'))

# Boxplots for number of populations
# MSAT and DNA
boxplot(Total ~ migRate, data=N1200_DF[which(N1200_DF$marker == "MSAT"),], 
        ylab='Allelic representation (%)', xlab = 'Migration rate', pch=16, col='purple')

boxplot(Total ~ nPops, data=N1200_MSAT_DF, ylab='Allelic representation (%)',
        xlab = 'Number of populations', pch=16)
# %%% DNA
# Standard resampling curve
plot(Total ~ nSamples, data=N1200_DNA_DF, ylab='Allelic representation (%)',
     xlab='Number of samples', pch=16, col=alpha("darkgreen", 0.10))
# Title
mtext("N1200: SNP: Allelic representation across scenario parameters", line=2, cex=1.5)
# Boxplots for migration rate and number of populations
boxplot(Total ~ migRate, data=N1200_DNA_DF, ylab='Allelic representation (%)',
        xlab = 'Migration rate', pch=16)
boxplot(Total ~ nPops, data=N1200_DNA_DF, ylab='Allelic representation (%)',
        xlab = 'Number of populations', pch=16)

# %%% NIND = 4800 %%% ----
# %%% Read in resampling arrays and build results data.frames ----
# %%% MSAT
N4800_MSAT_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/MSAT_N4800_marker/data.MSAT/N0/MSAT_N4800_resampArr.Rdata')
N4800_MSAT_resamplingArrays <- readRDS(N4800_MSAT_resampArr_filepath)
# Transform resampling arrays into data.frames (for linear modeling)
N4800_MSAT_DF <- rapply(N4800_MSAT_resamplingArrays, resample_array2dataframe, how = 'list')
# Collapse multiple data.frames into single large data.frame (with all replicates of all scenarios;
# this step takes multiple calls because we have a list of lists, and rbindlist cannot be rapply-ed)
N4800_MSAT_DF <- lapply(N4800_MSAT_DF, rbindlist)
N4800_MSAT_DF <- rbindlist(N4800_MSAT_DF)
# %%% DNA
N4800_DNA_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/DNA_N4800_marker/data.DNA/N0/DNA_N4800_resampArr.Rdata')
N4800_DNA_resamplingArrays <-readRDS(N4800_DNA_resampArr_filepath)
# Transform resampling arrays into data.frames (for linear modeling)
N4800_DNA_DF <- rapply(N4800_DNA_resamplingArrays, resample_array2dataframe, how = 'list')
# Collapse multiple data.frames into single large data.frame (with all replicates of all scenarios;
# this step takes multiple calls because we have a list of lists, and rbindlist cannot be rapply-ed)
N4800_DNA_DF <- lapply(N4800_DNA_DF, rbindlist)
N4800_DNA_DF <- rbindlist(N4800_DNA_DF)
# %%% Combine DNA and MSAT results into single data.frame
N4800_DF <- rbind(N4800_MSAT_DF, N4800_DNA_DF)

# %%% Predict 95% MSSEs ----
# Declare a variable that's used as the prediction point (i.e. 95% "Total" allelic representation)
predictPoint <- data.frame(Total=95)
# MSAT: build model and predict
N4800_MSAT_MSSEmodel <- lm(nSamples ~ Total, data=N4800_MSAT_DF)
predict(N4800_MSAT_MSSEmodel, predictPoint, interval = "prediction")
# DNA: build model and predict
N4800_DNA_MSSEmodel <- lm(nSamples ~ Total, data=N4800_DNA_DF)
predict(N4800_DNA_MSSEmodel, predictPoint, interval = "prediction")

# %%% Build linear models to capture the impact of simulation parameters ----
# With markers
N4800_model <- lm(Total ~ nSamples + nPops + migRate + marker, data=N4800_DF)
summary(N4800_model)
summary.aov(N4800_model)
# Without markers
N4800_noMarker_model <- lm(Total ~ nSamples + nPops + migRate, data=N4800_DF)
summary(N4800_noMarker_model)
summary.aov(N4800_noMarker_model)
# ANOVA demonstrating the value of introducing markers (note p value and large F statistic value)
anova(N4800_model, N4800_noMarker_model)

# # %%% Plotting model results ----
# # Use layout to generate multiple plots in a single window
# layout(matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
# # %%% MSAT
# # Standard resampling curve
# plot(Total ~ nSamples, data=N4800_MSAT_DF, ylab='Allelic representation (%)',
#      xlab='Number of samples', pch=16, col=alpha("purple", 0.10))
# # Title
# mtext("N4800: MSAT: Allelic representation across scenario parameters", line=2, cex=1.5)
# # Boxplots for migration rate and number of populations
# boxplot(Total ~ migRate, data=N4800_MSAT_DF, ylab='Allelic representation (%)',
#         xlab = 'Migration rate', pch=16)
# boxplot(Total ~ nPops, data=N4800_MSAT_DF, ylab='Allelic representation (%)',
#         xlab = 'Number of populations', pch=16)
# # %%% DNA
# # Standard resampling curve
# plot(Total ~ nSamples, data=N4800_DNA_DF, ylab='Allelic representation (%)',
#      xlab='Number of samples', pch=16, col=alpha("darkgreen", 0.10))
# # Title
# mtext("N4800: SNP: Allelic representation across scenario parameters", line=2, cex=1.5)
# # Boxplots for migration rate and number of populations
# boxplot(Total ~ migRate, data=N4800_DNA_DF, ylab='Allelic representation (%)',
#         xlab = 'Migration rate', pch=16)
# boxplot(Total ~ nPops, data=N4800_DNA_DF, ylab='Allelic representation (%)',
#         xlab = 'Number of populations', pch=16)

# %%% DNA: LOW MUTATION RATE (1E-8) %%% ----
# %%% N1200 ----
# %%% Read in resampling arrays and build results data.frames
N1200_DNA_lowMut_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/DNA_N1200_lowMut/data.DNA/DNA_N1200_resampArr.Rdata')
N1200_DNA_lowMut_resamplingArrays <-readRDS(N1200_DNA_lowMut_resampArr_filepath)
# Transform resampling arrays into data.frames (for linear modeling)
N1200_DNA_lowMut_DF <- rapply(N1200_DNA_lowMut_resamplingArrays, resample_array2dataframe, how = 'list')
# Collapse multiple data.frames into single large data.frame (with all replicates of all scenarios;
# this step takes multiple calls because we have a list of lists, and rbindlist cannot be rapply-ed)
N1200_DNA_lowMut_DF <- lapply(N1200_DNA_lowMut_DF, rbindlist)
N1200_DNA_lowMut_DF <- rbindlist(N1200_DNA_lowMut_DF)

# %%% Predict 95% MSSEs
# Declare a variable that's used as the prediction point (i.e. 95% "Total" allelic representation)
predictPoint <- data.frame(Total=95)
# Build model and predict
N1200_DNA_lowMut_MSSEmodel <- lm(nSamples ~ Total, data=N1200_DNA_lowMut_DF)
predict(N1200_DNA_lowMut_MSSEmodel, predictPoint, interval = "prediction")

# %%% Build linear models to capture the impact of simulation parameters
N1200_DNA_lowMut_model <- lm(Total ~ nSamples + nPops + migRate, data=N1200_DNA_lowMut_DF)
summary(N1200_DNA_lowMut_model)
summary.aov(N1200_DNA_lowMut_model)

# %%% Plotting model results ----
# Use layout to generate multiple plots in a single window
layout(matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
# Standard resampling curve
plot(Total ~ nSamples, data=N1200_DNA_lowMut_DF, ylab='Allelic representation (%)',
     xlab='Number of samples', pch=16, col=alpha("darkblue", 0.10))
# Title
mtext("N1200: SNP: Low mutation (1e-8): Allelic representation across scenario parameters", line=2, cex=1.5)
# Boxplots for migration rate and number of populations
boxplot(Total ~ migRate, data=N1200_DNA_lowMut_DF, ylab='Allelic representation (%)',
        xlab = 'Migration rate', pch=16)
boxplot(Total ~ nPops, data=N1200_DNA_lowMut_DF, ylab='Allelic representation (%)',
        xlab = 'Number of populations', pch=16)

# %%% N4800 ----
# %%% Read in resampling arrays and build results data.frames
N4800_DNA_lowMut_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/DNA_N4800_lowMut/data.DNA/DNA_N4800_resampArr.Rdata')
N4800_DNA_lowMut_resamplingArrays <-readRDS(N4800_DNA_lowMut_resampArr_filepath)
# Transform resampling arrays into data.frames (for linear modeling)
N4800_DNA_lowMut_DF <- rapply(N4800_DNA_lowMut_resamplingArrays, resample_array2dataframe, how = 'list')
# Collapse multiple data.frames into single large data.frame (with all replicates of all scenarios;
# this step takes multiple calls because we have a list of lists, and rbindlist cannot be rapply-ed)
N4800_DNA_lowMut_DF <- lapply(N4800_DNA_lowMut_DF, rbindlist)
N4800_DNA_lowMut_DF <- rbindlist(N4800_DNA_lowMut_DF)

# %%% Predict 95% MSSEs
# Declare a variable that's used as the prediction point (i.e. 95% "Total" allelic representation)
predictPoint <- data.frame(Total=95)
# Build model and predict
N4800_DNA_lowMut_MSSEmodel <- lm(nSamples ~ Total, data=N4800_DNA_lowMut_DF)
predict(N4800_DNA_lowMut_MSSEmodel, predictPoint, interval = "prediction")

# %%% Build linear models to capture the impact of simulation parameters
N4800_DNA_lowMut_model <- lm(Total ~ nSamples + nPops + migRate, data=N4800_DNA_lowMut_DF)
summary(N4800_DNA_lowMut_model)
summary.aov(N4800_DNA_lowMut_model)

# %%% Plotting model results ----
# Use layout to generate multiple plots in a single window
layout(matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
# Standard resampling curve
plot(Total ~ nSamples, data=N4800_DNA_lowMut_DF, ylab='Allelic representation (%)',
     xlab='Number of samples', pch=16, col=alpha("darkblue", 0.10))
# Title
mtext("N4800: SNP: Low mutation (1e-8): Allelic representation across scenario parameters", line=2, cex=1.5)
# Boxplots for migration rate and number of populations
boxplot(Total ~ migRate, data=N4800_DNA_lowMut_DF, ylab='Allelic representation (%)',
        xlab = 'Migration rate', pch=16)
boxplot(Total ~ nPops, data=N4800_DNA_lowMut_DF, ylab='Allelic representation (%)',
        xlab = 'Number of populations', pch=16)
