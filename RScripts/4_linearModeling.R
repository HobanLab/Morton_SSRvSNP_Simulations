# %%%%%%%%%%%%%%%%%%%%%
# %%% LINEAR MODELS %%%
# %%%%%%%%%%%%%%%%%%%%%

# This script reads in resampling arrays generated for MSAT and SNP datasets (using nInd values
# of 1200 and 4800; see 3_resampling.R), and builds dataframes of the parameter values used for 
# simulations (number of populations, total population size, migration rate, and marker type, 
# MSAT or DNA).

# Using these data.frames, 2 types of linear models are constructed: 
# 1. One which makes the number of samples the response variable. These models are used, along with
#    the predict function, to estimate the number of samples required, on average, to represent 95%
#    of the total allelic diversity found within a population (the 95% minimum sample size estimate, or MSSE)
# 2. One which sets the Total allelic representation as the response variable. These models are used 
#    to examine how different simulation parameters impact estimates of allelic representation.

# The script is split into 3 different sections based on total population size and mutation rate: 
# N1200 and N4800 (MSAT and DNA; for DNA, mutation rate of 1e-7), and then low mutation rate 
# DNA datasets (1e-8, N1200 and N4800). The results of linear models are also calculated for each
# level of loci used.

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
sim.wd <- '/home/akoontz/Shared/SSRvSNP_Sim/Code/'
setwd(sim.wd)
# Read in relevant functions
source('RScripts/0_functions.R')
# data.table::rbindlist: used below to collapse lists of data.frames into data.frames
library(data.table)
# Libraries for plotting
library(RColorBrewer) ; library(scales) ; library(report)

# %%% NIND = 1200 %%% ----
# %%% Read in resampling arrays and build results dataframes ----
# %%% MSAT
N1200_MSAT_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/MSAT_N1200_marker/data.MSAT/N0/MSAT_N1200_resampArr.Rdata')
N1200_MSAT_resamplingArrays <- readRDS(N1200_MSAT_resampArr_filepath)
# Transform resampling arrays into dataframes (for linear modeling)
N1200_MSAT_DF <- rapply(N1200_MSAT_resamplingArrays, resample_array2dataframe, how = 'list')
# For each item in N1200_MSAT_DF, transform the list of lists of dataframes such that results
# are grouped by the number of loci (the dataframe name). This is achieved using the 
# do.call(Map, c(rbind, X)) command.
N1200_MSAT_DF <- lapply(N1200_MSAT_DF, function(X) do.call(Map, c(rbind, X)))
# The resulting list of lists is separated by simulation scenario. Repeat the do.call command
# to generate a single level list organized by number of loci
N1200_MSAT_DF <- do.call(Map, c(rbind, N1200_MSAT_DF))
# %%% DNA
N1200_DNA_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/DNA_N1200_marker/data.DNA/N0/DNA_N1200_resampArr.Rdata')
N1200_DNA_resamplingArrays <-readRDS(N1200_DNA_resampArr_filepath)
# Transform resampling arrays into dataframes (for linear modeling)
N1200_DNA_DF <- rapply(N1200_DNA_resamplingArrays, resample_array2dataframe, how = 'list')
# Group results by number of loci, for each scenario
N1200_DNA_DF <- lapply(N1200_DNA_DF, function(X) do.call(Map, c(rbind, X)))
# Combine results from across scenarios
N1200_DNA_DF <- do.call(Map, c(rbind, N1200_DNA_DF))

# %%% Predict 95% MSSEs ----
# Declare a variable that's used as the prediction point (i.e. 95% "Total" allelic representation)
predictPoint <- data.frame(Total=95)
# %%% MSAT
# Create table to capture the MSSE ("fit"), lower PI, upper PI, and PI width, for both markers
N1200_MSAT_results <- matrix(data=NA, nrow=length(N1200_MSAT_DF), ncol=4)
colnames(N1200_MSAT_results) <- c('MSSE (fit)', 'Lower PI','Upper PI', 'PI Width')
rownames(N1200_MSAT_results) <- names(N1200_MSAT_DF)
# For each level of loci (rows in the results matrix)...
for(i in 1:nrow(N1200_MSAT_results)){
  # Build a linear model and predict
  N1200_MSAT_MSSEmodel <- lm(nSamples ~ I(Total^3), data=N1200_MSAT_DF[[i]])
  predictPoints <- predict(N1200_MSAT_MSSEmodel, predictPoint, interval = "prediction")
  # Calculate the PI Width (upper minus lower prediction interval)
  predictPoints <- c(predictPoints, predictPoints[[3]]-predictPoints[[2]])
  # Round the MSSE (fit) values to the next whole number of samples
  predictPoints[[1]] <- ceiling(predictPoints[[1]])
  # Store the results into the matrix
  N1200_MSAT_results[i,] <- predictPoints
  # For models using maximal number of loci, print a model summary for the Supplement
  if(i == nrow(N1200_MSAT_results)) {print(summary(N1200_MSAT_MSSEmodel))}
}
# %%% DNA
# Create table to capture the MSSE ("fit"), lower PI, upper PI, and PI width, for both markers
N1200_DNA_results <- matrix(data=NA, nrow=length(N1200_DNA_DF), ncol=4)
colnames(N1200_DNA_results) <- c('MSSE (fit)', 'Lower PI','Upper PI', 'PI Width')
rownames(N1200_DNA_results) <- names(N1200_DNA_DF)
# For each level of loci (rows in the results matrix)...
for(i in 1:nrow(N1200_DNA_results)){
  # Build a linear model and predict
  N1200_DNA_MSSEmodel <- lm(nSamples ~ I(Total^3), data=N1200_DNA_DF[[i]])
  predictPoints <- predict(N1200_DNA_MSSEmodel, predictPoint, interval = "prediction")
  # Calculate the PI Width (upper minus lower prediction interval)
  predictPoints <- c(predictPoints, predictPoints[[3]]-predictPoints[[2]])
  # Round the MSSE (fit) values to the next whole number of samples
  predictPoints[[1]] <- ceiling(predictPoints[[1]])
  # Store the results into the matrix
  N1200_DNA_results[i,] <- predictPoints
  # For models using maximal number of loci, print a model summary for the Supplement
  if(i == nrow(N1200_DNA_results)) {print(summary(N1200_DNA_MSSEmodel))}
}

# %%% Build linear models to capture the impact of simulation parameters ----
# Create dataframe combining MSAT and DNA results with highest loci level (last list item)
N1200_DF <- rbind(N1200_MSAT_DF[[length(N1200_MSAT_DF)]],N1200_DNA_DF[[length(N1200_DNA_DF)]])
# With markers
N1200_marker_model <- lm(I(Total^3) ~ nSamples + nPops + migRate + marker, data=N1200_DF)
summary(N1200_marker_model)
# Without markers
N1200_noMarker_model <- lm(I(Total^3) ~ nSamples + nPops + migRate, data=N1200_DF)
summary(N1200_noMarker_model)
# ANOVA demonstrating the value of introducing markers (note p value and large F statistic value)
anova(N1200_marker_model, N1200_noMarker_model)

# %%% Plotting model results ----
plotColors <- c(magma(n=12)[[8]],magma(n=12)[[11]])
# Use layout to generate multiple plots in a single window
layout(matrix(c(1,1,2,2,3,3), nrow=3, byrow=TRUE))
par(oma=rep(0.06,4), mar=c(1.5,4.5,3,1.5)+0.1)
# Standard resampling curve, MSAT data
plot(Total ~ nSamples, data=N1200_DF[which(N1200_DF$marker == "MSAT"),], 
     ylab='', xlab='', pch=16, col=alpha(plotColors[[1]], 0.01), 
     ylim=c(30,100.5), cex.axis=1.9)
# Standard resampling curve, DNA data
points(Total ~ nSamples, data=N1200_DF[which(N1200_DF$marker == "DNA"),],
       pch=16, col=alpha(plotColors[[2]], 0.01))
# Horizontal line indicating 95% coverage
abline(h=95, col="black", lty=3)
# Legend, title, and x-axis
par(cex=1.2)
legend(x=700, y=86, inset = 0.05, xpd=TRUE,
       legend = c('MSAT', 'DNA'), col=c(plotColors[[1]],plotColors[[2]]), 
       pch = c(20,20), pt.cex = 2, bty='n', y.intersp = 0.7)
par(cex=1.0)
mtext("N1200: Allelic representation across simulation parameters", line=0.7, cex=1.5)
mtext("Number of samples", side=1, line=2.4, cex=1)
# Boxplots
# Migration rate
par(cex.axis=1.2)
boxplot(Total ~ migRate+marker, data=N1200_DF, ylab='',
        xlab = '', col=c(rep(plotColors[[1]],2),rep(plotColors[[2]],2)),
        names=rep(c('Low migration (0.001)','High migration (0.01)'),2),
        medlwd=1, ylim=c(68,101))
par(cex.axis=1)
# Labels
abline(v=2.5)
mtext('MSAT', cex=1.2, side=3, line=0.08, at=1.5)
mtext('DNA', cex=1.2, side=3, line=0.08, at=3.5)
# Y axis label (all plots)
mtext(text="Allelic representation (%)", side=2, line=2.7, cex=1.2, srt=90)
# Increase bottom margin
par(mar=c(2,4.5,3,1.5)+0.1)
# Number of populations
par(cex.axis=1.2)
boxplot(Total ~ nPops+marker, data=N1200_DF, ylab='',
        xlab = '', col=c(rep(plotColors[[1]],3),rep(plotColors[[2]],3)),
        names=rep(c('1 population','4 populations','16 populations'),2),
        medlwd=1, ylim=c(68,101))
par(cex.axis=1)
# Labels
abline(v=3.5)
mtext('MSAT', cex=1.2, side=3, line=0.08, at=2)
mtext('DNA', cex=1.2, side=3, line=0.08, at=5)

# %%% NIND = 4800 %%% ----
# %%% Read in resampling arrays and build results dataframes ----
# %%% MSAT
N4800_MSAT_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/MSAT_N4800_marker/data.MSAT/N0/MSAT_N4800_resampArr.Rdata')
N4800_MSAT_resamplingArrays <- readRDS(N4800_MSAT_resampArr_filepath)
# Transform resampling arrays into dataframes (for linear modeling)
N4800_MSAT_DF <- rapply(N4800_MSAT_resamplingArrays, resample_array2dataframe, how = 'list')
# For each item in N4800_MSAT_DF, transform the list of lists of dataframes such that results
# are grouped by the number of loci (the dataframe name). This is achieved using the 
# do.call(Map, c(rbind, X)) command.
N4800_MSAT_DF <- lapply(N4800_MSAT_DF, function(X) do.call(Map, c(rbind, X)))
# The resulting list of lists is separated by simulation scenario. Repeat the do.call command
# to generate a single level list organized by number of loci
N4800_MSAT_DF <- do.call(Map, c(rbind, N4800_MSAT_DF))
# %%% DNA
N4800_DNA_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/DNA_N4800_marker/data.DNA/N0/DNA_N4800_resampArr.Rdata')
N4800_DNA_resamplingArrays <-readRDS(N4800_DNA_resampArr_filepath)
# Transform resampling arrays into dataframes (for linear modeling)
N4800_DNA_DF <- rapply(N4800_DNA_resamplingArrays, resample_array2dataframe, how = 'list')
# Group results by number of loci, for each scenario
N4800_DNA_DF <- lapply(N4800_DNA_DF, function(X) do.call(Map, c(rbind, X)))
# Combine results from across scenarios
N4800_DNA_DF <- do.call(Map, c(rbind, N4800_DNA_DF))

# %%% Predict 95% MSSEs ----
# Declare a variable that's used as the prediction point (i.e. 95% "Total" allelic representation)
predictPoint <- data.frame(Total=95)
# %%% MSAT
# Create table to capture the MSSE ("fit"), lower PI, upper PI, and PI width, for both markers
N4800_MSAT_results <- matrix(data=NA, nrow=length(N4800_MSAT_DF), ncol=4)
colnames(N4800_MSAT_results) <- c('MSSE (fit)', 'Lower PI','Upper PI', 'PI Width')
rownames(N4800_MSAT_results) <- names(N4800_MSAT_DF)
# For each level of loci (rows in the results matrix)...
for(i in 1:nrow(N4800_MSAT_results)){
  # Build a linear model and predict
  N4800_MSAT_MSSEmodel <- lm(nSamples ~ I(Total^3), data=N4800_MSAT_DF[[i]])
  predictPoints <- predict(N4800_MSAT_MSSEmodel, predictPoint, interval = "prediction")
  # Calculate the PI Width (upper minus lower prediction interval)
  predictPoints <- c(predictPoints, predictPoints[[3]]-predictPoints[[2]])
  # Round the MSSE (fit) values to the next whole number of samples
  predictPoints[[1]] <- ceiling(predictPoints[[1]])
  # Store the results into the matrix
  N4800_MSAT_results[i,] <- predictPoints
  # For models using maximal number of loci, print a model summary for the Supplement
  if(i == nrow(N4800_MSAT_results)) {print(summary(N4800_MSAT_MSSEmodel))}
}
# %%% DNA
# Create table to capture the MSSE ("fit"), lower PI, upper PI, and PI width, for both markers
N4800_DNA_results <- matrix(data=NA, nrow=length(N4800_DNA_DF), ncol=4)
colnames(N4800_DNA_results) <- c('MSSE (fit)', 'Lower PI','Upper PI', 'PI Width')
rownames(N4800_DNA_results) <- names(N4800_DNA_DF)
# For each level of loci (rows in the results matrix)...
for(i in 1:nrow(N4800_DNA_results)){
  # Build a linear model and predict
  N4800_DNA_MSSEmodel <- lm(nSamples ~ I(Total^3), data=N4800_DNA_DF[[i]])
  predictPoints <- predict(N4800_DNA_MSSEmodel, predictPoint, interval = "prediction")
  # Calculate the PI Width (upper minus lower prediction interval)
  predictPoints <- c(predictPoints, predictPoints[[3]]-predictPoints[[2]])
  # Round the MSSE (fit) values to the next whole number of samples
  predictPoints[[1]] <- ceiling(predictPoints[[1]])
  # Store the results into the matrix
  N4800_DNA_results[i,] <- predictPoints
  # For models using maximal number of loci, print a model summary for the Supplement
  if(i == nrow(N4800_DNA_results)) {print(summary(N4800_DNA_MSSEmodel))}
}

# %%% Build linear models to capture the impact of simulation parameters ----
# Create dataframe combining MSAT and DNA results with highest loci level (last list item)
N4800_DF <- rbind(N4800_MSAT_DF[[length(N4800_MSAT_DF)]],N4800_DNA_DF[[length(N4800_DNA_DF)]])
# With markers
N4800_marker_model <- lm(I(Total^3) ~ nSamples + nPops + migRate + marker, data=N4800_DF)
summary(N4800_marker_model)
# Without markers
N4800_noMarker_model <- lm(I(Total^3) ~ nSamples + nPops + migRate, data=N4800_DF)
summary(N4800_noMarker_model)
# ANOVA demonstrating the value of introducing markers (note p value and large F statistic value)
anova(N4800_marker_model, N4800_noMarker_model)

# %%% Plotting model results ----
plotColors <- c(magma(n=12)[[8]],magma(n=12)[[11]])
# Use layout to generate multiple plots in a single window
layout(matrix(c(1,1,2,2,3,3), nrow=3, byrow=TRUE))
par(oma=rep(0.06,4), mar=c(1.5,4.5,3,1.5)+0.1)
# Standard resampling curve, MSAT data
plot(Total ~ nSamples, data=N4800_DF[which(N4800_DF$marker == "MSAT"),], 
     ylab='', xlab='', pch=16, col=alpha(plotColors[[1]], 0.01), 
     ylim=c(26,100.5), cex.axis=1.9)
# Standard resampling curve, DNA data
points(Total ~ nSamples, data=N4800_DF[which(N4800_DF$marker == "DNA"),],
       pch=16, col=alpha(plotColors[[2]], 0.01))
# Horizontal line indicating 95% coverage
abline(h=95, col="black", lty=3)
# Legend, title, and x-axis
par(cex=1.2)
legend(x=3000, y=86, inset = 0.05, xpd=TRUE,
       legend = c('MSAT', 'DNA'), col=c(plotColors[[1]],plotColors[[2]]), 
       pch = c(20,20), pt.cex = 2, bty='n', y.intersp = 0.7)
par(cex=1.0)
mtext("N4800: Allelic representation across simulation parameters", line=0.7, cex=1.5)
mtext("Number of samples", side=1, line=2.4, cex=1)
# Boxplots
# Migration rate
par(cex.axis=1.2)
boxplot(Total ~ migRate+marker, data=N4800_DF, ylab='',
        xlab = '', col=c(rep(plotColors[[1]],2),rep(plotColors[[2]],2)),
        names=rep(c('Low migration (0.001)','High migration (0.01)'),2),
        medlwd=1, ylim=c(68,101))
par(cex.axis=1)
# Labels
abline(v=2.5)
mtext('MSAT', cex=1.2, side=3, line=0.08, at=1.5)
mtext('DNA', cex=1.2, side=3, line=0.08, at=3.5)
# Y axis label (all plots)
mtext(text="Allelic representation (%)", side=2, line=2.5, cex=1.2, srt=90)
# Increase bottom margin
par(mar=c(2,4.5,3,1.5)+0.1)
# Number of populations
par(cex.axis=1.2)
boxplot(Total ~ nPops+marker, data=N4800_DF, ylab='',
        xlab = '', col=c(rep(plotColors[[1]],3),rep(plotColors[[2]],3)),
        names=rep(c('1 population','4 populations','16 populations'),2),
        medlwd=1, ylim=c(68,101))
par(cex.axis=1)
# Labels
abline(v=3.5)
mtext('MSAT', cex=1.2, side=3, line=0.08, at=2)
mtext('DNA', cex=1.2, side=3, line=0.08, at=5)

# %%% DNA: LOW MUTATION RATE (1E-8) %%% ----
# %%% N1200 ----
# %%% Read in resampling arrays and build results data.frames
N1200_DNA_lowMut_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/DNA_N1200_lowMut/data.DNA/20231112/DNA_N1200_resampArr.Rdata')
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
N1200_DNA_lowMut_MSSEmodel <- lm(nSamples ~ I(Total^3), data=N1200_DNA_lowMut_DF)
predict(N1200_DNA_lowMut_MSSEmodel, predictPoint, interval = "prediction")
summary(N1200_DNA_lowMut_MSSEmodel)

# %%% Build linear models to capture the impact of simulation parameters
N1200_DNA_lowMut_model <- lm(I(Total^3) ~ nSamples + nPops + migRate, data=N1200_DNA_lowMut_DF)
summary(N1200_DNA_lowMut_model)
summary.aov(N1200_DNA_lowMut_model)

# %%% Plotting model results ----
# Use layout to generate multiple plots in a single window
layout(matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
# Standard resampling curve
plot(Total ~ nSamples, data=N1200_DNA_lowMut_DF, ylab='Allelic representation (%)',
     xlab='Number of samples', pch=16, col=alpha("darkblue", 0.10))
# Title
mtext("N1200: SNP: Low mutation (1e-8): Allelic representation across simulation parameters", line=2, cex=1.5)
# Boxplots for migration rate and number of populations
par(cex.axis=1.2)
boxplot(Total ~ migRate, data=N1200_DNA_lowMut_DF, ylab='Allelic representation (%)',
        xlab = 'Migration rate', pch=16)
boxplot(Total ~ nPops, data=N1200_DNA_lowMut_DF, ylab='Allelic representation (%)',
        xlab = 'Number of populations', pch=16)
par(cex.axis=1)

# %%% N4800 ----
# %%% Read in resampling arrays and build results data.frames
N4800_DNA_lowMut_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/DNA_N4800_lowMut/data.DNA//20231112/DNA_N4800_resampArr.Rdata')
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
N4800_DNA_lowMut_MSSEmodel <- lm(nSamples ~ I(Total^3), data=N4800_DNA_lowMut_DF)
predict(N4800_DNA_lowMut_MSSEmodel, predictPoint, interval = "prediction")
summary(N4800_DNA_lowMut_MSSEmodel)

# %%% Build linear models to capture the impact of simulation parameters
N4800_DNA_lowMut_model <- lm(I(Total^3) ~ nSamples + nPops + migRate, data=N4800_DNA_lowMut_DF)
summary(N4800_DNA_lowMut_model)
summary.aov(N4800_DNA_lowMut_model)

# %%% Plotting model results ----
# Use layout to generate multiple plots in a single window
layout(matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
# Standard resampling curve
plot(Total ~ nSamples, data=N4800_DNA_lowMut_DF, ylab='Allelic representation (%)',
     xlab='Number of samples', pch=16, col=alpha("darkblue", 0.10))
# Title
mtext("N4800: SNP: Low mutation (1e-8): Allelic representation across simulation parameters", line=2, cex=1.5)
# Boxplots for migration rate and number of populations
par(cex.axis=1.2)
boxplot(Total ~ migRate, data=N4800_DNA_lowMut_DF, ylab='Allelic representation (%)',
        xlab = 'Migration rate', pch=16)
boxplot(Total ~ nPops, data=N4800_DNA_lowMut_DF, ylab='Allelic representation (%)',
        xlab = 'Number of populations', pch=16)
par(cex.axis=1)
