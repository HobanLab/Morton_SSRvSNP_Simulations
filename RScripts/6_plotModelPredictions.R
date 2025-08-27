# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% MSSE PREDICTION MODEL PLOTTING %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script plots the raw data from the MSATvSNP simulation studies (4 datasets total),
# along with the lines describing the models generated to make predictions of 95% MSSEs.

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
sim.wd <- '/home/akoontz/Shared/SSRvSNP_Sim/Code/'
setwd(sim.wd)
# Read in relevant functions
source('RScripts/0_functions.R')
# data.table::rbindlist: used below to collapse lists of data.frames into data.frames
library(data.table)
# Libraries for plotting
library(RColorBrewer) ; library(scales) ; library(report) ; library(viridis)
# Specify directory to save graph outputs to
imageOutDir <- 
  '/home/akoontz/Documents//SSRvSNP_Simulations/Documentation/Images/MS_Revisions_20250409/'

# %%% NIND = 1200 %%% ----
# %%% Read in resampling arrays and build results dataframes ----
# %%% MSAT
N1200_MSAT_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/MSAT_N1200_marker/data.MSAT/N0/MSAT_N1200_resampArr.Rdata')
N1200_MSAT_resamplingArrays <- readRDS(N1200_MSAT_resampArr_filepath)
# Transform resampling arrays into dataframes (for modeling)
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
# Transform resampling arrays into dataframes (for modeling)
N1200_DNA_DF <- rapply(N1200_DNA_resamplingArrays, resample_array2dataframe, how = 'list')
# Group results by number of loci, for each scenario
N1200_DNA_DF <- lapply(N1200_DNA_DF, function(X) do.call(Map, c(rbind, X)))
# Combine results from across scenarios
N1200_DNA_DF <- do.call(Map, c(rbind, N1200_DNA_DF))
# Create dataframe combining MSAT and DNA results with highest loci level (last list item)
N1200_DF <- rbind(N1200_MSAT_DF[[length(N1200_MSAT_DF)]],N1200_DNA_DF[[length(N1200_DNA_DF)]])

# %%% Predict 95% MSSEs ----
# %%% MSAT
# Create a sequence of x values for smooth prediction line, and predict
N1200_MSAT_X <- seq(min(N1200_DF$Total), max(N1200_DF$Total), length.out = 200)
# CUBIC
# Build a model, using the largest number of loci
N1200_MSAT_cubModel <- lm(nSamples ~ poly(Total, 3, raw=TRUE), data=N1200_MSAT_DF[[5]])
N1200_MSAT_cubPredict <- predict(N1200_MSAT_cubModel, newdata = data.frame(Total=N1200_MSAT_X),
                                 interval = "prediction", level = 0.95)
# EXPONENTIAL
# Build a model
N1200_MSAT_expModel <- lm(log(nSamples) ~ Total, data=N1200_MSAT_DF[[5]])
N1200_MSAT_expPredict_log <- predict(N1200_MSAT_expModel, newdata = data.frame(Total=N1200_MSAT_X),
                                     interval = "prediction", level = 0.95)
# Back-transform predictions
N1200_MSAT_expPredict <- exp(N1200_MSAT_expPredict_log)

# %%% DNA
# Create a sequence of x values for smooth prediction line, and predict
N1200_DNA_X <- seq(min(N1200_DF$Total), max(N1200_DF$Total), length.out = 200)
# CUBIC
# Build a model, using the largest number of loci
N1200_DNA_cubModel <- lm(nSamples ~ poly(Total, 3, raw=TRUE), data=N1200_DNA_DF[[5]])
N1200_DNA_cubPredict <- predict(N1200_DNA_cubModel, newdata = data.frame(Total=N1200_DNA_X),
                                 interval = "prediction", level = 0.95)
# EXPONENTIAL
# Build a model
N1200_DNA_expModel <- lm(log(nSamples) ~ Total, data=N1200_DNA_DF[[5]])
N1200_DNA_expPredict_log <- predict(N1200_DNA_expModel, newdata = data.frame(Total=N1200_DNA_X),
                                     interval = "prediction", level = 0.95)
# Back-transform predictions
N1200_DNA_expPredict <- exp(N1200_DNA_expPredict_log)

# %%% Plotting models ----
# Specify colors, and stack 2 plots vertically
plotColors <- c(magma(n=12)[[8]], magma(n=12)[[11]])
par(mfcol=c(2,1), mar=c(2,3,4,2)+0.1)
# %%% MSAT
# CUBIC
# Plot samples versus allelic representation
plot(nSamples ~ Total, data=N1200_DF[which(N1200_DF$marker == "MSAT"),], 
     ylab='Number of Samples', xlab='Allelic Representation (%)', 
     pch=16, col=alpha(plotColors[[1]], 0.01), 
     xlim=c(30,100.5), ylim=c(0,1300), cex.axis=1.2, 
     main='MSAT 1200 MSSE Prediction Model: Cubic')
# Add fitted prediction line and intervals
lines(N1200_MSAT_X, N1200_MSAT_cubPredict[,'fit'], col = 'orange', lwd = 2)
lines(N1200_MSAT_X, N1200_MSAT_cubPredict[,'lwr'], col = 'red', lty = 2)
lines(N1200_MSAT_X, N1200_MSAT_cubPredict[,'upr'], col = 'red', lty = 2)
# Add legend
legend(x=30, y=1300, legend = c("Fitted line", "95% PI"), 
       col = c("blue","red"), lwd = 2, lty = c(1,2), bty = "n")
# EXPONENTIAL
# Plot samples versus allelic representation
plot(nSamples ~ Total, data=N1200_DF[which(N1200_DF$marker == "MSAT"),], 
     ylab='Number of Samples', xlab='Allelic Representation (%)', 
     pch=16, col=alpha(plotColors[[1]], 0.01), 
     xlim=c(30,100.5), ylim=c(0,1300), cex.axis=1.2, 
     main='MSAT 1200 MSSE Prediction Model: Exponential')
# Add fitted prediction line and intervals
lines(N1200_MSAT_X, N1200_MSAT_expPredict[,'fit'], col = 'blue', lwd = 2)
lines(N1200_MSAT_X, N1200_MSAT_expPredict[,'lwr'], col = 'purple', lty = 2)
lines(N1200_MSAT_X, N1200_MSAT_expPredict[,'upr'], col = 'purple', lty = 2)
# Add legend
legend(x=30, y=1300, legend = c("Fitted line", "95% PI"), 
       col = c("blue","red"), lwd = 2, lty = c(1,2), bty = "n")

# %%% DNA
# CUBIC
# Plot samples versus allelic representation
plot(nSamples ~ Total, data=N1200_DF[which(N1200_DF$marker == "DNA"),], 
     ylab='Number of Samples', xlab='Allelic Representation (%)', 
     pch=16, col=alpha(plotColors[[1]], 0.01), 
     xlim=c(30,100.5), ylim=c(0,1300), cex.axis=1.2, 
     main='DNA 1200 MSSE Prediction Model: Cubic')
# Add fitted prediction line and intervals
lines(N1200_DNA_X, N1200_DNA_cubPredict[,'fit'], col = 'orange', lwd = 2)
lines(N1200_DNA_X, N1200_DNA_cubPredict[,'lwr'], col = 'red', lty = 2)
lines(N1200_DNA_X, N1200_DNA_cubPredict[,'upr'], col = 'red', lty = 2)
# Add legend
legend(x=30, y=1300, legend = c("Fitted line", "95% PI"), 
       col = c("blue","red"), lwd = 2, lty = c(1,2), bty = "n")
# EXPONENTIAL
# Plot samples versus allelic representation
plot(nSamples ~ Total, data=N1200_DF[which(N1200_DF$marker == "DNA"),], 
     ylab='Number of Samples', xlab='Allelic Representation (%)', 
     pch=16, col=alpha(plotColors[[1]], 0.01), 
     xlim=c(30,100.5), ylim=c(0,1300), cex.axis=1.2, 
     main='DNA 1200 MSSE Prediction Model: Exponential')
# Add fitted prediction line and intervals
lines(N1200_DNA_X, N1200_DNA_expPredict[,'fit'], col = 'blue', lwd = 2)
lines(N1200_DNA_X, N1200_DNA_expPredict[,'lwr'], col = 'purple', lty = 2)
lines(N1200_DNA_X, N1200_DNA_expPredict[,'upr'], col = 'purple', lty = 2)
# Add legend
legend(x=30, y=1300, legend = c("Fitted line", "95% PI"), 
       col = c("blue","red"), lwd = 2, lty = c(1,2), bty = "n")

# %%% NIND = 4800 %%% ----
# %%% Read in resampling arrays and build results dataframes ----
# %%% MSAT
N4800_MSAT_resampArr_filepath <- 
  paste0(sim.wd, 'SimulationOutputs/MSAT_N4800_marker/data.MSAT/N0/MSAT_N4800_resampArr.Rdata')
N4800_MSAT_resamplingArrays <- readRDS(N4800_MSAT_resampArr_filepath)
# Transform resampling arrays into dataframes (for modeling)
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
# Transform resampling arrays into dataframes (for modeling)
N4800_DNA_DF <- rapply(N4800_DNA_resamplingArrays, resample_array2dataframe, how = 'list')
# Group results by number of loci, for each scenario
N4800_DNA_DF <- lapply(N4800_DNA_DF, function(X) do.call(Map, c(rbind, X)))
# Combine results from across scenarios
N4800_DNA_DF <- do.call(Map, c(rbind, N4800_DNA_DF))
# Create dataframe combining MSAT and DNA results with highest loci level (last list item)
N4800_DF <- rbind(N4800_MSAT_DF[[length(N4800_MSAT_DF)]],N4800_DNA_DF[[length(N4800_DNA_DF)]])

# %%% Predict 95% MSSEs ----
# %%% MSAT
# Create a sequence of x values for smooth prediction line, and predict
N4800_MSAT_X <- seq(min(N4800_DF$Total), max(N4800_DF$Total), length.out = 200)
# CUBIC
# Build a model, using the largest number of loci
N4800_MSAT_cubModel <- lm(nSamples ~ poly(Total, 3, raw=TRUE), data=N4800_MSAT_DF[[5]])
N4800_MSAT_cubPredict <- predict(N4800_MSAT_cubModel, newdata = data.frame(Total=N4800_MSAT_X),
                                 interval = "prediction", level = 0.95)
# EXPONENTIAL
# Build a model
N4800_MSAT_expModel <- lm(log(nSamples) ~ Total, data=N4800_MSAT_DF[[5]])
N4800_MSAT_expPredict_log <- predict(N4800_MSAT_expModel, newdata = data.frame(Total=N4800_MSAT_X),
                                     interval = "prediction", level = 0.95)
# Back-transform predictions
N4800_MSAT_expPredict <- exp(N4800_MSAT_expPredict_log)

# %%% DNA
# Create a sequence of x values for smooth prediction line, and predict
N4800_DNA_X <- seq(min(N4800_DF$Total), max(N4800_DF$Total), length.out = 200)
# CUBIC
# Build a model, using the largest number of loci
N4800_DNA_cubModel <- lm(nSamples ~ poly(Total, 3, raw=TRUE), data=N4800_DNA_DF[[5]])
N4800_DNA_cubPredict <- predict(N4800_DNA_cubModel, newdata = data.frame(Total=N4800_DNA_X),
                                interval = "prediction", level = 0.95)
# EXPONENTIAL
# Build a model
N4800_DNA_expModel <- lm(log(nSamples) ~ Total, data=N4800_DNA_DF[[5]])
N4800_DNA_expPredict_log <- predict(N4800_DNA_expModel, newdata = data.frame(Total=N4800_DNA_X),
                                    interval = "prediction", level = 0.95)
# Back-transform predictions
N4800_DNA_expPredict <- exp(N4800_DNA_expPredict_log)

# %%% Plotting models ----
# Specify colors, and stack 2 plots vertically
plotColors <- c(magma(n=12)[[8]], magma(n=12)[[11]])
par(mfcol=c(2,1), mar=c(2,3,4,2)+0.1)
# %%% MSAT
# CUBIC
# Plot samples versus allelic representation
plot(nSamples ~ Total, data=N4800_DF[which(N4800_DF$marker == "MSAT"),], 
     ylab='Number of Samples', xlab='Allelic Representation (%)', 
     pch=16, col=alpha(plotColors[[1]], 0.01), 
     xlim=c(30,100.5), ylim=c(0,4900), cex.axis=1.2, 
     main='MSAT 4800 MSSE Prediction Model: Cubic')
# Add fitted prediction line and intervals
lines(N4800_MSAT_X, N4800_MSAT_cubPredict[,'fit'], col = 'orange', lwd = 2)
lines(N4800_MSAT_X, N4800_MSAT_cubPredict[,'lwr'], col = 'red', lty = 2)
lines(N4800_MSAT_X, N4800_MSAT_cubPredict[,'upr'], col = 'red', lty = 2)
# Add legend
legend(x=30, y=4900, legend = c("Fitted line", "95% PI"), 
       col = c("blue","red"), lwd = 2, lty = c(1,2), bty = "n")
# EXPONENTIAL
# Plot samples versus allelic representation
plot(nSamples ~ Total, data=N4800_DF[which(N4800_DF$marker == "MSAT"),], 
     ylab='Number of Samples', xlab='Allelic Representation (%)', 
     pch=16, col=alpha(plotColors[[1]], 0.01), 
     xlim=c(30,100.5), ylim=c(0,4900), cex.axis=1.2, 
     main='MSAT 4800 MSSE Prediction Model: Exponential')
# Add fitted prediction line and intervals
lines(N4800_MSAT_X, N4800_MSAT_expPredict[,'fit'], col = 'blue', lwd = 2)
lines(N4800_MSAT_X, N4800_MSAT_expPredict[,'lwr'], col = 'purple', lty = 2)
lines(N4800_MSAT_X, N4800_MSAT_expPredict[,'upr'], col = 'purple', lty = 2)
# Add legend
legend(x=30, y=4900, legend = c("Fitted line", "95% PI"), 
       col = c("blue","red"), lwd = 2, lty = c(1,2), bty = "n")

# %%% DNA
# CUBIC
# Plot samples versus allelic representation
plot(nSamples ~ Total, data=N4800_DF[which(N4800_DF$marker == "DNA"),], 
     ylab='Number of Samples', xlab='Allelic Representation (%)', 
     pch=16, col=alpha(plotColors[[1]], 0.01), 
     xlim=c(30,100.5), ylim=c(0,4900), cex.axis=1.2, 
     main='DNA 4800 MSSE Prediction Model: Cubic')
# Add fitted prediction line and intervals
lines(N4800_DNA_X, N4800_DNA_cubPredict[,'fit'], col = 'orange', lwd = 2)
lines(N4800_DNA_X, N4800_DNA_cubPredict[,'lwr'], col = 'red', lty = 2)
lines(N4800_DNA_X, N4800_DNA_cubPredict[,'upr'], col = 'red', lty = 2)
# Add legend
legend(x=30, y=4900, legend = c("Fitted line", "95% PI"), 
       col = c("blue","red"), lwd = 2, lty = c(1,2), bty = "n")
# EXPONENTIAL
# Plot samples versus allelic representation
plot(nSamples ~ Total, data=N4800_DF[which(N4800_DF$marker == "DNA"),], 
     ylab='Number of Samples', xlab='Allelic Representation (%)', 
     pch=16, col=alpha(plotColors[[1]], 0.01), 
     xlim=c(30,100.5), ylim=c(0,4900), cex.axis=1.2, 
     main='DNA 4800 MSSE Prediction Model: Exponential')
# Add fitted prediction line and intervals
lines(N4800_DNA_X, N4800_DNA_expPredict[,'fit'], col = 'blue', lwd = 2)
lines(N4800_DNA_X, N4800_DNA_expPredict[,'lwr'], col = 'purple', lty = 2)
lines(N4800_DNA_X, N4800_DNA_expPredict[,'upr'], col = 'purple', lty = 2)
# Add legend
legend(x=30, y=4900, legend = c("Fitted line", "95% PI"), 
       col = c("blue","red"), lwd = 2, lty = c(1,2), bty = "n")
