# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE MSAT AND SNP SIMULATED DATASETS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in genind files generated from previously run fastSimcoal simulations,
# then subsets each file to specify a group of "ex situ" (garden) individuals. The ex situ
# representation (i.e. how well do these garden individuals represent total allelic diversity) is calculated.

# Then, the remaining individuals ("wild") are resampled iteratively, and the allelic diversity
# of sample subsets (in comparison to the whole of wild allelic diversity) is calculated, then plotted

# In order to function iteratively over large objects (i.e. lists of genind objects), the steps
# in this script use many apply family functions

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
# Parallelism: specify number of cores to use
num_cores <- detectCores() - 4
# Specify number of resampling replicates, to use for all scenarios below
num_reps <- 5
# Pick plot colors (for all plots!)
plotColors <- c("red","red4","darkorange3","coral","purple")
# Flags for processing different datasets (nInd=1200 and nInd=4800; different n_to_drop flags)
nInd_1200_Flag <- TRUE
nInd_4800_Flag <- TRUE
N0_Flag <- TRUE
N1_Flag <- TRUE
N2_Flag <- TRUE

# %%% ORIGINAL SIMULATIONS: NIND 1200 %%% ----
# Based on flag value, analyze NIND 1200 dataset
if(nInd_1200_Flag==TRUE){
  print("%%% ANALYZING NIND=1200 DATASET %%%")
  # %%% Read in simulations and process results ----
  # Run the simulations
  # source("RScripts/GenerateFSCparams_N1200.R")
  # Alternatively, source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N1200_marker/data.MSAT/"))
  readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N1200_marker/data.DNA/"))
  
  # Combine all scenarios for each marker into a list of genind lists
  MSAT_geninds <- list(MSAT_01pop_migLow.genList, MSAT_01pop_migHigh.genList, MSAT_04pop_migLow.genList,
                       MSAT_04pop_migHigh.genList, MSAT_16pop_migLow.genList, MSAT_16pop_migHigh.genList)
  DNA_geninds <- list(DNA_01pop_migLow.genList, DNA_01pop_migHigh.genList, DNA_04pop_migLow.genList,
                      DNA_04pop_migHigh.genList, DNA_16pop_migLow.genList, DNA_16pop_migHigh.genList)
  
  # %%% Assign garden samples ----
  # Specify the proportion of total individuals assigned to gardens
  gardenRate <- 0.05
  # Assign individuals to garden populations. Since we're dealing with a list of lists, we use rapply
  # (We could also use a nested lapply: genind_list <- lapply (genind_list, lapply, function, args))
  MSAT_geninds <- rapply(MSAT_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
  DNA_geninds <- rapply(DNA_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
  
  if(N0_Flag==TRUE){
    # %%% Summarize simulations ----
    # MSAT
    summarize_simulations(MSAT_geninds)
    # # DNA
    summarize_simulations(DNA_geninds)
    
    # %%% Resampling ----
    # MSAT 
    # Declare filepath to save resampling array to
    N1200_MSAT_N0_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/MSAT_N1200_marker/data.MSAT/N0/MSAT_N1200_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N1200_MSAT_N0_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores)
    saveRDS(N1200_MSAT_N0_resamplingArrays, file=N1200_MSAT_resampArr_N0_filepath)
    # DNA
    # Declare filepath to save resampling array to
    N1200_DNA_N0_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/DNA_N1200_marker/data.DNA/N0/DNA_N1200_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N1200_DNA_N0_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores)
    saveRDS(N1200_DNA_N0_resamplingArrays, file=N1200_DNA_N0_resampArr_filepath)
    
    # %%% Summarize resampling results ----
    # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
    # MSAT
    N1200_MSAT_N0_min95_Means <- rapply(N1200_MSAT_N0_resamplingArrays, resample_min95_mean, how = "list")
    N1200_MSAT_N0_meanValues <- rapply(N1200_MSAT_N0_resamplingArrays, resample_meanValues, how = "list")
    # DNA
    N1200_DNA_min95Values <- rapply(N1200_DNA_N0_resamplingArrays, resample_min95_mean, how = "list")
    N1200_DNA_meanValues <- rapply(N1200_DNA_N0_resamplingArrays, resample_meanValues, how = "list")
    
    # %%% Plotting ----
    # Specify the directories to save plots to
    N1200_MSAT_N0_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N1200/MSAT/N0/"
    N1200_DNA_N0_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N1200/DNA/N0/"
    # Plotting commands nested in invisible function, to prevent text from being printed
    # MSAT
    invisible(rapply(N1200_MSAT_N0_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_MSAT_N0_plotDir))
    # DNA
    invisible(rapply(N1200_DNA_N0_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_DNA_N0_plotDir))
  }
  if(n1_Flag==TRUE){
    # %%% Summarize simulations ----
    # MSAT
    summarize_simulations(MSAT_geninds, n_to_drop=1)
    # # DNA
    summarize_simulations(DNA_geninds, n_to_drop=1)
    
    # %%% Resampling ----
    # MSAT 
    # Declare filepath to save resampling array to
    N1200_MSAT_N1_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/MSAT_N1200_marker/data.MSAT/N1/MSAT_N1200_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N1200_MSAT_N1_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=1)
    saveRDS(N1200_MSAT_N1_resamplingArrays, file=N1200_MSAT_resampArr_N1_filepath)
    # DNA
    # Declare filepath to save resampling array to
    N1200_DNA_N1_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/DNA_N1200_marker/data.DNA/N1/DNA_N1200_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N1200_DNA_N1_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=1)
    saveRDS(N1200_DNA_N1_resamplingArrays, file=N1200_DNA_N1_resampArr_filepath)
    
    # %%% Summarize resampling results ----
    # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
    # MSAT
    N1200_MSAT_N1_min95_Means <- rapply(N1200_MSAT_N1_resamplingArrays, resample_min95_mean, how = "list")
    N1200_MSAT_N1_meanValues <- rapply(N1200_MSAT_N1_resamplingArrays, resample_meanValues, how = "list")
    # DNA
    N1200_DNA_min95Values <- rapply(N1200_DNA_N1_resamplingArrays, resample_min95_mean, how = "list")
    N1200_DNA_meanValues <- rapply(N1200_DNA_N1_resamplingArrays, resample_meanValues, how = "list")
    
    # %%% Plotting ----
    # Specify the directories to save plots to
    N1200_MSAT_N1_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N1200/MSAT/N1/"
    N1200_DNA_N1_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N1200/DNA/N1/"
    # Plotting commands nested in invisible function, to prevent text from being printed
    # MSAT
    invisible(rapply(N1200_MSAT_N1_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_MSAT_N1_plotDir))
    # DNA
    invisible(rapply(N1200_DNA_N1_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_DNA_N1_plotDir))
  }
  if(n2_Flag==TRUE){
    # %%% Summarize simulations ----
    # MSAT
    summarize_simulations(MSAT_geninds, n_to_drop=2)
    # # DNA
    summarize_simulations(DNA_geninds, n_to_drop=2)
    
    # %%% Resampling ----
    # MSAT 
    # Declare filepath to save resampling array to
    N1200_MSAT_N2_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/MSAT_N1200_marker/data.MSAT/N1/MSAT_N1200_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N1200_MSAT_N2_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=2)
    saveRDS(N1200_MSAT_N2_resamplingArrays, file=N1200_MSAT_resampArr_N2_filepath)
    # DNA
    # Declare filepath to save resampling array to
    N1200_DNA_N2_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/DNA_N1200_marker/data.DNA/N1/DNA_N1200_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N1200_DNA_N2_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=2)
    saveRDS(N1200_DNA_N2_resamplingArrays, file=N1200_DNA_N2_resampArr_filepath)
    
    # %%% Summarize resampling results ----
    # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
    # MSAT
    N1200_MSAT_N2_min95_Means <- rapply(N1200_MSAT_N2_resamplingArrays, resample_min95_mean, how = "list")
    N1200_MSAT_N2_meanValues <- rapply(N1200_MSAT_N2_resamplingArrays, resample_meanValues, how = "list")
    # DNA
    N1200_DNA_min95Values <- rapply(N1200_DNA_N2_resamplingArrays, resample_min95_mean, how = "list")
    N1200_DNA_meanValues <- rapply(N1200_DNA_N2_resamplingArrays, resample_meanValues, how = "list")
    
    # %%% Plotting ----
    # Specify the directories to save plots to
    N1200_MSAT_N2_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N1200/MSAT/N1/"
    N1200_DNA_N2_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N1200/DNA/N1/"
    # Plotting commands nested in invisible function, to prevent text from being printed
    # MSAT
    invisible(rapply(N1200_MSAT_N2_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_MSAT_N2_plotDir))
    # DNA
    invisible(rapply(N1200_DNA_N2_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_DNA_N2_plotDir))
  }
}

# %%% ORIGINAL SIMULATIONS: NIND 4800 %%% ----
# Based on flag value, analyze NIND 4800 dataset
if(nInd_4800_Flag==TRUE){
  print("%%% ANALYZING NIND=4800 DATASET %%%")
  # %%% Read in simulations and process results ----
  # Run the simulations
  # source("RScripts/GenerateFSCparams_N4800.R")
  # Alternatively, source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N4800_marker/data.MSAT/"))
  readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N4800_marker/data.DNA/"))
  
  # Combine all scenarios for each marker into a list of genind lists
  MSAT_geninds <- list(MSAT_01pop_migLow.genList, MSAT_01pop_migHigh.genList, MSAT_04pop_migLow.genList,
                       MSAT_04pop_migHigh.genList, MSAT_16pop_migLow.genList, MSAT_16pop_migHigh.genList)
  DNA_geninds <- list(DNA_01pop_migLow.genList, DNA_01pop_migHigh.genList, DNA_04pop_migLow.genList,
                      DNA_04pop_migHigh.genList, DNA_16pop_migLow.genList, DNA_16pop_migHigh.genList)
  
  # %%% Assign garden samples ----
  # Specify the proportion of total individuals assigned to gardens
  gardenRate <- 0.05
  # Assign individuals to garden populations. Since we're dealing with a list of lists, we use rapply
  # (We could also use a nested lapply: genind_list <- lapply (genind_list, lapply, function, args))
  MSAT_geninds <- rapply(MSAT_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
  DNA_geninds <- rapply(DNA_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
  
  if(N0_Flag==TRUE){
    # %%% Summarize simulations ----
    # MSAT
    summarize_simulations(MSAT_geninds)
    # # DNA
    summarize_simulations(DNA_geninds)
    
    # %%% Resampling ----
    # MSAT 
    # Declare filepath to save resampling array to
    N4800_MSAT_N0_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/MSAT_N4800_marker/data.MSAT/N0/MSAT_N4800_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N4800_MSAT_N0_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores)
    saveRDS(N4800_MSAT_N0_resamplingArrays, file=N4800_MSAT_resampArr_N0_filepath)
    # DNA
    # Declare filepath to save resampling array to
    N4800_DNA_N0_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/DNA_N4800_marker/data.DNA/N0/DNA_N4800_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N4800_DNA_N0_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores)
    saveRDS(N4800_DNA_N0_resamplingArrays, file=N4800_DNA_N0_resampArr_filepath)
    
    # %%% Summarize resampling results ----
    # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
    # MSAT
    N4800_MSAT_N0_min95_Means <- rapply(N4800_MSAT_N0_resamplingArrays, resample_min95_mean, how = "list")
    N4800_MSAT_N0_meanValues <- rapply(N4800_MSAT_N0_resamplingArrays, resample_meanValues, how = "list")
    # DNA
    N4800_DNA_min95Values <- rapply(N4800_DNA_N0_resamplingArrays, resample_min95_mean, how = "list")
    N4800_DNA_meanValues <- rapply(N4800_DNA_N0_resamplingArrays, resample_meanValues, how = "list")
    
    # %%% Plotting ----
    # Specify the directories to save plots to
    N4800_MSAT_N0_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N4800/MSAT/N0/"
    N4800_DNA_N0_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N4800/DNA/N0/"
    # Plotting commands nested in invisible function, to prevent text from being printed
    # MSAT
    invisible(rapply(N4800_MSAT_N0_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N4800_MSAT_N0_plotDir))
    # DNA
    invisible(rapply(N4800_DNA_N0_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N4800_DNA_N0_plotDir))
  }
  if(n1_Flag==TRUE){
    # %%% Summarize simulations ----
    # MSAT
    summarize_simulations(MSAT_geninds, n_to_drop=1)
    # # DNA
    summarize_simulations(DNA_geninds, n_to_drop=1)
    
    # %%% Resampling ----
    # MSAT 
    # Declare filepath to save resampling array to
    N4800_MSAT_N1_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/MSAT_N4800_marker/data.MSAT/N1/MSAT_N4800_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N4800_MSAT_N1_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=1)
    saveRDS(N4800_MSAT_N1_resamplingArrays, file=N4800_MSAT_resampArr_N1_filepath)
    # DNA
    # Declare filepath to save resampling array to
    N4800_DNA_N1_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/DNA_N4800_marker/data.DNA/N1/DNA_N4800_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N4800_DNA_N1_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=1)
    saveRDS(N4800_DNA_N1_resamplingArrays, file=N4800_DNA_N1_resampArr_filepath)
    
    # %%% Summarize resampling results ----
    # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
    # MSAT
    N4800_MSAT_N1_min95_Means <- rapply(N4800_MSAT_N1_resamplingArrays, resample_min95_mean, how = "list")
    N4800_MSAT_N1_meanValues <- rapply(N4800_MSAT_N1_resamplingArrays, resample_meanValues, how = "list")
    # DNA
    N4800_DNA_min95Values <- rapply(N4800_DNA_N1_resamplingArrays, resample_min95_mean, how = "list")
    N4800_DNA_meanValues <- rapply(N4800_DNA_N1_resamplingArrays, resample_meanValues, how = "list")
    
    # %%% Plotting ----
    # Specify the directories to save plots to
    N4800_MSAT_N1_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N4800/MSAT/N1/"
    N4800_DNA_N1_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N4800/DNA/N1/"
    # Plotting commands nested in invisible function, to prevent text from being printed
    # MSAT
    invisible(rapply(N4800_MSAT_N1_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N4800_MSAT_N1_plotDir))
    # DNA
    invisible(rapply(N4800_DNA_N1_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N4800_DNA_N1_plotDir))
  }
  if(n2_Flag==TRUE){
    # %%% Summarize simulations ----
    # MSAT
    summarize_simulations(MSAT_geninds, n_to_drop=2)
    # # DNA
    summarize_simulations(DNA_geninds, n_to_drop=2)
    
    # %%% Resampling ----
    # MSAT 
    # Declare filepath to save resampling array to
    N4800_MSAT_N2_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/MSAT_N4800_marker/data.MSAT/N1/MSAT_N4800_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N4800_MSAT_N2_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=2)
    saveRDS(N4800_MSAT_N2_resamplingArrays, file=N4800_MSAT_resampArr_N2_filepath)
    # DNA
    # Declare filepath to save resampling array to
    N4800_DNA_N2_resampArr_filepath <- 
      paste0(sim.wd, "SimulationOutputs/DNA_N4800_marker/data.DNA/N1/DNA_N4800_resampArr.Rdata")
    # Run resampling in parallel, and save the resampling array result to specified location
    N4800_DNA_N2_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=2)
    saveRDS(N4800_DNA_N2_resamplingArrays, file=N4800_DNA_N2_resampArr_filepath)
    
    # %%% Summarize resampling results ----
    # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
    # MSAT
    N4800_MSAT_N2_min95_Means <- rapply(N4800_MSAT_N2_resamplingArrays, resample_min95_mean, how = "list")
    N4800_MSAT_N2_meanValues <- rapply(N4800_MSAT_N2_resamplingArrays, resample_meanValues, how = "list")
    # DNA
    N4800_DNA_min95Values <- rapply(N4800_DNA_N2_resamplingArrays, resample_min95_mean, how = "list")
    N4800_DNA_meanValues <- rapply(N4800_DNA_N2_resamplingArrays, resample_meanValues, how = "list")
    
    # %%% Plotting ----
    # Specify the directories to save plots to
    N4800_MSAT_N2_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N4800/MSAT/N1/"
    N4800_DNA_N2_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/nToDropTests_052023/N4800/DNA/N1/"
    # Plotting commands nested in invisible function, to prevent text from being printed
    # MSAT
    invisible(rapply(N4800_MSAT_N2_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N4800_MSAT_N2_plotDir))
    # DNA
    invisible(rapply(N4800_DNA_N2_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N4800_DNA_N2_plotDir))
  }
}
