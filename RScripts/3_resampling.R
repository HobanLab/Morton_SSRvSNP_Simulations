# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE MSAT AND SNP SIMULATED DATASETS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in genind files generated from previously run fastSimcoal simulations,
# then assigns all individuals in that genind file to a 'wild' population, in order to include
# them in resampling analyses. 

# Then, individuals are resampled iteratively, and the allelic diversity of sample subsets 
# (in comparison to the whole of wild allelic diversity) is calculated, then plotted.

# In order to function iteratively over large objects (i.e. lists of genind objects), the steps
# in this script use many apply family functions. To use this script multiple times without overwriting
# data, we frame code sections in if() statements using flags, to modularly resample different datasets.

library(strataG)
library(adegenet)
library(stringr)
library(parallel)
library(RColorBrewer)
library(scales)

# %%% FUNCTIONS AND VARIABLES %%% ----
# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
sim.wd <- '/home/akoontz/Shared/SSRvSNP_Sim/Code/'
setwd(sim.wd)
# Read in relevant functions
source('RScripts/0_functions.R')
# Parallelism: specify number of cores to use
num_cores <- detectCores() - 4
# Specify number of resampling replicates, to use for all scenarios below
num_reps <- 5
# Pick plot colors (for all plots!)
plotColors <- c('red','red4','darkorange3','coral','purple')
# Flags for processing different datasets (marker=msat, dna; nInd=1200, 4800; different n_to_drop flags; DNA low mutation)
MSAT_Flag <- TRUE
DNA_Flag <- FALSE
nInd_1200_Flag <- TRUE
nInd_4800_Flag <- TRUE
N0_Flag <- TRUE
N1_Flag <- FALSE
N2_Flag <- FALSE
dnaLowMut_Flag <- FALSE

# Testing loci bootstrapping


# %%% ORIGINAL SIMULATIONS: NIND 1200 %%% ----
# Based on flag value, analyze NIND 1200 dataset
if(nInd_1200_Flag==TRUE){
  print('%%% ANALYZING NIND=1200 DATASET %%%')
  # %%% Read in simulations and process results ----
  # Source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_MSAT(paste0(sim.wd,'SimulationOutputs/MSAT_N1200_marker/data.MSAT/'))
  readGeninds_DNA(paste0(sim.wd,'SimulationOutputs/DNA_N1200_marker/data.DNA/'))
  
  # Combine all scenarios for each marker into a list of genind lists
  MSAT_geninds <- list(MSAT_01pop_migLow.genList, MSAT_01pop_migHigh.genList, MSAT_04pop_migLow.genList,
                       MSAT_04pop_migHigh.genList, MSAT_16pop_migLow.genList, MSAT_16pop_migHigh.genList)
  DNA_geninds <- list(DNA_01pop_migLow.genList, DNA_01pop_migHigh.genList, DNA_04pop_migLow.genList,
                      DNA_04pop_migHigh.genList, DNA_16pop_migLow.genList, DNA_16pop_migHigh.genList)
  
  # %%% Assign samples to 'wild' population ----
  MSAT_geninds <- rapply(MSAT_geninds, assignSamplePopulations, how = 'list')
  DNA_geninds <- rapply(DNA_geninds, assignSamplePopulations, how = 'list')
  
  if(N0_Flag==TRUE){
    print('%%% N TO DROP VALUE: 0')
    if(MSAT_Flag==TRUE){
      # %%% MSAT ----
      print('--- MSAT ---')
      # Summarize simulations
      summarize_simulations(MSAT_geninds)
      # Declare filepath to save resampling array to
      N1200_MSAT_N0_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/MSAT_N1200_marker/data.MSAT/N0/MSAT_N1200_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N1200_MSAT_N0_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores)
      saveRDS(N1200_MSAT_N0_resamplingArrays, file=N1200_MSAT_N0_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N1200_MSAT_N0_min95_Means <- rapply(N1200_MSAT_N0_resamplingArrays, resample_min95_mean, how = 'list')
      N1200_MSAT_N0_meanValues <- rapply(N1200_MSAT_N0_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N1200_MSAT_N0_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N1200/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N1200_MSAT_N0_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_MSAT_N0_plotDir))
    }
    if(DNA_Flag==TRUE){
      # %%% DNA ----
      print('--- DNA ---')
      # Summarize simulations
      summarize_simulations(DNA_geninds)
      # Declare filepath to save resampling array to
      N1200_DNA_N0_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/DNA_N1200_marker/data.DNA/N0/DNA_N1200_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N1200_DNA_N0_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores)
      saveRDS(N1200_DNA_N0_resamplingArrays, file=N1200_DNA_N0_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N1200_DNA_N0_min95_Means <- rapply(N1200_DNA_N0_resamplingArrays, resample_min95_mean, how = 'list')
      N1200_DNA_N0_meanValues <- rapply(N1200_DNA_N0_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N1200_DNA_N0_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N1200/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N1200_DNA_N0_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_DNA_N0_plotDir))
    }
  }
  if(N1_Flag==TRUE){
    print('%%% N TO DROP VALUE: 1')
    if(MSAT_Flag==TRUE){
      # %%% MSAT ----
      print('--- MSAT ---')
      # Summarize simulations
      summarize_simulations(MSAT_geninds, n_to_drop=1)
      # Declare filepath to save resampling array to
      N1200_MSAT_N1_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/MSAT_N1200_marker/data.MSAT/N1/MSAT_N1200_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N1200_MSAT_N1_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=1)
      saveRDS(N1200_MSAT_N1_resamplingArrays, file=N1200_MSAT_N1_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N1200_MSAT_N1_min95_Means <- rapply(N1200_MSAT_N1_resamplingArrays, resample_min95_mean, how = 'list')
      N1200_MSAT_N1_meanValues <- rapply(N1200_MSAT_N1_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N1200_MSAT_N1_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N1200/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N1200_MSAT_N1_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_MSAT_N1_plotDir))
    }
    if(DNA_Flag==TRUE){
      # %%% DNA ----
      print('--- DNA ---')
      # Summarize simulations
      summarize_simulations(DNA_geninds, n_to_drop=1)
      # Declare filepath to save resampling array to
      N1200_DNA_N1_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/DNA_N1200_marker/data.DNA/N1/DNA_N1200_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N1200_DNA_N1_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=1)
      saveRDS(N1200_DNA_N1_resamplingArrays, file=N1200_DNA_N1_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N1200_DNA_N1_min95_Means <- rapply(N1200_DNA_N1_resamplingArrays, resample_min95_mean, how = 'list')
      N1200_DNA_N1_meanValues <- rapply(N1200_DNA_N1_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N1200_DNA_N1_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N1200/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N1200_DNA_N1_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_DNA_N1_plotDir))
    }
  }
  if(N2_Flag==TRUE){
    print('%%% N TO DROP VALUE: 2')
    if(MSAT_Flag==TRUE){
      # %%% MSAT ----
      print('--- MSAT ---')
      # Summarize simulations
      summarize_simulations(MSAT_geninds, n_to_drop=2)
      # Declare filepath to save resampling array to
      N1200_MSAT_N2_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/MSAT_N1200_marker/data.MSAT/N2/MSAT_N1200_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N1200_MSAT_N2_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=2)
      saveRDS(N1200_MSAT_N2_resamplingArrays, file=N1200_MSAT_N2_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N1200_MSAT_N2_min95_Means <- rapply(N1200_MSAT_N2_resamplingArrays, resample_min95_mean, how = 'list')
      N1200_MSAT_N2_meanValues <- rapply(N1200_MSAT_N2_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N1200_MSAT_N2_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N1200/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N1200_MSAT_N2_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_MSAT_N2_plotDir))
    }
    if(DNA_Flag==TRUE){
      # %%% DNA ----
      print('--- DNA ---')
      # Summarize simulations
      summarize_simulations(DNA_geninds, n_to_drop=2)
      # Declare filepath to save resampling array to
      N1200_DNA_N2_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/DNA_N1200_marker/data.DNA/N2/DNA_N1200_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N1200_DNA_N2_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=2)
      saveRDS(N1200_DNA_N2_resamplingArrays, file=N1200_DNA_N2_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N1200_DNA_N2_min95_Means <- rapply(N1200_DNA_N2_resamplingArrays, resample_min95_mean, how = 'list')
      N1200_DNA_N2_meanValues <- rapply(N1200_DNA_N2_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N1200_DNA_N2_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N1200/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N1200_DNA_N2_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_DNA_N2_plotDir))
    }
  }
  print('%%% N1200 ANALYSES COMPLETE! %%%')
}

# %%% ORIGINAL SIMULATIONS: NIND 4800 %%% ----
# Based on flag value, analyze NIND 4800 dataset
if(nInd_4800_Flag==TRUE){
  print('%%% ANALYZING NIND=4800 DATASET %%%')
  # %%% Read in simulations and process results ----
  # Source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_MSAT(paste0(sim.wd,'SimulationOutputs/MSAT_N4800_marker/data.MSAT/'))
  readGeninds_DNA(paste0(sim.wd,'SimulationOutputs/DNA_N4800_marker/data.DNA/'))
  
  # Combine all scenarios for each marker into a list of genind lists
  MSAT_geninds <- list(MSAT_01pop_migLow.genList, MSAT_01pop_migHigh.genList, MSAT_04pop_migLow.genList,
                       MSAT_04pop_migHigh.genList, MSAT_16pop_migLow.genList, MSAT_16pop_migHigh.genList)
  DNA_geninds <- list(DNA_01pop_migLow.genList, DNA_01pop_migHigh.genList, DNA_04pop_migLow.genList,
                      DNA_04pop_migHigh.genList, DNA_16pop_migLow.genList, DNA_16pop_migHigh.genList)
  
  # %%% Assign samples to 'wild' population ----
  MSAT_geninds <- rapply(MSAT_geninds, assignSamplePopulations, how = 'list')
  DNA_geninds <- rapply(DNA_geninds, assignSamplePopulations, how = 'list')
  
  if(N0_Flag==TRUE){
    print('%%% N TO DROP VALUE: 0')
    if(MSAT_Flag==TRUE){
      # %%% MSAT ----
      print('--- MSAT ---')
      # Summarize simulations
      summarize_simulations(MSAT_geninds)
      # Declare filepath to save resampling array to
      N4800_MSAT_N0_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/MSAT_N4800_marker/data.MSAT/N0/MSAT_N4800_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N4800_MSAT_N0_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores)
      saveRDS(N4800_MSAT_N0_resamplingArrays, file=N4800_MSAT_N0_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N4800_MSAT_N0_min95_Means <- rapply(N4800_MSAT_N0_resamplingArrays, resample_min95_mean, how = 'list')
      N4800_MSAT_N0_meanValues <- rapply(N4800_MSAT_N0_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N4800_MSAT_N0_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N4800/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N4800_MSAT_N0_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE,
                       colors=plotColors, data.dir=N4800_MSAT_N0_plotDir))
    }
    if(DNA_Flag==TRUE){
      # %%% DNA ----
      print('--- DNA ---')
      # Summarize simulations
      summarize_simulations(DNA_geninds)
      # Declare filepath to save resampling array to
      N4800_DNA_N0_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/DNA_N4800_marker/data.DNA/N0/DNA_N4800_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N4800_DNA_N0_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores)
      saveRDS(N4800_DNA_N0_resamplingArrays, file=N4800_DNA_N0_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N4800_DNA_N0_min95_Means <- rapply(N4800_DNA_N0_resamplingArrays, resample_min95_mean, how = 'list')
      N4800_DNA_N0_meanValues <- rapply(N4800_DNA_N0_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N4800_DNA_N0_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N4800/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N4800_DNA_N0_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE,
                       colors=plotColors, data.dir=N4800_DNA_N0_plotDir))
    }
  }
  if(N1_Flag==TRUE){
    print('%%% N TO DROP VALUE: 1')
    if(MSAT_Flag==TRUE){
      # %%% MSAT ----
      print('--- MSAT ---')
      # Summarize simulations
      summarize_simulations(MSAT_geninds, n_to_drop=1)
      # Declare filepath to save resampling array to
      N4800_MSAT_N1_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/MSAT_N4800_marker/data.MSAT/N1/MSAT_N4800_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N4800_MSAT_N1_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=1)
      saveRDS(N4800_MSAT_N1_resamplingArrays, file=N4800_MSAT_N1_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N4800_MSAT_N1_min95_Means <- rapply(N4800_MSAT_N1_resamplingArrays, resample_min95_mean, how = 'list')
      N4800_MSAT_N1_meanValues <- rapply(N4800_MSAT_N1_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N4800_MSAT_N1_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N4800/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N4800_MSAT_N1_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE,
                       colors=plotColors, data.dir=N4800_MSAT_N1_plotDir))
    }
    if(DNA_Flag==TRUE){
      # %%% DNA ----
      print('--- DNA ---')
      # Summarize simulations
      summarize_simulations(DNA_geninds, n_to_drop=1)
      # Declare filepath to save resampling array to
      N4800_DNA_N1_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/DNA_N4800_marker/data.DNA/N1/DNA_N4800_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N4800_DNA_N1_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=1)
      saveRDS(N4800_DNA_N1_resamplingArrays, file=N4800_DNA_N1_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N4800_DNA_N1_min95_Means <- rapply(N4800_DNA_N1_resamplingArrays, resample_min95_mean, how = 'list')
      N4800_DNA_N1_meanValues <- rapply(N4800_DNA_N1_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N4800_DNA_N1_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N4800/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N4800_DNA_N1_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE,
                       colors=plotColors, data.dir=N4800_DNA_N1_plotDir))
    }
  }
  if(N2_Flag==TRUE){
    print('%%% N TO DROP VALUE: 2')
    if(MSAT_Flag==TRUE){
      # %%% MSAT ----
      print('--- MSAT ---')
      # Summarize simulations
      summarize_simulations(MSAT_geninds, n_to_drop=2)
      # Declare filepath to save resampling array to
      N4800_MSAT_N2_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/MSAT_N4800_marker/data.MSAT/N2/MSAT_N4800_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N4800_MSAT_N2_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=2)
      saveRDS(N4800_MSAT_N2_resamplingArrays, file=N4800_MSAT_N2_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N4800_MSAT_N2_min95_Means <- rapply(N4800_MSAT_N2_resamplingArrays, resample_min95_mean, how = 'list')
      N4800_MSAT_N2_meanValues <- rapply(N4800_MSAT_N2_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N4800_MSAT_N2_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N4800/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N4800_MSAT_N2_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE,
                       colors=plotColors, data.dir=N4800_MSAT_N2_plotDir))
    }
    if(DNA_Flag==TRUE){
      # %%% DNA ----
      print('--- DNA ---')
      # Summarize simulations
      summarize_simulations(DNA_geninds, n_to_drop=2)
      # Declare filepath to save resampling array to
      N4800_DNA_N2_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/DNA_N4800_marker/data.DNA/N2/DNA_N4800_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N4800_DNA_N2_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores, n_to_drop=2)
      saveRDS(N4800_DNA_N2_resamplingArrays, file=N4800_DNA_N2_resampArr_filepath)
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      N4800_DNA_N2_min95_Means <- rapply(N4800_DNA_N2_resamplingArrays, resample_min95_mean, how = 'list')
      N4800_DNA_N2_meanValues <- rapply(N4800_DNA_N2_resamplingArrays, resample_meanValues, how = 'list')
      # Specify the directories to save plots to
      N4800_DNA_N2_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/LociBootstrappingUpdates_20240220/N4800/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      invisible(rapply(N4800_DNA_N2_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE,
                       colors=plotColors, data.dir=N4800_DNA_N2_plotDir))
    }
  }
  print('%%% N4800 ANALYSES COMPLETE! %%%')
}

# %%% DNA LOW MUTATION SCENARIOS %%% ----
# Based on flag value, analyze DNA datasets with lower mutation rates (1e-8) 
# We don't utilize n_to_drop flags for this analysis, but we do look at different population sizes
nInd_1200_Flag <- TRUE
nInd_4800_Flag <- TRUE

if(DNA_Flag==TRUE){
  if(dnaLowMut_Flag==TRUE){
    print('%%% DNA LOW MUTATION RESAMPLING %%%')
    # %%% Read in simulations and process results ----
    if(nInd_1200_Flag==TRUE){
      print('%%% ANALYZING N1200 DATASETS')
      # Source the genind objects from previously run simulations, using readGeninds functions
      readGeninds_DNA(paste0(sim.wd,'SimulationOutputs/DNA_N1200_lowMut/data.DNA/'))
      
      # Combine all scenarios for each marker into a list of genind lists
      DNA_geninds <- list(DNA_01pop_migLow.genList, DNA_01pop_migHigh.genList, DNA_04pop_migLow.genList,
                          DNA_04pop_migHigh.genList, DNA_16pop_migLow.genList, DNA_16pop_migHigh.genList)
      
      # %%% Assign samples to 'wild' population ----
      DNA_geninds <- rapply(DNA_geninds, assignSamplePopulations, how = 'list')
      
      # %%% Resampling ----
      # DNA
      # Declare filepath to save resampling array to
      N1200_DNA_lowMut_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/DNA_N1200_lowMut/data.DNA/DNA_N1200_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N1200_DNA_lowMut_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores)
      saveRDS(N1200_DNA_lowMut_resamplingArrays, file=N1200_DNA_lowMut_resampArr_filepath)
      
      # %%% Summarize resampling results ----
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      # DNA
      N1200_DNA_min95Values <- rapply(N1200_DNA_lowMut_resamplingArrays, resample_min95_mean, how = 'list')
      N1200_DNA_meanValues <- rapply(N1200_DNA_lowMut_resamplingArrays, resample_meanValues, how = 'list')
      
      # %%% Plotting ----
      # Specify the directories to save plots to
      N1200_DNA_lowMut_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/dnaMutationRateTests_052023/N1200/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      # DNA
      invisible(rapply(N1200_DNA_lowMut_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE, 
                       colors=plotColors, data.dir=N1200_DNA_lowMut_plotDir))
    }
    if(nInd_4800_Flag==TRUE){
      print('%%% ANALYZING N4800 DATASETS')
      # Source the genind objects from previously run simulations, using readGeninds functions
      readGeninds_DNA(paste0(sim.wd,'SimulationOutputs/DNA_N4800_lowMut/data.DNA/'))
      
      # Combine all scenarios for each marker into a list of genind lists
      DNA_geninds <- list(DNA_01pop_migLow.genList, DNA_01pop_migHigh.genList, DNA_04pop_migLow.genList,
                          DNA_04pop_migHigh.genList, DNA_16pop_migLow.genList, DNA_16pop_migHigh.genList)
      
      # %%% Assign samples to 'wild' population ----
      DNA_geninds <- rapply(DNA_geninds, assignSamplePopulations, how = 'list')
      
      # %%% Resampling ----
      # DNA
      # Declare filepath to save resampling array to
      N4800_DNA_lowMut_resampArr_filepath <- 
        paste0(sim.wd, 'SimulationOutputs/DNA_N4800_lowMut/data.DNA/DNA_N4800_resampArr.Rdata')
      # Run resampling in parallel, and save the resampling array result to specified location
      N4800_DNA_lowMut_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores)
      saveRDS(N4800_DNA_lowMut_resamplingArrays, file=N4800_DNA_lowMut_resampArr_filepath)
      
      # %%% Summarize resampling results ----
      # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
      # DNA
      N4800_DNA_min95Values <- rapply(N4800_DNA_lowMut_resamplingArrays, resample_min95_mean, how = 'list')
      N4800_DNA_meanValues <- rapply(N4800_DNA_lowMut_resamplingArrays, resample_meanValues, how = 'list')
      
      # %%% Plotting ----
      # Specify the directories to save plots to
      N4800_DNA_lowMut_plotDir <- 
        '/home/akoontz/Documents/SSRvSNP/Simulations/Documentation/Images/dnaMutationRateTests_052023/N4800/'
      # Plotting commands nested in invisible function, to prevent text from being printed
      # DNA
      invisible(rapply(N4800_DNA_lowMut_resamplingArrays, resample_Plot_PNG, largePopFlag=TRUE, 
                       colors=plotColors, data.dir=N4800_DNA_lowMut_plotDir))
    }
    print('%%% DNA LOW MUTATION ANALYSES COMPLETE! %%%')
  }
}
