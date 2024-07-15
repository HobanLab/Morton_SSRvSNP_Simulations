# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% 03/26/2024 Simulated MSAT Loci bootstrapping %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script generates a list "QUAC_array_list" of 5 arrays. Each array stores 
# representation values across replicates from each scenario based on a number of loci
# (5-25 loci are bootstrapped at intervals of 5)

library(adegenet)
sim.wd <- 'C:/Users/gsalas/Documents/resampling_CIs/Code/'
setwd(sim.wd)
# ---- FUNCTIONS ----
# Building resampling array ----
gm_resamp_array_function <- function(insert_genind, num_loci, num_reps){
  # Create an empty array named 'resamp_category5loc' to store results
  resamp_category <- array(dim = c(nrow(insert_genind@tab),4,num_reps))
  # loop 25 times for sets (5 loci in one set) of randomly selected loci
  for (i in 1:num_reps) {
    # Randomly sample an amount of loci from the genind object based on user input
    samp_loc <- sample(locNames(insert_genind), size = num_loci, replace = FALSE)
    # Subset the genind object to include only the columns corresponding to the sampled loci
    gm.Wild.genind <- insert_genind[, loc = samp_loc]
    # declare objects
    # access the genind matrix that shows the type of alleles and quantity present among wild individuals
    wildSamp <- gm.Wild.genind@tab
    # calculate the sum of each column in the matrix, ignoring NA values.
    # identify the indices where the sum is not equal to zero, this indicates columns with 
    # variation in allele counts. Subset the original matrix by selecting only the columns identified
    # in the previous step, removing columns with no variation in allele counts
    wildSamp <- wildSamp[, which(colSums(wildSamp, na.rm = TRUE) != 0)]
    # calculating a vector of allele frequencies (sum of allele counts divided by number of haplotypes, or individuals * 2)
    wildComplete <- colSums(wildSamp, na.rm = TRUE) / (nrow(wildSamp) * 2)
    # Subset 'wildComplete' to include only the alleles with non-zero frequency
    wildSubset <- wildComplete[wildComplete > 0]
    # initialize vectors to store results 
    total <- vector(length = nrow(wildSamp))
    common <- vector(length = nrow(wildSamp))
    lowfreq <- vector(length = nrow(wildSamp))
    rare <- vector(length = nrow(wildSamp))
    # Loop for each sample size, ranging from 1 tto the number of loci in the wild population
    for (j in 1:nrow(wildSamp)) {
      # randomly select a subset of rows from the matrix, without replacement
      samp <- sample(nrow(wildSamp), size = j, replace = FALSE)
      # subset the original matrix, to inlcude only the rows randomly selected
      samp <- wildSamp[samp,]
      # Measure the proportion of allelic representation in that sample
      # Check if its the first iteration of the inner loop
      if (j == 1) {
        # Identify the names of alleles in the sample that are also present in 'wildSubset'. 
        # Calculate the proportion of alleles in the sample that are also present in 'wildSubset'.
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        # Identify the names of the alleles in the sample that are also present in the 'wildSubset' for alleles with frequency >0.1
        # Calculate the proportion of common alleles in the sample compared to the 'wildSubset' for alleles with frequency >0.1
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        # Identify the names of alleles in the sample that are also present in the 'wildSubset' for alleles with frequency >0.01 and frequency <0.1
        # Calculate the proportion of low frequency alleles in the sample compared to the 'wildSubset' for allels with frequency >0.01 and frequency <0.1
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        # Identify the names of alleles in the sample that are also present in the 'wildSubset' for allleles with frequency <0.02.
        # Calcualte the proportion of rare alleles in the sample compared to the 'wildSubset' for alleles with frequency <0.01.
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        # Calculate the proportion of alleles in the sample that are also present in the 'wildSubset' for the case where j is not 1
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        # Calculate the proportion of alleles of common alleles in the sample compared toe the 'wildSubset' where j is not 1
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        # Calculate the proportion of alleles of low frequency alleles in the sample compared to the 'wildSubset' for the case where j is not 1
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        # Calculate the proportion of rare alleles in the sample compared ot the 'wildSubset' for the case where j is not 1
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      # Combine the results (proportions) for each sample size into a matrix named 'categorymat_10loc'.
      categorymat <- cbind(total, common, lowfreq, rare)
      # extract the column names (categories) from the matrix
      category <- colnames(categorymat)
      # set row and column names for 'resamp_category10loc'
      dimnames(resamp_category) <- list(paste0("sample ", 1:nrow(categorymat)), category, paste0("replicate ", 1:num_reps))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category[, , i] <- categorymat
    }
  }
  return(resamp_category)
}

# Linear model ----
# pass all arrays to a dataframe using the function. define the function that takes an input 'data_array'
analyze_resampling_array <- function(data_array, allele_cat) {
  # linear model of resampling array. Extract the column named 'total' from the array.
  # concatenate the extracted column values into the vector.
  totalsVector <- c(data_array[,allele_cat,]) 
  
  # Specify sample numbers column.
  # Create a vector of sample numbers from 1 to the number of rows in the 'total' column 
  gm_sampleNumbers <- 1:(nrow(data_array[,allele_cat,]))
  # Repeat the sample numbers vector for the number of replicates
  gm_sampleNumbers <- rep(gm_sampleNumbers, dim(data_array)[[3]])
  
  # Create data frame from resampling array values
  gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=totalsVector)
  
  # Build and analyze linear models
  gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
  # Create a new data fram 'gm_newData with a single column 'totalValues' containing the value 0.95
  gm_newData <- data.frame(totalValues=0.95)
  # Use the linear model to predict the response for the new data frame. Specify 'interval = prediction to obtain a prediction interval
  gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
  
  # Pass the gm_95MSSEprediction to the object storing our results 
  # Store the predicted values and predictino interval in the object named 'result' 
  result <- gm_95MSSEprediction
  # Calculate the width of the prediction interval by substracting the lower limit from the upper limit
  piWidth <- gm_95MSSEprediction[3] - gm_95MSSEprediction[2]
  # Return a list containing the predicted values and the width of the prediction interval
  return(list(result = result, piWidth = piWidth))
}

# Filling in matrix ----
# Iterate through the arrays and store results in the matrix
# Initiate loop that iterates over the indices of 'array_list'
build_matrix_func <- function(array_list, input_matrix){
  for (i in 1:length(array_list)) {
    # Store results and piWidth values in the ith row of the matrix
    input_matrix[i, ] <- c(array_list[[i]]$result, 
                           array_list[[i]]$piWidth)
  }
  return(input_matrix)
}
# ---- Loci bootstrapping representation value replicates across scenarios ----
# %%% MSAT N1200 %%% ----
# declare all objects for the loop
# this array stores a resampling array from a simulated genind object of a scenario, of which there are 5 simulated genind objects.
# this array stores all the replicates from each scenario at a specific loci level
combinedArray <- array(dim = c(1200, 4, 6*5*1))
MSAT_levels <- c(5, 10, 15, 20, 25)
MSAT_N1200_arrayList <- vector("list", length = length(MSAT_levels))
# reading in the names of each of the scenarios to be processed
MSATscenarios <- list.files(path = 'Datasets/Simulated/MSAT_N1200/', pattern = "genind.MSAT_", full.names = TRUE)
# number of replicates we want the resampling array to have
numReps <- 1
# this for loop creates a list that contains all the simulation replicates for every scenario
scenariorepsList <- list()
for(i in 1:length(MSATscenarios)){
  # browser()
  currentScenario <- readRDS(MSATscenarios[[i]])
  scenariorepsList <- c(currentScenario,scenariorepsList)
}


# Print starting time
startTime <- Sys.time()
print(paste0('%%% RUNTIME START: ', startTime))
for (i in 1:length(MSAT_levels)) {
  for (j in 1:length(scenariorepsList)){
    # browser()
    combinedArray[,,j] <- gm_resamp_array_function(scenariorepsList[[j]],MSAT_levels[i], numReps)
  }
  MSAT_N1200_arrayList[[i]] <- combinedArray
  combinedArray <- array(dim = c(1200,4,6*5*1))
}
# Print ending time and total runtime
endTime <- Sys.time()
print(paste0('%%% RUNTIME END: ', endTime))
cat(paste0('\n', '%%% TOTAL RUNTIME: ', endTime-startTime))

saveRDS(MSAT_N1200_arrayList, "MSAT_N1200_arrayList.Rdata")

# %%% MSAT N4800 %%% ----
# declare all objects for the loop
# this array stores the arrays of all the scenarios, of which there are 5 simulation replicates, and 1 replicate per simulation replicate
combinedArray <- array(dim = c(4800, 4, 6*5*1))
MSAT_levels <- c(5, 10, 15, 20, 25)
MSAT_N4800_arrayList <- vector("list", length = length(MSAT_levels))
# reading in the names of each of the scenarios to be processed
MSATscenarios <- list.files(path = 'Datasets/Simulated/MSAT_N4800/', pattern = "genind.MSAT_", full.names = TRUE)
# number of times you do the resampling of each resampling replicate
numReps <- 1

scenariorepsList <- list()
for(i in 1:length(MSATscenarios)){
  # browser()
  currentScenario <- readRDS(MSATscenarios[[i]])
  scenariorepsList <- c(currentScenario,scenariorepsList)
}

# Print starting time
startTime <- Sys.time()
print(paste0('%%% RUNTIME START: ', startTime))
for (i in 1:length(MSAT_levels)) {
  for (j in 1:length(scenariorepsList)){
    # browser()
    combinedArray[,,j] <- gm_resamp_array_function(scenariorepsList[[j]],MSAT_levels[i], numReps)
  }
  MSAT_N4800_arrayList[[i]] <- combinedArray
  combinedArray <- array(dim = c(4800,4,6*5*1))
}
# Print ending time and total runtime
endTime <- Sys.time()
print(paste0('%%% RUNTIME END: ', endTime))
cat(paste0('\n', '%%% TOTAL RUNTIME: ', endTime-startTime))

saveRDS(MSAT_N4800_arrayList, "MSAT_N4800_arrayList.Rdata")

# %%% DNA N1200 %%% ----
combinedArray <- array(dim = c(1200, 4, 6*5*1))
DNA_levels <- c(100,250,500,750,1000)
DNA_N1200_arrayList <- vector("list", length = length(DNA_levels))
QUAC_array_list = list(length(MSAT_levels))
QUAC_array_list <- vector("list", length = length(MSAT_levels))
# reading in the names of each of the scenarios to be processed
DNAscenarios <- list.files(path = 'Datasets/Simulated/DNA_N1200/', pattern = "genind.DNA_", full.names = TRUE)
# a list that stores all the simulation scenarios after being added to the environment

numReps <- 1

scenariorepsList <- list()
for(i in 1:length(DNAscenarios)){
  # browser()
  currentScenario <- readRDS(DNAscenarios[[i]])
  scenariorepsList <- c(currentScenario,scenariorepsList)
}



# Print starting time
startTime <- Sys.time()
print(paste0('%%% RUNTIME START: ', startTime))
for (i in 1:length(DNA_levels)) {
  for (j in 1:length(scenariorepsList)){
    # browser()
    combinedArray[,,j] <- gm_resamp_array_function(scenariorepsList[[j]],DNA_levels[i], numReps)
  }
  DNA_N1200_arrayList[[i]] <- combinedArray
  combinedArray <- array(dim = c(1200,4,6*5*1))
}
# Print ending time and total runtime
endTime <- Sys.time()
print(paste0('%%% RUNTIME END: ', endTime))
cat(paste0('\n', '%%% TOTAL RUNTIME: ', endTime-startTime))
saveRDS(DNA_N1200_arrayList, "DNA_N1200_arrayList.Rdata")

# %%% DNA N4800 %%% ----
combinedArray <- array(dim = c(4800, 4, 6*5*1))
DNA_levels <- c(100,250,500,750,1000)
DNA_N4800_arrayList <- vector("list", length = length(DNA_levels))
# reading in the names of each of the scenarios to be processed
DNAscenarios <- list.files(path = 'Datasets/Simulated/DNA_N4800/', pattern = "genind.DNA_", full.names = TRUE)
# a list that stores all the simulation scenarios after being added to the environment

numReps <- 1

scenariorepsList <- list()
for(i in 1:length(DNAscenarios)){
  # browser()
  currentScenario <- readRDS(DNAscenarios[[i]])
  scenariorepsList <- c(currentScenario,scenariorepsList)
}



# Print starting time
startTime <- Sys.time()
print(paste0('%%% RUNTIME START: ', startTime))
for (i in 1:length(DNA_levels)) {
  for (j in 1:length(scenariorepsList)){
    # browser()
    combinedArray[,,j] <- gm_resamp_array_function(scenariorepsList[[j]], DNA_levels[i], numReps)
  }
  DNA_N4800_arrayList[[i]] <- combinedArray
  combinedArray <- array(dim = c(4800,4,6*5*1))
}
# Print ending time and total runtime
endTime <- Sys.time()
print(paste0('%%% RUNTIME END: ', endTime))
cat(paste0('\n', '%%% TOTAL RUNTIME: ', endTime-startTime))

saveRDS(DNA_N4800_arrayList, "DNA_N4800_arrayList.Rdata")

# ---- ALLELE CATEGORIES ----
sim.wd <- "C:/Users/gsalas/Documents/Morton_SSRvSNP_Simulations/SimulationOutputs"
setwd(sim.wd)
MSAT_N1200 <- readRDS(paste0(sim.wd,"/MSAT_N1200_marker/lociBootstrapping/MSAT_N1200_arrayList.Rdata"))
MSAT_N4800 <- readRDS(paste0(sim.wd, "/MSAT_N4800_marker/lociBootstrapping/MSAT_N4800_arrayList.Rdata"))
DNA_N1200 <- readRDS(paste0(sim.wd,"/DNA_N1200_marker/lociBootstrapping/DNA_N1200_arrayList.Rdata"))
DNA_N4800 <- readRDS(paste0(sim.wd, "/DNA_N4800_marker/lociBootstrapping/DNA_N4800_arrayList.Rdata"))
MSAT_levels <- c("5 loci", "10 loci", "15 loci", "20 loci", "25 loci")
DNA_levels <- c("100 loci","250 loci","500 loci","750 loci","1000 loci")
predict_outputs <- c("MSSE", "lower", "upper", "piWidth")
# empty list with length of 5 that will contain the outputs of the predict analysis
MSAT_N1200rareCatpredictResults <-vector("list", length(MSAT_levels))
MSAT_N1200lowfreqpredictResults <-vector("list", length(MSAT_levels))
MSAT_N1200commonpredictResults <-vector("list", length(MSAT_levels))
MSAT_N1200_predict_results <- vector("list", length(MSAT_levels))
MSAT_N4800rareCatpredictResults <-vector("list", length(MSAT_levels))
MSAT_N4800lowfreqpredictResults <-vector("list", length(MSAT_levels))
MSAT_N4800commonpredictResults <-vector("list", length(MSAT_levels))
MSAT_N4800_predict_results <- vector("list", length(MSAT_levels))
DNA_N1200rareCatpredictResults <-vector("list", length(DNA_levels))
DNA_N1200lowfreqpredictResults <-vector("list", length(DNA_levels))
DNA_N1200commonpredictResults <-vector("list", length(DNA_levels))
DNA_N1200_predict_results <- vector("list", length(DNA_levels))
DNA_N4800rareCatpredictResults <-vector("list", length(DNA_levels))
DNA_N4800lowfreqpredictResults <-vector("list", length(DNA_levels))
DNA_N4800commonpredictResults <-vector("list", length(DNA_levels))
DNA_N4800_predict_results <- vector("list", length(DNA_levels))
# %%% Create an empty matrix to store the pi results ----
# Set column names for 'results_Matrix'
# Set row names for 'resutls_Matrix'
MSAT_N1200rareCat_matrix <- matrix(nrow = length(MSAT_levels), ncol = length(predict_outputs))
rownames(MSAT_N1200rareCat_matrix) <- MSAT_levels
colnames(MSAT_N1200rareCat_matrix) <- predict_outputs

MSAT_N1200lowfreqCat_matrix <- matrix(nrow = length(MSAT_levels), ncol = length(predict_outputs))
colnames(MSAT_N1200lowfreqCat_matrix) <- predict_outputs
rownames(MSAT_N1200lowfreqCat_matrix) <- MSAT_levels

MSAT_N1200commonCat_matrix <- matrix(nrow = length(MSAT_levels), ncol = length(predict_outputs))
colnames(MSAT_N1200commonCat_matrix) <- predict_outputs
rownames(MSAT_N1200commonCat_matrix) <- MSAT_levels

MSAT_N1200_matrix <- matrix(nrow = length(MSAT_levels), ncol = length(predict_outputs))
colnames(MSAT_N1200_matrix) <- predict_outputs
rownames(MSAT_N1200_matrix) <- MSAT_levels

MSAT_N4800rareCat_matrix <- matrix(nrow = length(MSAT_levels), ncol = length(predict_outputs))
colnames(MSAT_N4800rareCat_matrix) <- predict_outputs
rownames(MSAT_N4800rareCat_matrix) <- MSAT_levels

MSAT_N4800lowfreqCat_matrix <- matrix(nrow = length(MSAT_levels), ncol = length(predict_outputs))
colnames(MSAT_N4800lowfreqCat_matrix) <- predict_outputs
rownames(MSAT_N4800lowfreqCat_matrix) <- MSAT_levels

MSAT_N4800commonCat_matrix <- matrix(nrow = length(MSAT_levels), ncol = length(predict_outputs))
colnames(MSAT_N4800commonCat_matrix) <- predict_outputs
rownames(MSAT_N4800commonCat_matrix) <- MSAT_levels

MSAT_N4800_matrix <- matrix(nrow = length(MSAT_levels), ncol = length(predict_outputs))
colnames(MSAT_N4800_matrix) <- predict_outputs
rownames(MSAT_N4800_matrix) <- MSAT_levels

DNA_N1200rareCat_matrix <- matrix(nrow = length(DNA_levels), ncol = length(predict_outputs))
colnames(DNA_N1200rareCat_matrix) <- predict_outputs
rownames(DNA_N1200rareCat_matrix) <- DNA_levels

DNA_N1200lowfreqCat_matrix <- matrix(nrow = length(DNA_levels), ncol = length(predict_outputs))
colnames(DNA_N1200lowfreqCat_matrix) <- predict_outputs
rownames(DNA_N1200lowfreqCat_matrix) <- DNA_levels

DNA_N1200commonCat_matrix <- matrix(nrow = length(DNA_levels), ncol = length(predict_outputs))
colnames(DNA_N1200commonCat_matrix) <- predict_outputs
rownames(DNA_N1200commonCat_matrix) <- DNA_levels

DNA_N1200_matrix <- matrix(nrow = length(DNA_levels), ncol = length(predict_outputs))
colnames(DNA_N1200_matrix) <- predict_outputs
rownames(DNA_N1200_matrix) <- DNA_levels

DNA_N4800rareCat_matrix <- matrix(nrow = length(DNA_levels), ncol = length(predict_outputs))
colnames(DNA_N4800rareCat_matrix) <- predict_outputs
rownames(DNA_N4800rareCat_matrix) <- DNA_levels

DNA_N4800lowfreqCat_matrix <- matrix(nrow = length(DNA_levels), ncol = length(predict_outputs))
colnames(DNA_N4800lowfreqCat_matrix) <- predict_outputs
rownames(DNA_N4800lowfreqCat_matrix) <- DNA_levels

DNA_N4800commonCat_matrix <- matrix(nrow = length(DNA_levels), ncol = length(predict_outputs))
colnames(DNA_N4800commonCat_matrix) <- predict_outputs
rownames(DNA_N4800commonCat_matrix) <- DNA_levels

DNA_N4800_matrix <- matrix(nrow = length(DNA_levels), ncol = length(predict_outputs))
colnames(DNA_N4800_matrix) <- predict_outputs
rownames(DNA_N4800_matrix) <- DNA_levels

# %%% Prediction outputs ---- 
# prints out the prediction width results and the fit for each array with the same number of loci
# %%% MSAT N1200 %%% ----
for (i in 1:length(MSAT_levels)) {
  MSAT_N1200rareCatpredictResults[[i]] <- analyze_resampling_array(MSAT_N1200[[i]],4)
}

for (i in 1:length(MSAT_levels)) {
  MSAT_N1200lowfreqpredictResults[[i]] <- analyze_resampling_array(MSAT_N1200[[i]],3)
}

for (i in 1:length(MSAT_levels)) {
  MSAT_N1200commonpredictResults[[i]] <- analyze_resampling_array(MSAT_N1200[[i]],2)
}

for (i in 1:length(MSAT_levels)) {
  MSAT_N1200_predict_results[[i]] <- analyze_resampling_array(MSAT_N1200[[i]],1)
}

# %%% MSAT N4800 %%% ----
for (i in 1:5) {
  MSAT_N4800rareCatpredictResults[[i]] <- analyze_resampling_array(MSAT_N4800[[i]],4)
}

for (i in 1:5) {
  MSAT_N4800lowfreqpredictResults[[i]] <- analyze_resampling_array(MSAT_N4800[[i]],3)
}

for (i in 1:5) {
  MSAT_N4800commonpredictResults[[i]] <- analyze_resampling_array(MSAT_N4800[[i]],2)
}

for (i in 1:length(MSAT_levels)) {
  MSAT_N4800_predict_results[[i]] <- analyze_resampling_array(MSAT_N4800[[i]],1)
}

# %%% DNA 1200 %%% ----
for (i in 1:length(DNA_levels)) {
  DNA_N1200rareCatpredictResults[[i]] <- analyze_resampling_array(DNA_N1200[[i]],4)
}

for (i in 1:length(DNA_levels)) {
  DNA_N1200lowfreqpredictResults[[i]] <- analyze_resampling_array(DNA_N1200[[i]],3)
}

for (i in 1:length(DNA_levels)) {
  DNA_N1200commonpredictResults[[i]] <- analyze_resampling_array(DNA_N1200[[i]],2)
}

for (i in 1:length(DNA_levels)) {
  DNA_N1200_predict_results[[i]] <- analyze_resampling_array(DNA_N1200[[i]],1)
}

# %%% DNA N4800 %%% ----

for (i in 1:length(DNA_levels)) {
  DNA_N4800rareCatpredictResults[[i]] <- analyze_resampling_array(DNA_N4800[[i]],4)
}

for (i in 1:length(DNA_levels)) {
  DNA_N4800lowfreqpredictResults[[i]] <- analyze_resampling_array(DNA_N4800[[i]],3)
}

for (i in 1:length(DNA_levels)) {
  DNA_N4800commonpredictResults[[i]] <- analyze_resampling_array(DNA_N4800[[i]],2)
}

for (i in 1:length(DNA_levels)) {
  DNA_N4800_predict_results[[i]] <- analyze_resampling_array(DNA_N4800[[i]],1)
}

# %%% Insert list and empty matrix to the build matrix function & rounding ----
MSAT_N1200rareCat_matrix <- round(build_matrix_func(MSAT_N1200rareCatpredictResults, MSAT_N1200rareCat_matrix),digits = 1)
MSAT_N1200lowfreqCat_matrix <- round(build_matrix_func(MSAT_N1200lowfreqpredictResults, MSAT_N1200lowfreqCat_matrix),digits = 1)
MSAT_N1200commonCat_matrix <- round(build_matrix_func(MSAT_N1200commonpredictResults, MSAT_N1200commonCat_matrix),digits=1)
MSAT_N1200_matrix <- round(build_matrix_func(MSAT_N1200_predict_results, MSAT_N1200_matrix),digits = 1)

MSAT_N4800rareCat_matrix <- round(build_matrix_func(MSAT_N4800rareCatpredictResults, MSAT_N4800rareCat_matrix),digits = 1)
MSAT_N4800lowfreqCat_matrix <- round(build_matrix_func(MSAT_N4800lowfreqpredictResults, MSAT_N4800lowfreqCat_matrix),digits=1)
MSAT_N4800commonCat_matrix <- round(build_matrix_func(MSAT_N4800commonpredictResults, MSAT_N4800commonCat_matrix),digits = 1)
MSAT_N4800_matrix <- round(build_matrix_func(MSAT_N4800_predict_results, MSAT_N4800_matrix),digits=1)

DNA_N1200rareCat_matrix <- round(build_matrix_func(DNA_N1200rareCatpredictResults, DNA_N1200rareCat_matrix),digits=1)
DNA_N1200lowfreqCat_matrix <- round(build_matrix_func(DNA_N1200lowfreqpredictResults, DNA_N1200lowfreqCat_matrix),digits = 1)
DNA_N1200commonCat_matrix <- round(build_matrix_func(DNA_N1200commonpredictResults, DNA_N1200commonCat_matrix), digits = 1)
DNA_N1200_matrix <- round(build_matrix_func(DNA_N1200_predict_results, DNA_N1200_matrix),digits = 1)

DNA_N4800rareCat_matrix <- round(build_matrix_func(DNA_N4800rareCatpredictResults, DNA_N4800rareCat_matrix),digits = 1)
DNA_N4800lowfreqCat_matrix <- round(build_matrix_func(DNA_N4800lowfreqpredictResults, DNA_N4800lowfreqCat_matrix), digits = 1)
DNA_N4800commonCat_matrix <- round(build_matrix_func(DNA_N4800commonpredictResults, DNA_N4800commonCat_matrix), digits = 1)
DNA_N4800_matrix <- round(build_matrix_func(DNA_N4800_predict_results, DNA_N4800_matrix),digits = 1)

# %%% write csv to outputs folder ----
write.csv(MSAT_N1200rareCat_matrix, file = paste0(sim.wd, "/MSAT_N1200_marker/lociBootstrapping/MSAT_N1200_Rare.csv"))
write.csv(MSAT_N1200lowfreqCat_matrix, file = paste0(sim.wd, "/MSAT_N1200_marker/lociBootstrapping/MSAT_N1200_LowFreq.csv"))
write.csv(MSAT_N1200commonCat_matrix, file = paste0(sim.wd, "/MSAT_N1200_marker/lociBootstrapping/MSAT_N1200_Common.csv"))
write.csv(MSAT_N1200_matrix, file = paste0(sim.wd, "/MSAT_N1200_marker/lociBootstrapping/MSAT_N1200_Total.csv"))

write.csv(MSAT_N4800rareCat_matrix, file = paste0(sim.wd, "/MSAT_N4800_marker/lociBootstrapping/MSAT_N4800_Rare.csv"))
write.csv(MSAT_N4800lowfreqCat_matrix, file = paste0(sim.wd, "/MSAT_N4800_marker/lociBootstrapping/MSAT_N4800_LowFreq.csv"))
write.csv(MSAT_N4800commonCat_matrix, file = paste0(sim.wd, "/MSAT_N4800_marker/lociBootstrapping/MSAT_N4800_Common.csv"))
write.csv(MSAT_N4800_matrix, file = paste0(sim.wd, "/MSAT_N4800_marker/lociBootstrapping/MSAT_N4800_Total.csv"))

write.csv(DNA_N1200rareCat_matrix, file = paste0(sim.wd, "/DNA_N1200_marker/lociBootstrapping/DNA_N1200_Rare.csv"))
write.csv(DNA_N1200lowfreqCat_matrix, file = paste0(sim.wd, "/DNA_N1200_marker/lociBootstrapping/DNA_N1200_LowFreq.csv"))
write.csv(DNA_N1200commonCat_matrix, file = paste0(sim.wd, "/DNA_N1200_marker/lociBootstrapping/DNA_N1200_Common.csv"))
write.csv(DNA_N1200_matrix, file = paste0(sim.wd, "/DNA_N1200_marker/lociBootstrapping/DNA_N1200_Total.csv"))

write.csv(DNA_N4800rareCat_matrix, file = paste0(sim.wd, "/DNA_N4800_marker/lociBootstrapping/DNA_N4800_Rare.csv"))
write.csv(DNA_N4800lowfreqCat_matrix, file = paste0(sim.wd, "/DNA_N4800_marker/lociBootstrapping/DNA_N4800_LowFreq.csv"))
write.csv(DNA_N4800commonCat_matrix, file = paste0(sim.wd, "/DNA_N4800_marker/lociBootstrapping/DNA_N4800_Common.csv"))
write.csv(DNA_N4800_matrix, file = paste0(sim.wd, "/DNA_N4800_marker/lociBootstrapping/DNA_N4800_Total.csv"))