# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% TEST CASE: 14 IUCN RED LIST THREATENED US OAKS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in the resampling arrays generated in Rosenberger et al. 2022, which examines
# simulated datasets tailored to model 14 species of IUCN Red List Threatened oaks in the U.S.
# It then calculates the minimum sample size estimates (MSSEs) and the prediction intervals
# surrounding those estimates. 

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)

# Set working directory to the GitHub repo
sim.wd <- '/home/akoontz/Shared/SSRvSNP_Sim/Code/'
setwd(sim.wd)
# Read in relevant functions
source('RScripts/0_functions.R')
# Load .Rdata file containing resampling arrays for 14 oak species
# This command will create and object in the environment named 'final_quercus_results', which is
# an array in which the rows are number of samples (1 to 500), the columns are resampling replicates,
# and the slices (3rd dimension) are the different oak species.
load('SimulationOutputs/Rosenberger2022_IUCN14Oaks.Rdata')

# Predeclare a vector capturing allelic representation values
quercus14_genDiv <- vector(length = (nrow(final_quercus_results)))
# Predeclare a matrix storing the results of the predict function
quercus14_predict_Matrix <- matrix(ncol = 3, nrow = 14)
colnames(quercus14_predict_Matrix) <- c("fit", "lwr", "upr")
rownames(quercus14_predict_Matrix) <- c("QUAC","QUAR","QUAU", "QUBO","QUCA","QUCE","QUEN","QUGE","QUGR","QUHA","QUHI","QUOG","QUPA", "QUTO")
# Loop through species (slices in array)
for (i in 1:(dim(final_quercus_results)[3])){
  # Capture Total resampling values
  quercus14_genDiv <- vector(length = (nrow(final_quercus_results[,,i])) * (dim(final_quercus_results)[2]))
  quercus14_genDiv <- c(final_quercus_results[,,i])
  # Specify sample numbers column
  quercus14_sampleNumbers <- 2:(nrow(final_quercus_results[,,i])+1)
  quercus14_sampleNumbers <- rep(quercus14_sampleNumbers, dim(final_quercus_results[,,i])[2])
  # Create data.frame from resampling array values
  quercus14_DF <- data.frame(sampleNumbers=quercus14_sampleNumbers, totalValues=quercus14_genDiv)
  # Build linear models, then calculate the 95% MSSE using predict
  quercus14_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = quercus14_DF)
  quercus14_newData <- data.frame(totalValues=0.95)
  quercus14_95MSSEprediction <- predict(quercus14_Model, quercus14_newData, interval = "prediction")
  # Round the MSSE fit value to the next highest value; round PI values to 1 decimal place
  quercus14_95MSSEprediction[[1]] <- ceiling(quercus14_95MSSEprediction[[1]])
  quercus14_95MSSEprediction[2:3] <- round(quercus14_95MSSEprediction[2:3], digits = 1) 
  # Store the predict results to a matrix
  quercus14_predict_Matrix[i,] <- quercus14_95MSSEprediction
  # Calculate the prediction interval width (upper PI - lower PI), and store the result
  piWidth <- vector()
  for (j in 1:(nrow(quercus14_predict_Matrix))){
    piWidth[j] <- quercus14_predict_Matrix[j,3] - quercus14_predict_Matrix[j,2]
  } 
}
# Create and print summary matrix, with results from all species
q14Matrix <- cbind(quercus14_predict_Matrix, piWidth) ; print(q14Matrix)
