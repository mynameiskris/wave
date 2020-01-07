# development script

library(dplyr)
# Source functions from other files
source('simvee.R')
source("readParams.R")

set.seed(1234)
# Read parameters from input files
params <- readParams("input/SimVEE_input__Test_03.csv")

# run simulation
sim_results <- run_sim(params, path = "output/")
