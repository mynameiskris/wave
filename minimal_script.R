# SIMVEE minimal script

library(dplyr)
# Source functions from other files
source('simvee.R')
source("readParams.R")

set.seed(1234)
# Read parameters from input files
params <- readParams("SimVEE_input__Test_03.csv")

# run simulation
mysims <- run_sim(params)
