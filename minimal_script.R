# SIMVEE minimal script

library(dplyr)
# Source functions from other files
source('simvee.R')
source("readParams.R")

set.seed(1234)
# Read parameters from input files
params <- readParams("SimVEE_input.csv")

# run simulation
mysims <- run_sim(params)

# write output to file
write_output(params, mysims, path = "output/")

# post-processing (still needs to be coded...)
