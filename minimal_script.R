# SIMVEE minimal script

library(dplyr)
# Source functions from other files
source('simvee.R')
source("readParams.R")

set.seed(1234)
# Read parameters from input files
params <- readParams("input/SimVEE_input__Test_03.csv")

# run simulation
#   there is an optional path argument for run_sim(params, path = )
#   if no path is specified, it will default to current working directory
mysims <- run_sim(params)

# manipulate data to determine VE over time
flu_dat <- mysims %>% select(-c(Sim, ID))
# apply method from Durham et al. 1988
source('ve_methods.R')
# fit ordinary Cox propotional hazards model
flu_coxmod <- coxph(Surv(DINF) ~ X + V, data=flu_dat)
# test the proportional hazards assumption and compute the Schoenfeld residuals ($y)
flu_tdhaz <- cox.zph(flu_coxmod, transform = "identity")
# calculate VE
ve_est <- durham_ve(flu_tdhaz, nsmo = 10, var = "V")


