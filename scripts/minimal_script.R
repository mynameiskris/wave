# SIMVEE minimal script
# update: 11 May, 2020
# time: 18:02
# install.packages("timereg")
# install.packages("DEoptim")
library(dplyr)
library(foreign)
library(survival)
library(splines)
library(timereg)
library(ggplot2)
# library(DEoptim)

### Source functions from other files
source('R/simvee.R')
source('R/readParams.R')
source('R/ve_methods.R')

### Read parameters from input files
#   you can specify the folder and file names of the input file within the ""
#   if no path is specified a window will pop up and allow you to choose a file from your computer
params <- readParams("input/Input_ban_403.csv")

### run simulation
#   there is an optional path argument for run_simvee(params, path = )
#   if no path is specified, it will default to current working directory
outcomes_dat <- run_simvee(params)

### read in outcomes file
# you can specify the file name/path of the output file inside ""
# outcomes_dat <- read.csv(file.choose())

# add FARI indicator variable
outcomes_dat <- outcomes_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                                        DINF_new = ifelse(DINF == 0, 999, DINF))

# apply VE estimation methods
ve_estimates <- estimate_ve(dat = outcomes_dat, params)

# print proportion of null hypotheses rejected
ve_estimates$prop_reject_h0

# print mean mle parameter estimates
ve_estimates$mean_mle_params




