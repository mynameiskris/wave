# SIMVEE minimal script
library(dplyr)
library(foreign)
library(survival)
library(splines)
library(timereg)
library(ggplot2)
# Source functions from other files
source('simvee.R')
source('readParams.R')
source('ve_methods.R')
# Read parameters from input files
# you can specify the folder and file names of the infput file within the ""
params <- readParams("input/SimVEE_input__Test_04.csv")

# run simulation
#   there is an optional path argument for run_sim(params, path = )
#   if no path is specified, it will default to current working directory
mysims <- run_simvee(params)

### read in outcomes file
# you can specify the file name of the output file inside ""
outcomes_dat <- read.csv("Outcomes__Test_03.csv")
# add FARI indicator variable
outcomes_dat <- outcomes_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                                        DINF_new = ifelse(DINF == 0, 999, DINF))

### apply VE estimation methods 

# loop through simulations and apply each method
for (i in 1:max(outcomes_dat$Sim)){
  # subset data for each simulation
  outcomes_dat1 <- outcomes_dat %>% filter(Sim == i)
  ######################################
  ### method from Durham et al. 1988 ###
  ######################################
  # fit ordinary Cox propotional hazards model
    flu_coxmod <- coxph(Surv(DINF_new,FARI) ~ V, data=outcomes_dat1)
  # test the proportional hazards assumption and compute the Schoenfeld residuals ($y)
    flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
  # calculate VE
    temp <- durham_ve(flu_zph, var = "V") %>% mutate(Sim = i, Method = "Durham")
    if (i > 1){
    ve_est <- bind_rows(ve_est,temp)
    } else {ve_est <- temp}
    
    ######################################
    ### method from Tian et al. 2005 ###
    ######################################
    # calculate VE
    temp2 <- tian_ve(outcomes_dat1) %>% mutate(Sim = i, Method = "Tian")
      ve_est <- bind_rows(ve_est,temp2)
}

# mean VE from simulations
mean_ve <- ve_est %>% group_by(Method, time) %>% summarise_at("ve", c(mean, sd)) %>% rename(ve_mean = fn1, ve_sd = fn2)
mean_ve


### ignore this for now ###

# plot of mean VE estimates over time
# ve_est1 <- ve_est %>% filter(Method == "Durham")
# # plot VE
# y_min <- ifelse(min(ve_est1$y_lower) < 0, min(ve_est1$y_lower), 0)
# y_max <- ifelse(max(ve_est1$y_upper) > 1, max(ve_est1$y_upper), 1)
# 
# p_durham <- ggplot(data = ve_est1, aes(x = pred_x, y = y_hat)) + geom_line() +
#             geom_ribbon(aes(x = pred_x, ymin = y_lower, ymax = y_upper, linetype = NA), alpha=0.2) +
#             xlab("Day") + ylab("VE") + scale_y_continuous(limits = c(y_min, y_max)) +
#             theme(panel.grid.major = element_blank(),
#                   panel.grid.minor = element_blank(),
#                   panel.background = element_blank(),
#                   axis.line = element_line(colour = "black"))
# p_durham





