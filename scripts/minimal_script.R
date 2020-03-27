# SIMVEE minimal script
#install.packages("timereg")
library(dplyr)
library(foreign)
library(survival)
library(splines)
library(timereg)
library(ggplot2)

### Source functions from other files
source('R/simvee.R')
source('R/readParams.R')
source('R/ve_methods.R')

### Read parameters from input files
#   you can specify the folder and file names of the input file within the ""
#   if no path is specified a window will pop up and allow you to choose a file from your computer
params <- readParams()

### run simulation
#   there is an optional path argument for run_simvee(params, path = )
#   if no path is specified, it will default to current working directory
outcomes_dat <- run_simvee(params)

### read in outcomes file
#   you can specify the file name/path of the output file inside ""
# outcomes_dat <- read.csv("Outcomes_Test_05.csv")

# add FARI indicator variable
outcomes_dat <- outcomes_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                                        DINF_new = ifelse(DINF == 0, 999, DINF))

### apply VE estimation methods 
reject_h0_durham <- 0
reject_h0_tian <- 0
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
    reject_h0_durham <- reject_h0_durham + ifelse(flu_zph$table[1,3] < 0.05, 1, 0)
  # calculate VE
  # the nsmo argument indicates the number of time points to calculate VE at 
    temp <- durham_ve(flu_zph, nsmo = 20, var = "V") %>% mutate(Sim = i, Method = "Durham")
    if (i > 1){
    ve_est <- bind_rows(ve_est,temp)
    } else {ve_est <- temp}
    
    ######################################
    ### method from Tian et al. 2005 ###
    ######################################
    # calculate VE
    # n_timepoint_breaks argument specifies the number of time points to calculate VE for
    temp2 <- tian_ve(outcomes_dat1, n_timepoint_breaks = 20) 
    temp2a <- temp2$output %>% mutate(Sim = i, Method = "Tian")
    ve_est <- bind_rows(ve_est,temp2a)
    # proportion of sims where null hypothesis is rejected
    reject_h0_tian <- reject_h0_tian + temp2$reject_h0
    
    #######################################
    ### method from Ainslie et al. 2017 ###
    #######################################
    temp3 <- ainslie_ve(outcomes_dat1, n_days = params$ND)
    temp3a <- temp3 %>% mutate(Sim = i, Method = "Ainslie")
    ve_est <- bind_rows(ve_est,temp3a)
}

### proportion of sims where H0 rejected
prop_reject_h0 <- tibble(method = c('Durham','Tian'),
                         proportion = c(reject_h0_durham/params$sim, reject_h0_tian/params$sim))
write.csv(prop_reject_h0,file = paste0("Reject_H0_Prop_",params$title,".csv"))
### mean VE from simulations at each timepoint
mean_ve <- ve_est %>% group_by(Method, time) %>% summarise_at("ve", c(mean, sd)) %>% rename(ve_mean = fn1, ve_sd = fn2)
write.csv(mean_ve,file = paste0("Mean_VE_Estimates_",params$title,".csv"))


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





