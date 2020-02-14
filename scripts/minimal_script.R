# SIMVEE minimal script

library(dplyr)
library(foreign)
library(survival)
library(splines)
library(ggplot2)
# Source functions from other files
source('simvee.R')
source("readParams.R")

# Read parameters from input files
params <- readParams("input/SimVEE_input.csv")

# run simulation
#   there is an optional path argument for run_sim(params, path = )
#   if no path is specified, it will default to current working directory
mysims <- run_simvee(params)

### read in outcomes file
outcomes_dat <- read.csv("Outcomes__Test_04.csv")

### apply method from Durham et al. 1988
source('ve_methods.R')
# loop through simulations
for (i in 1:max(outcomes_dat$Sim)){
  outcomes_dat1 <- outcomes_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1)) %>% filter(Sim == i)
  # fit ordinary Cox propotional hazards model
    flu_coxmod <- coxph(Surv(DINF,FARI) ~ V, data=outcomes_dat1)
  # test the proportional hazards assumption and compute the Schoenfeld residuals ($y)
    flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
  # calculate VE
    temp <- durham_ve(flu_zph, var = "V") %>% mutate(Sim = i)
    if (i > 1){
    ve_est <- bind_rows(ve_est,temp)
    } else {ve_est <- temp}
}
# subset for 1 sim to plot
ve_est1 <- ve_est %>% filter(Sim == 7)
# plot VE
y_min <- ifelse(min(ve_est1$y_lower) < 0, min(ve_est1$y_lower), 0)
y_max <- ifelse(max(ve_est1$y_upper) > 1, max(ve_est1$y_upper), 1)

p_durham <- ggplot(data = ve_est1, aes(x = pred_x, y = y_hat)) + geom_line() +
            geom_ribbon(aes(x = pred_x, ymin = y_lower, ymax = y_upper, linetype = NA), alpha=0.2) +
            xlab("Day") + ylab("VE") + scale_y_continuous(limits = c(y_min, y_max)) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))
p_durham
