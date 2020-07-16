# ----------------------------------------------------------
#' Estimate VE using a maximum likelihood (ML) method with MCMC
#'
#' This method is similar to that of Ainslie et al. 2017 (SIM). The ML method is based on calculating the
#' contribution to the likelihood function of each study participant. See XX for more details.
#' @param data data set
#' @param params input parameters
#' @param my_prior prior function for MCMC
#' @param latent_period length of latent period
#' @param infectious_period length of infectious period
#' @return list with a tibble of VE estimates for each period (the estimate for each period is the average VE over the days
#' within each period) and the maximu likelihood estimates of the parameters.
#' @keywords wave
#' @import dplyr
#' @import tidyr
#' @import lazymcmc
#' @export
ml_ve2 <- function(data, params, my_prior = NULL, latent_period = 1, infectious_period = 4,
                   file_name = "test"){

phi <- params$theta_d[1] - params$theta_d[2] + 1

parTab <- data.frame(values=c(params$alpha_0, params$theta_d[1], phi,
                              params$ND, params$NJ, params$NDJ, latent_period,
                              infectious_period),
                     names=c("alpha","theta_0","phi", "n_days", "n_periods", "n_days_period",
                             "latent_period", "infectious_period"),
                     fixed=c(0,0,0,rep(1,5)),
                     steps=c(rep(0.01,8)),
                     lower_bound=c(rep(0.0001,3),rep(0,5)),
                     upper_bound=c(1,1,2, rep(1000, 5)),
                     stringsAsFactors=FALSE)

## Quick test
#real_pars <- c(parTab$values[1],parTab$values[2],parTab$values[3])
## Likelihood of true pars
#print(likelihood_func(pars = real_pars))

## Update the proposal step size every 1000 iterations (opt_freq) for the first 5000 iterations
## (adaptive_period) to aim for an acceptance rate of 0.44 (popt). After the adaptive period,
## run for a further 10000 (iterations) steps. Save every nth rows, where n is "thin" (ie. 1 here).
## Write to disk every 100 iterations (post thinning). Note that the adaptive period is also saved
## to disk
mcmcPars <- c("iterations"=10000,"popt"=0.44,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=5000,"save_block"=100)

## The MCMC code uses the parameter table. Here, we should specify some random starting
## points in the "values" column.
startTab <- parTab
startTab$values[1:3] <- c(0.5, 0.2, 1.3)

## You could have done something like this:
## startTab$values <- runif(nrow(startTab), startTab$lower_bound, startTab$upper_bound)

output <- run_MCMC(parTab = startTab, data = data, mcmcPars = mcmcPars,
                   filename = file_name,
                   CREATE_POSTERIOR_FUNC = my_creation_function, mvrPars = NULL,
                   PRIOR_FUNC = my_prior, OPT_TUNING = 0.2)

# plot results (exclude adaptive period)
chain <- read.csv(output$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))

## Determine mle
chain1 <- chain[chain$sampno > mcmcPars["adaptive_period"],]
max_loglik <- which(chain1$lnlike == max(chain1$lnlike))
mles <- unique(chain1[max_loglik, c("alpha", "theta_0", "phi")])

periods <- rep(1:params$NJ, each = params$NDJ)
ve_dat <- tibble(day = 1:params$ND, period = periods, ve = 1 - (mles$theta_0 + ((mles$phi - 1) * .data$day))) %>%
  select(-.data$day) %>%
  group_by(.data$period) %>%
  summarise_all(.funs = mean) %>%
  mutate(ve = ifelse(ve > 1, 1, ve))

# output
rtn <- list(param_est = mles, ve_dat = ve_dat)

return(rtn)
}
