parTab <- data.frame(values=c(0.4,0.1,0, params$ND, params$NJ, params$NDJ, 1, 4),
                     names=c("alpha","theta_0","phi", "n_days", "n_periods", "n_days_period",
                             "latent_period", "infectious_period"),
                     fixed=c(0.0001,0.0001,0.0001,rep(1,5)),
                     lower_bound=c(rep(0,8)),
                     upper_bound=c(1,1,2, rep(1000, 5)),
                     steps=c(0.01,0.01, 0.01, rep(0,5)),
                     stringsAsFactors=FALSE)


## You MUST specify the arguments to this function as parTab, data then PRIOR_FUNC.
## Use the `...` to pass additional arguments through.
my_creation_function <- function(parTab, data, PRIOR_FUNC, ...){
  ##############################
  ## This is where you would manipulate all
  ## of your model arguments. For example,
  ## unpackaging the stuff from `...`,
  ## or transforming model parameters
  ##############################

  ## Somewhat redundant example
  parameter_names <- parTab$names

  N <- length(unique(data$ID))
  n_days <- parTab[parTab$names == "n_days",]$values
  n_periods <- parTab[parTab$names == "n_periods",]$values
  n_days_period <- parTab[parTab$names == "n_days_period",]$values
  latent_period <- parTab[parTab$names == "latent_period",]$values
  infectious_period <- parTab[parTab$names == "infectious_period",]$values

  prev <- numeric(n_days)

  for (d in 1:n_days){
    # calculate which days individuals got infected to be infectious on day d
    possible_day_of_infection <- (d  - latent_period - infectious_period):(d - latent_period)
    prev[d] <- length(which(data$DINF_new %in% possible_day_of_infection))/N
  }
  prev <- ifelse(prev == 0, 0.001, prev)

  x <- list(n = N,
            n_days = n_days,
            n_days_period = n_days_period,
            prev = prev,
            dinf = data$DINF_new,
            v = data$V
  )

  ##############################
  ## This is where you would put your own model code
  ##############################
  likelihood_func <- function(pars){
      #names(pars) <- parameter_names
      alpha = pars[1]
      theta_0 = pars[2]
      phi = pars[3]

      pi_0u <- pi_0v <- pi_1u <- pi_1v <- c(1,rep(0,x$n_days-1))
      psi_0u <- psi_0v <- psi_1u <- psi_1v <- pi_0u

      #initialise period
      period <- 1
      period_start_days <- seq(1, x$n_days, by = x$n_days_period)
      # loop over days
      for (d in 2:x$n_days){
        if(d %in% period_start_days){period <- period + 1}
        #print(period)
        lambda <- phi - 1
        theta_d <- theta_0 + lambda * period
        alpha_d <- alpha * x$prev[d]
        # conditional probabilities: pi_ju & pi_jv, where
        #   j = infection status,
        #   u = unvaccinated,
        #   v = vaccinated
        pi_0u[d] <- ifelse(1 - alpha_d < 0, 0.0001, 1 - alpha_d)
        pi_0v[d] <- ifelse(1 - alpha_d * theta_d < 0, 0.0001, 1 - alpha_d * theta_d)
        pi_1u[d] <- alpha_d
        pi_1v[d] <- ifelse(alpha_d * theta_d > 1, 1, alpha_d * theta_d)
        # unconditional probabilities: psi_ju & psi_jv, where
        #   j = infection status,
        #   u = unvaccinated,
        #   v = vaccinated
        if (d == 1){
          psi_0u[d] <- pi_0u[d]
          psi_0v[d] <- pi_0v[d]
          psi_1u[d] <- pi_1u[d]
          psi_1v[d] <- pi_1v[d]
        } else {
          psi_0u[d] <- pi_0u[d] * psi_0u[d-1]          # not infected, unvaccinated
          psi_0v[d] <- pi_0v[d] * psi_0v[d-1]          # not infected, vaccinated
          psi_1u[d] <- pi_1u[d] * psi_0u[d-1]          # infected, unvaccinated
          psi_1v[d] <- pi_1v[d] * psi_0v[d-1]          # infected, vaccinated
        }
      }
      # personal contribution to the likelihood ---------------------------------
      lik <- numeric(x$n)

      for (i in 1:x$n){
        if( x$dinf[i] == 999 ){
          if( x$v[i] == 0 ){ lik[i] <- psi_0u[x$n_days] }
          if( x$v[i] == 1 ){ lik[i] <- psi_0v[x$n_days] }
        }
        if ( x$dinf[i] != 999 ){
          if( x$v[i] == 0 ){ lik[i] <- psi_1u[x$dinf[i]] }
          if( x$v[i] == 1 ){ lik[i] <- psi_1v[x$dinf[i]] }
        }
      }

    ## Note the use of closures to capture the `data` argument
    loglik <- sum(log(lik)) #sum(dnorm(data, mu, sd, TRUE))
    if(!is.null(PRIOR_FUNC)) loglik <- loglik + PRIOR_FUNC(pars)
    loglik
  }
  return(likelihood_func)
}

## Note that we've given it a prior function too. lazymcmc will
## combine this function pointer with the likelihood function.
# my_prior <- function(pars){
#   a <- dnorm(pars["mu"],5,10,1)
#   b <- dnorm(pars["sd"],2,10,1)
#   return(a + b)
# }
## I will use the default option of a uniform prior, so PRIOR_FUNC = NULL
## To test that it's working
posterior <- my_creation_function(parTab, data = outcomes_dat, PRIOR_FUNC = NULL)
print(posterior(parTab$values[1:3]))

real_pars <- c(parTab$values[1],parTab$values[2],parTab$values[3])

## Likelihood of true pars
print(likelihood_func(x,real_pars))
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
startTab$values[1:3] <- c(0.3, 0.2, 1.3)

## You could have done something like this:
## startTab$values <- runif(nrow(startTab), startTab$lower_bound, startTab$upper_bound)

output <- run_MCMC(parTab=startTab, data=data, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=my_creation_function, mvrPars=NULL,
                   PRIOR_FUNC = NULL, OPT_TUNING=0.2)

# plot results (exclude adaptive period)
chain <- read.csv(output$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))

