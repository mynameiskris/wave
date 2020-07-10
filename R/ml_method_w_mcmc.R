# parTab <- data.frame(values=c(0.4,0.1,1, params$),
#                      names=c("alpha","theta_0","phi", "n_days", "n_periods", "n_days_period",
#                              "latent_period", "infectious_period"),
#                      fixed=c(0,0,0,rep(1,5)),
#                      lower_bound=c(rep(0,8)),
#                      upper_bound=c(1,1,2, rep(1000, 5)),
#                      steps=c(0.01,0.01, 0.01, rep(0,5)),
#                      stringsAsFactors=FALSE)
#
# ## You MUST specify the arguments to this function as parTab, data then PRIOR_FUNC.
# ## Use the `...` to pass additional arguments through.
# my_creation_function <- function(parTab, data, PRIOR_FUNC, ...){
#   ##############################
#   ## This is where you would manipulate all
#   ## of your model arguments. For example,
#   ## unpackaging the stuff from `...`,
#   ## or transforming model parameters
#   ##############################
#
#   ## Somewhat redundant example
#   parameter_names <- parTab$names
#
#   ##############################
#   ## This is where you would put your own model code
#   ##############################
# ## in my code this function is called loglik and is defined elsewhere
#   f <- function(pars){
#     # names(pars) <- parameter_names
#     # mu <- pars["mu"]
#     # sd <- pars["sd"]
#
#     ## Note the use of closures to capture the `data` argument
#     lik <- loglike #sum(dnorm(data, mu, sd, TRUE))
#     if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
#     lik
#   }
#   return(f)
# }
#
# ## Note that we've given it a prior function too. lazymcmc will
# ## combine this function pointer with the likelihood function.
# # my_prior <- function(pars){
# #   a <- dnorm(pars["mu"],5,10,1)
# #   b <- dnorm(pars["sd"],2,10,1)
# #   return(a + b)
# # }
# ## I will use the default option of a uniform prior, so PRIOR_FUNC = NULL
# ## To test that it's working
# posterior <- my_creation_function(parTab, data, my_prior)
# print(posterior(parTab$values))
